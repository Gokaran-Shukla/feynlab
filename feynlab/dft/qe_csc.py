"""Fully charge self-consistent DFT+DMFT with Quantum ESPRESSO (>= 7.0) and Wannier90.

Interface (Beck, Hampel, Parcollet, Ederer, Georges, J. Phys.: Condens. Matter 34, 235601 (2022)):
pw.x with ``dmft = .true.`` and ``dmft_prefix`` reads ``<dmft_prefix>.h5`` containing
    dft_misc_input/band_window : (1, nk, 2)          first/last band of the correlated window (1-based)
    dft_update/delta_N         : (nk, nb, nb, 2)     Delta N_ij(k) = N^DMFT_ij(k) - delta_ij f^DFT_i(k)  (re, im)
in the Bloch-band basis, per spin; QE adds Delta N to its Fermi weights, diagonalises, rotates the
wavefunctions and builds the new density (PW/src/add_dmft_occ.f90).

One CSC cycle (all steps run by `QEDMFT.cycle`):
    1. nscf     pw.x, explicit k list (nosym, noinv), dmft=.true., verbosity='high'
    2. W90      wannier90.x -pp, pw2wannier90.x, wannier90.x  ->  U(k) (and U_dis(k)), seed_hr.dat
    3. DMFT     H_W(k) = V(k)^+ diag(eps_k) V(k),  V = U_dis U  (window -> Wannier), impurity solver in feynlab
    4. update   N^DMFT(k) = T sum_n [i w + mu - eps_k - V (Sigma - Sigma_dc) V^+]^{-1} e^{i w 0+}  (window, exact tail)
                Delta N = N^DMFT - diag(f^DFT),  mu fixed so that the window charge is conserved
    5. scf      pw.x restart_mode='restart', startingwfc/pot='file', electron_maxstep=1, dmft=.true.
Steps 1-5 repeat until the interacting charge density (and the Wannier occupations) converge.
Energies in eV (QE eigenvalues are converted), DMFT temperature beta in 1/eV; use Fermi-Dirac smearing in QE with
degauss = k_B T (Ry), so that with Sigma = 0 the DMFT density equals the DFT density exactly (Delta N = 0).
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET

import numpy as np

HA = 27.211386245988
RY = HA / 2


# =============================================================================== QE / W90 file helpers
def read_qe_xml(path):
    """eigenvalues (eV, absolute), occupations f (per spin, 0..1), k points (cart, 2pi/alat), Fermi energy (eV)."""
    r = ET.parse(path).getroot()
    ks = r.findall(".//ks_energies")
    E = np.array([[float(x) for x in k.find("eigenvalues").text.split()] for k in ks]) * HA
    f = np.array([[float(x) for x in k.find("occupations").text.split()] for k in ks])
    kp = np.array([[float(x) for x in k.find("k_point").text.split()] for k in ks])
    ef = float(r.find(".//band_structure/fermi_energy").text) * HA
    etot = r.find(".//total_energy/etot")
    return dict(E=E, f=f, k_cart=kp, ef=ef, etot=float(etot.text) * HA if etot is not None else np.nan)


def read_u_matrix(path):
    """Wannier90 seed_u.mat / seed_u_dis.mat -> array (nk, nrow, ncol) with U[k, m, n] (m: band, n: Wannier)."""
    with open(path) as fh:
        fh.readline()
        nk, a, b = (int(x) for x in fh.readline().split())
        vals = []
        for line in fh:
            sp = line.split()
            if len(sp) == 2:
                vals.append(float(sp[0]) + 1j * float(sp[1]))
    vals = np.array(vals).reshape(nk, b, a)          # written with the row index fastest
    return np.transpose(vals, (0, 2, 1))


def mp_kgrid(n):
    n1, n2, n3 = n
    return np.array([(i / n1, j / n2, l / n3) for i in range(n1) for j in range(n2) for l in range(n3)])


def write_dmft_h5(path, band_window, delta_N):
    """Write the archive read by pw.x (band_window 1-based (lo, hi); delta_N (nk, nb, nb) complex, per spin)."""
    import h5py
    nk, nb, _ = delta_N.shape
    with h5py.File(path, "w") as h:
        bw = np.zeros((1, nk, 2))
        bw[0, :, 0], bw[0, :, 1] = band_window
        h.create_dataset("dft_misc_input/band_window", data=bw)
        d = np.zeros((nk, nb, nb, 2))
        d[..., 0], d[..., 1] = delta_N.real, delta_N.imag
        h.create_dataset("dft_update/delta_N", data=d)


# =============================================================================== Bloch-window DMFT update
def window_density_matrix(eps, V, sigma_w, mu, beta, n_iw):
    """N^DMFT(k) in the band window (per spin), exact 1/(iw)^2 tail.
    eps: (nk, nb) band energies (eV, relative to the reference); V: (nk, nb, nw) window->Wannier;
    sigma_w: (nw, n_iw) Wannier self-energy *including* the double-counting shift (diagonal)."""
    iw = 1j * np.pi * (2 * np.arange(n_iw) + 1) / beta
    nk, nb = eps.shape
    eye = np.eye(nb)
    wn = iw.imag
    tail = np.sum(1.0 / wn**2)
    sinf = np.real(sigma_w[:, -1])
    out = np.empty((nk, nb, nb), complex)
    for k in range(nk):
        Vk = V[k]
        Sb = np.einsum("ma,aw,na->wmn", Vk, sigma_w, Vk.conj())           # (n_iw, nb, nb)
        G = np.linalg.inv((iw[:, None, None] + mu) * eye - np.diag(eps[k])[None] - Sb)
        c2 = np.diag(eps[k]) - mu * eye + Vk @ np.diag(sinf) @ Vk.conj().T
        Gs = G + np.conj(np.swapaxes(G, -1, -2))
        N = 0.5 * eye + (Gs.sum(0) + 2 * c2 * tail - c2 * beta**2 / 4) / beta
        out[k] = 0.5 * (N + N.conj().T)
    return out


def find_mu_window(eps, V, sigma_w, beta, n_iw, n_target, lo=-10.0, hi=10.0, tol=1e-8):
    """mu such that 2 <tr N(k)>_k = n_target (both spins)."""
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        n = 2 * np.mean(np.real(np.trace(window_density_matrix(eps, V, sigma_w, mid, beta, n_iw), axis1=1, axis2=2)))
        lo, hi = (mid, hi) if n < n_target else (lo, mid)
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


# =============================================================================== driver
class QEDMFT:
    """Charge self-consistent DFT+DMFT driver for pw.x + wannier90.

    scf_input   : text of a working pw.x scf input (any structure); it must use occupations='smearing',
                  smearing='fermi-dirac', degauss = k_B T in Ry (consistent with beta), outdir and prefix.
    kmesh       : (n1, n2, n3) uniform mesh for nscf / Wannier90 / DMFT
    band_window : (lo, hi) 1-based band indices of the correlated window
    projections : Wannier90 projection lines, e.g. ["V:dxz;dyz;dxy"]
    U, J, beta  : interaction (eV) and inverse temperature (1/eV); solver: 'hubbard1' | 'ed' | 'cthyb'
    pw_cmd      : e.g. "mpirun -np 4 pw.x -nk 4"
    """

    def __init__(self, workdir, scf_input, kmesh, band_window, projections, U=0.0, J=0.0, beta=40.0,
                 nbnd=None, solver="ed", dc="FLL", n_iw=1024, seed="wann", pw_cmd="pw.x", w90_cmd="wannier90.x",
                 pw2wan_cmd="pw2wannier90.x", mixing_beta=0.3, solver_kwargs=None, dis=None, verbose=True):
        self.wd = os.path.abspath(workdir)
        os.makedirs(self.wd, exist_ok=True)
        self.scf_input = scf_input
        self.kmesh = tuple(kmesh)
        self.kfrac = mp_kgrid(self.kmesh)
        self.lo, self.hi = band_window
        self.nb = self.hi - self.lo + 1
        self.proj = projections
        self.U, self.J, self.beta = U, J, beta
        self.solver, self.dc, self.n_iw = solver, dc, n_iw
        self.nbnd = nbnd or self.hi + 6
        self.seed = seed
        self.pw_cmd, self.w90_cmd, self.pw2wan_cmd = pw_cmd, w90_cmd, pw2wan_cmd
        self.mixing_beta = mixing_beta
        self.solver_kwargs = solver_kwargs or {}
        self.dis = dis
        self.verbose = verbose
        m = re.search(r"prefix\s*=\s*'([^']+)'", scf_input)
        self.prefix = m.group(1)
        m = re.search(r"outdir\s*=\s*'([^']+)'", scf_input)
        self.outdir = m.group(1)
        self.history = []
        self.sigma = None

    # --------------------------------------------------------------- running codes
    def _run(self, cmd, stdin=None, log=None):
        with open(os.path.join(self.wd, log), "w") if log else open(os.devnull, "w") as out:
            p = subprocess.run(cmd, shell=True, cwd=self.wd, stdin=open(os.path.join(self.wd, stdin)) if stdin else None,
                               stdout=out, stderr=subprocess.STDOUT)
        return p.returncode

    def _klist(self):
        nk = len(self.kfrac)
        return "K_POINTS crystal\n%d\n" % nk + "".join("  %.12f %.12f %.12f %.6e\n" % (a, b, c, 1.0 / nk) for a, b, c in self.kfrac)

    def _system_extra(self, s, extra):
        return re.sub(r"(&system[^/]*?)(\n\s*/)", lambda m: m.group(1) + ",\n    " + extra + m.group(2), s, count=1,
                      flags=re.IGNORECASE | re.DOTALL)

    def write_inputs(self):
        s = self.scf_input
        body, _ = s.split("K_POINTS")[0], None
        nscf = re.sub(r"calculation\s*=\s*'scf'", "calculation='nscf', verbosity='high'", body)
        nscf = self._system_extra(nscf, f"nbnd={self.nbnd}, nosym=.true., noinv=.true., dmft=.true.")
        open(os.path.join(self.wd, "nscf.in"), "w").write(nscf + self._klist())
        d = re.sub(r"calculation\s*=\s*'scf'", "calculation='scf', restart_mode='restart', verbosity='high'", body)
        d = self._system_extra(d, f"nbnd={self.nbnd}, nosym=.true., noinv=.true., dmft=.true., dmft_prefix='{self.seed}_dmft'")
        d = re.sub(r"&electrons[^/]*/", f"&electrons\n    electron_maxstep=1, mixing_beta={self.mixing_beta}, "
                                        f"startingwfc='file', startingpot='file'\n /", d, flags=re.IGNORECASE | re.DOTALL)
        open(os.path.join(self.wd, "dmft_scf.in"), "w").write(d + self._klist())
        open(os.path.join(self.wd, "scf.in"), "w").write(s)

    def write_wannier(self, lattice_rows_A, atoms_frac):
        """lattice in Angstrom (3x3 rows); atoms_frac: list of (symbol, (x, y, z))."""
        nb, nw = self.nb, self._nw()
        excl = [i for i in range(1, self.nbnd + 1) if not (self.lo <= i <= self.hi)]
        lines = [f"num_wann = {nw}", f"num_bands = {nb}", f"exclude_bands = {self._ranges(excl)}",
                 "num_iter = 200", "write_hr = true", "write_u_matrices = true",
                 f"mp_grid = {self.kmesh[0]} {self.kmesh[1]} {self.kmesh[2]}"]
        if self.dis:
            lines += [f"dis_win_min = {self.dis[0]}", f"dis_win_max = {self.dis[1]}",
                      f"dis_froz_min = {self.dis[2]}", f"dis_froz_max = {self.dis[3]}", "dis_num_iter = 500"]
        lines += ["begin unit_cell_cart", "ang"] + ["  %.10f %.10f %.10f" % tuple(r) for r in lattice_rows_A] + ["end unit_cell_cart"]
        lines += ["begin atoms_frac"] + ["  %s %.10f %.10f %.10f" % (s, *x) for s, x in atoms_frac] + ["end atoms_frac"]
        lines += ["begin projections"] + list(self.proj) + ["end projections"]
        lines += ["begin kpoints"] + ["  %.12f %.12f %.12f" % tuple(k) for k in self.kfrac] + ["end kpoints"]
        open(os.path.join(self.wd, f"{self.seed}.win"), "w").write("\n".join(lines) + "\n")
        open(os.path.join(self.wd, f"{self.seed}.pw2wan"), "w").write(
            f"&inputpp\n outdir='{self.outdir}'\n prefix='{self.prefix}'\n seedname='{self.seed}'\n"
            " write_amn=.true.\n write_mmn=.true.\n write_unk=.false.\n/\n")

    def _nw(self):
        n = 0
        orb = {"s": 1, "p": 3, "d": 5, "f": 7, "pz": 1, "px": 1, "py": 1, "dxy": 1, "dxz": 1, "dyz": 1, "dz2": 1,
               "dx2-y2": 1, "t2g": 3, "eg": 2}
        for p in self.proj:
            for o in p.split(":")[1].split(";"):
                n += orb[o.strip()]
        return n

    @staticmethod
    def _ranges(idx):
        out, i = [], 0
        while i < len(idx):
            j = i
            while j + 1 < len(idx) and idx[j + 1] == idx[j] + 1:
                j += 1
            out.append(f"{idx[i]}-{idx[j]}" if j > i else f"{idx[i]}")
            i = j + 1
        return ", ".join(out)

    def run_scf(self):
        return self._run(f"{self.pw_cmd} -in scf.in", log="scf.out")

    def run_nscf(self):
        self._run(f"{self.pw_cmd} -in nscf.in", log="nscf.out")
        return "interrupted for call of DMFT" in open(os.path.join(self.wd, "nscf.out")).read()

    def run_wannier(self):
        self._run(f"{self.w90_cmd} -pp {self.seed}", log="w90pp.out")
        self._run(f"{self.pw2wan_cmd} -in {self.seed}.pw2wan", log="pw2wan.out")
        self._run(f"{self.w90_cmd} {self.seed}", log="w90.out")
        return os.path.exists(os.path.join(self.wd, f"{self.seed}_u.mat"))

    # --------------------------------------------------------------- DMFT step
    def load_bands(self):
        x = read_qe_xml(os.path.join(self.wd, self.outdir, f"{self.prefix}.save", "data-file-schema.xml"))
        self.ef = x["ef"]
        eps = x["E"][:, self.lo - 1:self.hi] - self.ef
        f = x["f"][:, self.lo - 1:self.hi]
        U = read_u_matrix(os.path.join(self.wd, f"{self.seed}_u.mat"))
        p = os.path.join(self.wd, f"{self.seed}_u_dis.mat")
        V = np.einsum("kmi,kin->kmn", read_u_matrix(p), U) if os.path.exists(p) else U
        return eps, f, V, x

    def dmft_step(self, n_dmft_iter=20, mix=0.5):
        from ..dmft.dftdmft import DFTDMFT
        eps, f, V, x = self.load_bands()
        nk = eps.shape[0]
        Hw = np.einsum("kmi,km,kmj->kij", V.conj(), eps, V)                      # Wannier Hamiltonian H_W(k)
        n_win = 2 * f.sum() / nk
        m = DFTDMFT(Hk=Hw, corr=None, U=self.U, J=self.J, n_elec=n_win, beta=self.beta, n_iw=self.n_iw, dc=self.dc)
        res = m.solve(self.solver, n_iter=n_dmft_iter, mix=mix, sigma0=self.sigma, **self.solver_kwargs)
        self.sigma = res["sigma"]                                                 # Sigma_imp - Sigma_dc (Wannier)
        mu = find_mu_window(eps, V, self.sigma, self.beta, self.n_iw, n_win, lo=res["mu"] - 3, hi=res["mu"] + 3)
        N = window_density_matrix(eps, V, self.sigma, mu, self.beta, self.n_iw)
        dN = N - np.einsum("km,mn->kmn", f, np.eye(self.nb))
        write_dmft_h5(os.path.join(self.wd, f"{self.seed}_dmft.h5"), (self.lo, self.hi), dN)
        Nw = np.einsum("kmi,kmn,knj->kij", V.conj(), N, V).mean(0)
        info = dict(mu_dmft=mu, n_window=n_win, n_wannier=2 * np.real(np.diag(Nw)), dN_max=float(np.abs(dN).max()),
                    dN_trace=float(2 * np.real(np.trace(dN, axis1=1, axis2=2)).mean()), Z=res["Z"],
                    ef=self.ef, etot_dft=x["etot"])
        return info

    def run_dmft_scf(self):
        self._run(f"{self.pw_cmd} -in dmft_scf.in", log="dmft_scf.out")
        out = open(os.path.join(self.wd, "dmft_scf.out")).read()
        ef = re.findall(r"the Fermi energy is\s+([-\d.]+)", out)
        acc = re.findall(r"estimated scf accuracy\s+<\s+([-\d.Ee+]+)", out)
        et = re.findall(r"total energy\s+=\s+([-\d.]+)\s+Ry", out)
        return dict(ef=float(ef[-1]) if ef else np.nan, scf_accuracy_Ry=float(acc[-1]) if acc else np.nan,
                    etot_Ry=float(et[-1]) if et else np.nan, ok="expecting DMFT hdf5 archive" in out)

    # --------------------------------------------------------------- full cycle
    def cycle(self, n_csc=10, n_dmft_first=20, n_dmft_per=4, tol_occ=1e-4, lattice_rows_A=None, atoms_frac=None,
              skip_scf=False, tol_rho_Ry=1e-6, tol_E_Ry=2e-5, min_csc=3):
        """Converged when QE's estimated scf accuracy < tol_rho_Ry, |dE| < tol_E_Ry and the change of the
        correlated-window charge and Wannier occupations < tol_occ (symmetry can pin the latter, so all are used)."""
        self.write_inputs()
        if lattice_rows_A is not None:
            self.write_wannier(lattice_rows_A, atoms_frac)
        if not skip_scf:
            self.run_scf()
        prev = None
        prevE, prevN = None, None
        for it in range(n_csc):
            if not self.run_nscf():
                raise RuntimeError("nscf did not stop for DMFT (check dmft=.true. and the QE build with HDF5)")
            if not self.run_wannier():
                raise RuntimeError("Wannier90 failed (see w90.out)")
            info = self.dmft_step(n_dmft_first if it == 0 else n_dmft_per)
            sc = self.run_dmft_scf()
            info.update({"scf_" + k: v for k, v in sc.items()})
            occ = info["n_wannier"]
            info["d_occ"] = float(np.max(np.abs(occ - prev))) if prev is not None else np.nan
            info["d_nwin"] = float(abs(info["n_window"] - prevN)) if prevN is not None else np.nan
            info["dE_Ry"] = float(abs(sc["etot_Ry"] - prevE)) if prevE is not None else np.nan
            self.history.append(info)
            if self.verbose:
                print(f"CSC {it}: n_win={info['n_window']:.5f} n_W={np.round(occ, 5)}  Z={np.round(info['Z'], 3)}  "
                      f"|dN|max={info['dN_max']:.2e}  E_F={sc['ef']:.4f} eV  scf_acc={sc['scf_accuracy_Ry']:.2e} Ry  "
                      f"E={sc['etot_Ry']:.6f} Ry", flush=True)
            if (it + 1 >= min_csc and prev is not None and info["d_occ"] < tol_occ and info["d_nwin"] < tol_occ
                    and sc["scf_accuracy_Ry"] < tol_rho_Ry and info["dE_Ry"] < tol_E_Ry):
                break
            prev, prevE, prevN = occ, sc["etot_Ry"], info["n_window"]
        return self.history
