"""Lattice polarons (Holstein and Su-Schrieffer-Heeger) by exact diagonalization in a variational
(translation-invariant) Hilbert space - Bonca, Trugman & Batistic, PRB 60, 1633 (1999).

1D chain, one electron, dispersionless phonons omega:
    Holstein:  H = -t sum (c+_i c_{i+1} + h.c.) + w sum b+ b - g sum_i n_i (b+_i + b_i)
    SSH (phonon-modulated hopping, Marchand et al. PRL 105, 266605 (2010)):
               H_int = g sum_i (c+_i c_{i+1} + h.c.) (b+_{i+1} + b_{i+1} - b+_i - b_i)
The basis is generated from |electron at 0, no phonons> by n_gen applications of the off-diagonal part of H;
states are phonon configurations relative to the electron (translation symmetry -> momentum k).  Convergence
in n_gen is exponential.  Returns the ground-state energy E(k), quasiparticle weight Z(k) = |<k, 0 ph|psi>|^2,
mean phonon number and the effective mass m*/m from E(k) near k = 0.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh


def _norm(conf):
    """Canonical relative phonon configuration: tuple of (rel_pos, n) with n > 0, sorted."""
    return tuple(sorted((p, n) for p, n in conf.items() if n > 0))


class LatticePolaronVED:
    def __init__(self, model="holstein", t=1.0, w=1.0, g=1.0, n_gen=12, max_ph_site=None):
        self.model, self.t, self.w, self.g = model, float(t), float(w), float(g)
        self.n_gen = n_gen
        self.max_ph = max_ph_site
        self._build()

    # -------------------------------------------------------------- basis + matrix elements (k independent)
    def _moves(self, conf):
        """Off-diagonal moves from a relative configuration.  Returns list of (new_conf, amplitude, shift),
        where shift = electron displacement (the relative frame moves with the electron)."""
        t, g = self.t, self.g
        out = []
        c = dict(conf)
        # electron hopping by +-1 (relative phonon positions shift by -+1)
        for d in (1, -1):
            nc = {p - d: n for p, n in c.items()}
            out.append((_norm(nc), -t, d))
        if self.model == "holstein":
            n0 = c.get(0, 0)
            if self.max_ph is None or n0 < self.max_ph:
                nc = dict(c); nc[0] = n0 + 1
                out.append((_norm(nc), -g * np.sqrt(n0 + 1), 0))
            if n0 > 0:
                nc = dict(c); nc[0] = n0 - 1
                out.append((_norm(nc), -g * np.sqrt(n0), 0))
        elif self.model == "ssh":
            # hop by d with a phonon created/annihilated on the two sites of the bond (i, i+d) (relative to old site)
            for d in (1, -1):
                # bond between old site 0 and new site d: coupling g (b_{max} - b_{min}) on bond (i, i+1)
                lo, hi = (0, 1) if d == 1 else (-1, 0)
                for site, sgn in ((hi, 1.0), (lo, -1.0)):
                    nsite = c.get(site, 0)
                    for create in (True, False):
                        if create:
                            if self.max_ph is not None and nsite >= self.max_ph:
                                continue
                            amp = g * sgn * np.sqrt(nsite + 1)
                            nn = nsite + 1
                        else:
                            if nsite == 0:
                                continue
                            amp = g * sgn * np.sqrt(nsite)
                            nn = nsite - 1
                        nc = dict(c); nc[site] = nn
                        nc = {p - d: n for p, n in nc.items()}
                        out.append((_norm(nc), amp, d))
        else:
            raise ValueError(self.model)
        return out

    def _build(self):
        start = ()
        index = {start: 0}
        frontier = [start]
        for _ in range(self.n_gen):
            new = []
            for conf in frontier:
                for nc, _, _ in self._moves(conf):
                    if nc not in index:
                        index[nc] = len(index)
                        new.append(nc)
            frontier = new
        self.index = index
        self.states = list(index.keys())
        rows, cols, amps, shifts = [], [], [], []
        diag = np.zeros(len(index))
        for conf, i in index.items():
            diag[i] = self.w * sum(n for _, n in conf)
            for nc, a, d in self._moves(conf):
                j = index.get(nc)
                if j is None:
                    continue
                rows.append(j); cols.append(i); amps.append(a); shifts.append(d)
        self._diag = diag
        self._rows, self._cols = np.array(rows), np.array(cols)
        self._amps, self._shifts = np.array(amps), np.array(shifts)

    def hamiltonian(self, k):
        """H(k) in the translation-invariant basis: |k; conf> = N^{-1/2} sum_j e^{ikj} |electron at j; conf>.
        A move that displaces the electron by d acquires the phase e^{-ikd}."""
        vals = self._amps * np.exp(-1j * k * self._shifts)
        n = len(self._diag)
        H = sp.coo_matrix((vals, (self._rows, self._cols)), shape=(n, n)).tocsr()
        H = H + sp.diags(self._diag)
        return 0.5 * (H + H.conj().T)

    def ground_state(self, k=0.0):
        H = self.hamiltonian(k)
        if H.shape[0] < 400:
            E, V = np.linalg.eigh(H.toarray())
            e, v = E[0], V[:, 0]
        else:
            E, V = eigsh(H, k=1, which="SA", tol=1e-12)
            e, v = E[0], V[:, 0]
        nph = float(np.sum(np.abs(v) ** 2 * self._diag) / self.w)
        return dict(E=float(e), Z=float(abs(v[0]) ** 2), n_ph=nph, dim=H.shape[0])

    def effective_mass(self, dk=0.02):
        """m*/m with m = 1/(2t) the bare band mass at k = 0 (E(k) ~ E0 + k^2/(2 m*))."""
        e0 = self.ground_state(0.0)["E"]
        e1 = self.ground_state(dk)["E"]
        e2 = self.ground_state(2 * dk)["E"]
        # E(+-dk) are equal by symmetry; fourth-order accurate second derivative
        d2 = (16 * e1 - e2 - 15 * e0) / (6 * dk**2)
        return 2 * self.t / d2


def holstein_weak_coupling(t, w, g, k=0.0, nq=20000):
    """Exact second-order (Rayleigh-Schroedinger) energy of the 1D Holstein polaron at momentum k."""
    q = 2 * np.pi * (np.arange(nq) + 0.5) / nq
    ek = -2 * t * np.cos(k)
    return ek - g**2 * np.mean(1.0 / (w - 2 * t * np.cos(k + q) - ek))


def ssh_weak_coupling(t, w, g, k=0.0, nq=20000):
    """Second-order energy for the SSH coupling above: vertex M(k,q) = 2 i g [sin(k+q) - sin k] for site phonons
    (one sign convention; |M|^2 is what enters)."""
    q = 2 * np.pi * (np.arange(nq) + 0.5) / nq
    ek = -2 * t * np.cos(k)
    M2 = (2 * g * (np.sin(k + q) - np.sin(k))) ** 2
    return ek - np.mean(M2 / (w - 2 * t * np.cos(k + q) - ek))


def holstein_strong_coupling(t, w, g):
    """Lang-Firsov: E = -g^2/w - 2 t exp(-g^2/w^2) - (2 t^2 / w) [...] ; leading two terms."""
    return -g**2 / w - 2 * t * np.exp(-(g**2) / w**2)
