"""Band models for Hall-type transport.  Every model provides vectorised

    H(kx, ky)   -> (..., nb, nb) complex
    dH(kx, ky)  -> (dH/dkx, dH/dky), same shape
    operators   : dict name -> (nb, nb) matrix (spin / orbital / quadrupole ...)
    kind        : 'continuum' (integrate over the plane, polar grid) or
                  'lattice'   (integrate over the Brillouin zone)

Units: hbar = 1, lattice constant = 1.
"""
from __future__ import annotations

import numpy as np

s0 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def _bc(*arrs):
    return np.broadcast_arrays(*[np.asarray(a, float) for a in arrs])


class BandModel:
    kind = "continuum"
    nb = 2
    operators: dict = {}

    def H(self, kx, ky):
        raise NotImplementedError

    def dH(self, kx, ky):
        raise NotImplementedError

    def eig(self, kx, ky):
        return np.linalg.eigh(self.H(kx, ky))


# ----------------------------------------------------------------------------
class MassiveDirac(BandModel):
    """H = v (kx sx + ky sy) + Delta sz   (Sinitsyn et al., PRB 75, 045315 (2007))."""

    def __init__(self, v: float = 1.0, Delta: float = 0.5):
        self.v, self.Delta = float(v), float(Delta)
        self.operators = {"sz/2": sz / 2, "sx/2": sx / 2, "sy/2": sy / 2}

    def H(self, kx, ky):
        kx, ky = _bc(kx, ky)
        return (self.v * (kx[..., None, None] * sx + ky[..., None, None] * sy)
                + self.Delta * sz)

    def dH(self, kx, ky):
        kx, ky = _bc(kx, ky)
        sh = kx.shape + (2, 2)
        return np.broadcast_to(self.v * sx, sh).copy(), np.broadcast_to(self.v * sy, sh).copy()

    # analytic references --------------------------------------------------
    def kF(self, mu):
        return np.sqrt(mu**2 - self.Delta**2) / self.v

    def lower_band_tail(self, kmax):
        """Intrinsic AHE (e^2/h) of the filled lower band outside |k| > kmax, to be
        added to a cut-off continuum integral: -Delta / (2 sqrt(Delta^2 + v^2 kmax^2))."""
        return -self.Delta / (2.0 * np.sqrt(self.Delta**2 + (self.v * kmax) ** 2))

    # Sinitsyn et al. PRB 75, 045315 (2007): all in units of e^2/h (hbar = 1)
    def ref_intrinsic(self, mu):           # Eq. (58)
        return -self.Delta / (2 * mu)

    def ref_side_jump(self, mu):           # Eq. (65)  (side-jump velocity)
        k2 = (self.v * self.kF(mu)) ** 2
        return -self.Delta * k2 / (mu * (k2 + 4 * self.Delta**2))

    def ref_anomalous_distribution(self, mu):  # Eq. (67) = Eq. (65)
        return self.ref_side_jump(mu)

    def ref_intrinsic_skew(self, mu):      # Gaussian-disorder (4th order) skew part of Eq. (72)
        k2 = (self.v * self.kF(mu)) ** 2
        return -3 * self.Delta * k2**2 / (2 * mu * (4 * self.Delta**2 + k2) ** 2)

    def ref_skew_V1(self, mu, n_i, V0, V1):  # conventional skew, second term of Eq. (72)
        k2 = (self.v * self.kF(mu)) ** 2
        return -(V1**3 / (n_i * V0**4)) * k2**2 * self.Delta / (4 * self.Delta**2 + k2) ** 2

    def ref_total_gaussian(self, mu):      # Eq. (B8), non-crossing approximation
        c = self.Delta / mu
        return -4 * c * (1 + c**2) / (1 + 3 * c**2) ** 2


class DiracQuadraticMass(BandModel):
    """H = v (kx sx + ky sy) + (Delta + b k^2) sz  (continuum limit of QWZ / BHZ near Gamma, b = 1/2)."""

    def __init__(self, v: float = 1.0, Delta: float = 0.3, b: float = 0.5):
        self.v, self.Delta, self.b = float(v), float(Delta), float(b)
        self.operators = {"sz/2": sz / 2}

    def H(self, kx, ky):
        kx, ky = _bc(kx, ky)
        k2 = kx**2 + ky**2
        return (self.v * (kx[..., None, None] * sx + ky[..., None, None] * sy)
                + (self.Delta + self.b * k2)[..., None, None] * sz)

    def dH(self, kx, ky):
        kx, ky = _bc(kx, ky)
        return (self.v * sx + 2 * self.b * kx[..., None, None] * sz,
                self.v * sy + 2 * self.b * ky[..., None, None] * sz)


class TiltedDirac(BandModel):
    """H = t kx 1 + v (kx sx + ky sy) + Delta sz : a single tilted massive Dirac cone
    (minimal model with a Berry curvature dipole, Sodemann & Fu PRL 115, 216806 (2015))."""

    def __init__(self, v: float = 1.0, Delta: float = 0.1, t: float = 0.2):
        self.v, self.Delta, self.t = float(v), float(Delta), float(t)
        self.operators = {"sz/2": sz / 2}

    def H(self, kx, ky):
        kx, ky = _bc(kx, ky)
        return (self.t * kx[..., None, None] * np.eye(2) + self.v * (kx[..., None, None] * sx + ky[..., None, None] * sy)
                + self.Delta * sz)

    def dH(self, kx, ky):
        kx, ky = _bc(kx, ky)
        sh = kx.shape + (2, 2)
        return (np.broadcast_to(self.t * np.eye(2) + self.v * sx, sh).copy(),
                np.broadcast_to(self.v * sy, sh).copy())


class Rashba2DEG(BandModel):
    """H = k^2/2m + alpha (sx ky - sy kx) + Delta_ex sz.
    Delta_ex = 0: Rashba spin Hall model (Sinova et al. PRL 92, 126603 (2004);
    vertex cancellation: Inoue, Bauer, Molenkamp PRB 70, 041303 (2004)).
    Delta_ex != 0: ferromagnetic Rashba AHE model."""

    def __init__(self, m: float = 1.0, alpha: float = 0.3, Delta_ex: float = 0.0):
        self.m, self.alpha, self.Dex = float(m), float(alpha), float(Delta_ex)
        self.operators = {"sz/2": sz / 2, "sx/2": sx / 2, "sy/2": sy / 2}

    def H(self, kx, ky):
        kx, ky = _bc(kx, ky)
        e = (kx**2 + ky**2) / (2 * self.m)
        return (e[..., None, None] * s0 + self.alpha * (ky[..., None, None] * sx - kx[..., None, None] * sy)
                + self.Dex * sz)

    def dH(self, kx, ky):
        kx, ky = _bc(kx, ky)
        hx = (kx / self.m)[..., None, None] * s0 - self.alpha * sy
        hy = (ky / self.m)[..., None, None] * s0 + self.alpha * sx
        return hx, hy


# ----------------------------------------------------------------------------
class QWZ(BandModel):
    """Qi-Wu-Zhang model  d = (sin kx, sin ky, m + cos kx + cos ky).
    Chern number C = -sign(m)... for 0<|m|<2 and C=0 for |m|>2 (Qi, Wu, Zhang PRB 74, 085308 (2006))."""
    kind = "lattice"

    def __init__(self, m: float = 1.0):
        self.mass = float(m)
        self.operators = {"sz/2": sz / 2}

    def H(self, kx, ky):
        kx, ky = _bc(kx, ky)
        return (np.sin(kx)[..., None, None] * sx + np.sin(ky)[..., None, None] * sy
                + (self.mass + np.cos(kx) + np.cos(ky))[..., None, None] * sz)

    def dH(self, kx, ky):
        kx, ky = _bc(kx, ky)
        return (np.cos(kx)[..., None, None] * sx - np.sin(kx)[..., None, None] * sz,
                np.cos(ky)[..., None, None] * sy - np.sin(ky)[..., None, None] * sz)


class Haldane(BandModel):
    """Haldane model on the honeycomb lattice (PRL 61, 2015 (1988)).
    t1 NN hopping, t2 e^{i phi} NNN, M sublattice mass; nearest-neighbour distance 1."""
    kind = "lattice"

    def __init__(self, t1=1.0, t2=0.1, phi=np.pi / 2, M=0.0):
        self.t1, self.t2, self.phi, self.M = t1, t2, phi, M
        self.a = np.array([[1.0, 0.0], [-0.5, np.sqrt(3) / 2], [-0.5, -np.sqrt(3) / 2]])  # NN vectors
        self.b = np.array([[0.0, np.sqrt(3)], [-1.5, -np.sqrt(3) / 2], [1.5, -np.sqrt(3) / 2]])  # NNN
        self.operators = {"sz/2": sz / 2}
        # reciprocal lattice for BZ integration (primitive vectors a1, a2)
        a1 = self.a[0] - self.a[1]
        a2 = self.a[0] - self.a[2]
        A = np.array([a1, a2])
        self.recip = 2 * np.pi * np.linalg.inv(A).T  # rows b1, b2

    def _parts(self, kx, ky):
        kx, ky = _bc(kx, ky)
        ka = kx[..., None] * self.a[:, 0] + ky[..., None] * self.a[:, 1]
        kb = kx[..., None] * self.b[:, 0] + ky[..., None] * self.b[:, 1]
        return kx, ky, ka, kb

    def H(self, kx, ky):
        kx, ky, ka, kb = self._parts(kx, ky)
        d0 = 2 * self.t2 * np.cos(self.phi) * np.cos(kb).sum(-1)
        dx = self.t1 * np.cos(ka).sum(-1)
        dy = self.t1 * np.sin(ka).sum(-1)
        dz = self.M - 2 * self.t2 * np.sin(self.phi) * np.sin(kb).sum(-1)
        return (d0[..., None, None] * s0 + dx[..., None, None] * sx + dy[..., None, None] * sy
                + dz[..., None, None] * sz)

    def dH(self, kx, ky):
        kx, ky, ka, kb = self._parts(kx, ky)
        out = []
        for c in (0, 1):
            ac, bc = self.a[:, c], self.b[:, c]
            d0 = -2 * self.t2 * np.cos(self.phi) * (np.sin(kb) * bc).sum(-1)
            dx = -self.t1 * (np.sin(ka) * ac).sum(-1)
            dy = self.t1 * (np.cos(ka) * ac).sum(-1)
            dz = -2 * self.t2 * np.sin(self.phi) * (np.cos(kb) * bc).sum(-1)
            out.append(d0[..., None, None] * s0 + dx[..., None, None] * sx + dy[..., None, None] * sy
                       + dz[..., None, None] * sz)
        return tuple(out)


# ----------------------------------------------------------------------------
def p_orbital_L():
    """Angular-momentum matrices in the real (px, py, pz) basis: (L_a)_{bc} = -i eps_{abc}."""
    L = np.zeros((3, 3, 3), dtype=complex)
    eps = np.zeros((3, 3, 3))
    eps[0, 1, 2] = eps[1, 2, 0] = eps[2, 0, 1] = 1
    eps[0, 2, 1] = eps[2, 1, 0] = eps[1, 0, 2] = -1
    for a in range(3):
        L[a] = -1j * eps[a]
    return L


class POrbitalSquare(BandModel):
    """(px, py, pz) x spin on a square lattice (6 bands): sigma- and pi-bonding
    hoppings, atomic spin-orbit coupling lam L.S and an orbital-Rashba term
    gamma (inversion breaking, z-mirror broken) coupling p_z to p_x/p_y.

    Minimal model of the orbital Hall effect: the OHC is finite already for
    lam = 0, while the SHC ~ lam  (Bernevig, Hughes, Zhang PRL 95, 066601 (2005);
    Go et al. PRL 121, 086602 (2018)).  m_ex adds an exchange field (magnetic
    case): needed for T-odd responses such as the AHE and the orbital-quadrupole
    Hall current {Q, v}/2, which vanish identically with time reversal.

    operators: 'Lz', 'Sz', 'Jz', and the orbital-quadrupole components
    'Qxy' = (LxLy+LyLx)/2, 'Qx2-y2' = (Lx^2 - Ly^2)/2,
    'Q3z2-r2' = (3 Lz^2 - L^2)/(2 sqrt 3)."""
    kind = "lattice"
    nb = 6

    def __init__(self, t_sigma=1.0, t_pi=0.25, t_z=0.25, eps_z=0.0, lam=0.0, gamma=0.0, m_ex=0.0):
        self.ts, self.tp, self.tz, self.ez = t_sigma, t_pi, t_z, eps_z
        self.m_ex = m_ex  # exchange field  m_ex * 2 S_z  (breaks time reversal)
        self.lam, self.gamma = lam, gamma
        L = p_orbital_L()
        self.L = L
        I2, I3 = np.eye(2), np.eye(3)
        S = np.array([sx, sy, sz]) / 2
        self.LS = sum(np.kron(L[a], S[a]) for a in range(3))
        Lk = [np.kron(L[a], I2) for a in range(3)]
        Sk = [np.kron(I3, S[a]) for a in range(3)]
        L2 = sum(l @ l for l in L)
        Q = {"Qxy": (L[0] @ L[1] + L[1] @ L[0]) / 2,
             "Qx2-y2": (L[0] @ L[0] - L[1] @ L[1]) / 2,
             "Q3z2-r2": (3 * L[2] @ L[2] - L2) / (2 * np.sqrt(3))}
        self.operators = {"Lx": Lk[0], "Ly": Lk[1], "Lz": Lk[2], "Sz": Sk[2], "Sx": Sk[0],
                          "Sy": Sk[1], "Jz": Lk[2] + Sk[2]}
        for k, v in Q.items():
            self.operators[k] = np.kron(v, I2)

    def _h3(self, kx, ky):
        kx, ky = _bc(kx, ky)
        h = np.zeros(kx.shape + (3, 3), dtype=complex)
        cx, cy = np.cos(kx), np.cos(ky)
        h[..., 0, 0] = 2 * self.ts * cx - 2 * self.tp * cy
        h[..., 1, 1] = 2 * self.ts * cy - 2 * self.tp * cx
        h[..., 2, 2] = self.ez - 2 * self.tz * (cx + cy)
        # orbital Rashba: <pz|H|px> ~ 2i gamma sin kx , <pz|H|py> ~ 2i gamma sin ky
        h[..., 2, 0] = 2j * self.gamma * np.sin(kx)
        h[..., 0, 2] = np.conj(h[..., 2, 0])
        h[..., 2, 1] = 2j * self.gamma * np.sin(ky)
        h[..., 1, 2] = np.conj(h[..., 2, 1])
        return h

    def _dh3(self, kx, ky):
        kx, ky = _bc(kx, ky)
        sxk, syk = np.sin(kx), np.sin(ky)
        hx = np.zeros(kx.shape + (3, 3), dtype=complex)
        hy = np.zeros_like(hx)
        hx[..., 0, 0] = -2 * self.ts * sxk
        hx[..., 1, 1] = 2 * self.tp * sxk
        hx[..., 2, 2] = 2 * self.tz * sxk
        hx[..., 2, 0] = 2j * self.gamma * np.cos(kx)
        hx[..., 0, 2] = np.conj(hx[..., 2, 0])
        hy[..., 0, 0] = 2 * self.tp * syk
        hy[..., 1, 1] = -2 * self.ts * syk
        hy[..., 2, 2] = 2 * self.tz * syk
        hy[..., 2, 1] = 2j * self.gamma * np.cos(ky)
        hy[..., 1, 2] = np.conj(hy[..., 2, 1])
        return hx, hy

    def H(self, kx, ky):
        h3 = self._h3(kx, ky)
        return np.kron(h3, np.eye(2)) + self.lam * self.LS + self.m_ex * 2 * self.operators['Sz']

    def dH(self, kx, ky):
        hx, hy = self._dh3(kx, ky)
        return np.kron(hx, np.eye(2)), np.kron(hy, np.eye(2))


class WannierTB(BandModel):
    """Generic tight-binding model from a Wannier90 seedname_hr.dat (2D use:
    kz = 0).  Real-space lattice vectors 'lat' (3x3, rows) in Angstrom; k in 1/A.
    Operators may be supplied by the user (e.g. Lz in the Wannier basis)."""
    kind = "lattice"

    def __init__(self, hr_file: str, lat, operators: dict | None = None):
        self.lat = np.asarray(lat, float)
        with open(hr_file) as f:
            f.readline()
            nw = int(f.readline())
            nr = int(f.readline())
            deg = []
            while len(deg) < nr:
                deg += [int(x) for x in f.readline().split()]
            data = np.loadtxt(f)
        self.nb = nw
        R = data[:, :3].astype(int)
        m = data[:, 3].astype(int) - 1
        n = data[:, 4].astype(int) - 1
        h = data[:, 5] + 1j * data[:, 6]
        Rs = np.unique(R, axis=0)
        idx = {tuple(r): i for i, r in enumerate(Rs)}
        HR = np.zeros((len(Rs), nw, nw), complex)
        degs = np.repeat(np.array(deg), nw * nw)
        for r, a, b, val, d in zip(R, m, n, h, degs):
            HR[idx[tuple(r)], a, b] = val / d
        self.R = Rs @ self.lat  # Cartesian
        self.HR = HR
        self.operators = operators or {}
        self.recip = 2 * np.pi * np.linalg.inv(self.lat).T[:2, :2]

    def H(self, kx, ky):
        kx, ky = _bc(kx, ky)
        ph = np.exp(1j * (kx[..., None] * self.R[:, 0] + ky[..., None] * self.R[:, 1]))
        return np.einsum("...r,rab->...ab", ph, self.HR)

    def dH(self, kx, ky):
        kx, ky = _bc(kx, ky)
        ph = np.exp(1j * (kx[..., None] * self.R[:, 0] + ky[..., None] * self.R[:, 1]))
        return (np.einsum("...r,rab->...ab", 1j * self.R[:, 0] * ph, self.HR),
                np.einsum("...r,rab->...ab", 1j * self.R[:, 1] * ph, self.HR))
