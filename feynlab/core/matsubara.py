"""Matsubara-frequency machinery.

Conventions (used everywhere in feynlab)
----------------------------------------
* Imaginary-time Green's function  G(tau) = -<T c(tau) c^dagger(0)>,  0 <= tau < beta.
* Fermionic Matsubara frequencies  w_n = (2n+1) pi / beta.
* Fourier pair:  G(iw_n) = int_0^beta dtau e^{i w_n tau} G(tau)
                 G(tau)  = (1/beta) sum_n e^{-i w_n tau} G(iw_n).
* A free level eps (measured from mu):  G0(iw) = 1/(iw - eps),
  G0(tau) = -exp(-eps tau) / (1 + exp(-beta eps)).
* High-frequency tail of any fermionic G:  G(iw) -> 1/(iw) + c2/(iw)^2 + c3/(iw)^3.
"""
from __future__ import annotations

import numpy as np


class MatsubaraGrid:
    """Fermionic Matsubara grid with tail-corrected Fourier transforms."""

    def __init__(self, beta: float, n_iw: int = 1024, n_tau: int | None = None):
        if beta <= 0:
            raise ValueError("beta must be positive")
        self.beta = float(beta)
        self.n_iw = int(n_iw)
        self.n = np.arange(n_iw)
        self.iw = 1j * (2 * self.n + 1) * np.pi / self.beta  # positive frequencies only
        self.n_tau = int(n_tau) if n_tau is not None else 4 * n_iw + 1
        if (self.n_tau - 1) < 2 * self.n_iw:
            raise ValueError('n_tau - 1 must be >= 2*n_iw')
        self.tau = np.linspace(0.0, self.beta, self.n_tau)

    # ------------------------------------------------------------------
    def iw_full(self) -> np.ndarray:
        """Symmetric grid w_{-N}..w_{N-1}."""
        n = np.arange(-self.n_iw, self.n_iw)
        return 1j * (2 * n + 1) * np.pi / self.beta

    # ------------------------------------------------------------------
    def tau_to_iw(self, g_tau: np.ndarray) -> np.ndarray:
        """G(iw_n), n=0..n_iw-1, from G(tau) sampled on the uniform grid self.tau.

        The integral of the piecewise-linear interpolant is evaluated *exactly*
        and vectorised with FFTs (attenuation-factor method).  Works on the last
        axis, so g_tau may have shape (..., n_tau)."""
        g_tau = np.asarray(g_tau)
        M = self.n_tau - 1
        h = self.beta / M
        w = self.iw.imag
        n = self.n
        jj = np.arange(M)
        ph = np.exp(1j * np.pi * jj / M)  # e^{i pi j / M}
        # S1_n = sum_{j=0}^{M-1} g_j e^{i w_n t_j};  S2_n = sum_{j=0}^{M-1} g_{j+1} e^{i w_n t_j}
        b1 = g_tau[..., :M] * ph
        b2 = g_tau[..., 1:M + 1] * ph
        F1 = np.fft.ifft(b1, axis=-1) * M
        F2 = np.fft.ifft(b2, axis=-1) * M
        idx = n % M
        S1 = F1[..., idx]
        S2 = F2[..., idx]
        iwv = 1j * w
        th = w * h
        eith = np.exp(1j * th)
        cb = eith / iwv - (eith - 1.0) / (iwv**2 * h)
        ca = (eith - 1.0) / iwv - cb
        return S1 * ca + S2 * cb

    def iw_to_tau(self, g_iw: np.ndarray, c1: float = 1.0, c2: float = 0.0,
                  c3: float = 0.0, tau: np.ndarray | None = None) -> np.ndarray:
        """G(tau) from G(iw_n>=0) assuming G(-iw) = conj(G(iw)) (real G(tau)).
        Known tail c1/iw + c2/(iw)^2 + c3/(iw)^3 is subtracted and added back
        analytically.  If tau is None the uniform grid self.tau is used (FFT);
        otherwise an explicit sum is done.  Works on the last axis."""
        beta = self.beta
        iw = self.iw
        g_iw = np.asarray(g_iw)
        tail = c1 / iw + c2 / iw**2 + c3 / iw**3
        diff = g_iw - tail
        if tau is None:
            tau = self.tau
            M = self.n_tau - 1
            N = self.n_iw
            if M < 2 * N:
                raise ValueError("n_tau-1 must be >= 2 n_iw")
            a = np.zeros(diff.shape[:-1] + (M,), dtype=complex)
            a[..., :N] = diff
            a[..., M - N:] = np.conj(diff[..., ::-1])  # n = -N..-1  <->  conj of n'=N-1..0
            F = np.fft.fft(a, axis=-1)  # sum_m a_m e^{-2 pi i m j / M}
            jj = np.arange(M)
            g = np.empty(diff.shape[:-1] + (M + 1,))
            g[..., :M] = np.real(np.exp(-1j * np.pi * jj / M) * F) / beta
            # tau = beta^- point: explicit sum
            g[..., M] = (2.0 / beta) * np.real(np.sum(diff * np.exp(-iw * beta), axis=-1))
        else:
            tau = np.asarray(tau, dtype=float)
            phase = np.exp(-np.outer(tau, iw))
            g = (2.0 / beta) * np.real(diff @ phase.T)
        g = g + (-0.5 * c1) + c2 * (2 * tau - beta) / 4.0 + c3 * (beta * tau - tau**2) / 4.0
        return g

    # ------------------------------------------------------------------
    # bosonic transforms (nu_m = 2 pi m / beta, m = 0..n_iw-1)
    # ------------------------------------------------------------------
    @property
    def iv(self):
        return 1j * 2 * np.pi * self.n / self.beta

    def tau_to_iv(self, f_tau: np.ndarray) -> np.ndarray:
        """Exact integral of the piecewise-linear interpolant of a bosonic f(tau)."""
        f_tau = np.asarray(f_tau)
        M = self.n_tau - 1
        h = self.beta / M
        nu = self.iv.imag
        F1 = np.fft.ifft(f_tau[..., :M], axis=-1) * M
        F2 = np.fft.ifft(f_tau[..., 1:M + 1], axis=-1) * M
        idx = self.n % M
        S1 = F1[..., idx]
        S2 = F2[..., idx]
        out = np.empty(S1.shape, dtype=complex)
        nz = nu != 0
        iv = 1j * nu[nz]
        th = nu[nz] * h
        e = np.exp(1j * th)
        cb = e / iv - (e - 1.0) / (iv**2 * h)
        ca = (e - 1.0) / iv - cb
        out[..., nz] = S1[..., nz] * ca + S2[..., nz] * cb
        out[..., ~nz] = (S1[..., ~nz] + S2[..., ~nz]) * h / 2.0
        return out

    def iv_to_tau(self, f_iv: np.ndarray) -> np.ndarray:
        """f(tau_j) = (1/beta)[F_0 + 2 Re sum_{m>=1} e^{-i nu_m tau} F_m] on self.tau
        (no tail correction: use for functions decaying faster than 1/nu^2)."""
        f_iv = np.asarray(f_iv)
        M = self.n_tau - 1
        N = self.n_iw
        a = np.zeros(f_iv.shape[:-1] + (M,), dtype=complex)
        a[..., :N] = f_iv
        a[..., M - N + 1:] = np.conj(f_iv[..., 1:][..., ::-1])
        F = np.fft.fft(a, axis=-1)
        out = np.empty(f_iv.shape[:-1] + (M + 1,))
        out[..., :M] = F.real / self.beta
        out[..., M] = out[..., 0]
        return out

    # ------------------------------------------------------------------
    @staticmethod
    def tail_coefficients(iw: np.ndarray, g_iw: np.ndarray, n_fit: int = 40):
        """Least-squares estimate of c1,c2,c3 from the last n_fit frequencies."""
        x = iw[-n_fit:]
        y = g_iw[-n_fit:]
        A = np.vstack([1 / x, 1 / x**2, 1 / x**3]).T
        c, *_ = np.linalg.lstsq(np.vstack([A.real, A.imag]),
                                np.concatenate([y.real, y.imag]), rcond=None)
        return c

    def density(self, g_iw: np.ndarray, c2: float | None = None) -> float:
        """n = G(beta^-) * (-1) = -G(tau=beta) per spin."""
        if c2 is None:
            c2 = 0.0
        g_beta = self.iw_to_tau(g_iw, 1.0, c2, 0.0, tau=np.array([self.beta]))[0]
        return float(-g_beta)


def fermi(e, beta):
    """Numerically safe Fermi function."""
    x = np.clip(beta * np.asarray(e, dtype=float), -700, 700)
    return 1.0 / (np.exp(x) + 1.0)


def bose(e, beta):
    x = np.clip(beta * np.asarray(e, dtype=float), 1e-300, 700)
    return 1.0 / np.expm1(x)


def g0_tau_level(eps: float, tau: np.ndarray, beta: float) -> np.ndarray:
    """Free fermion G0(tau) for a single level eps (relative to mu), 0<=tau<beta."""
    tau = np.asarray(tau, dtype=float)
    if eps >= 0:
        return -np.exp(-eps * tau) / (1.0 + np.exp(-beta * eps))
    return -np.exp(eps * (beta - tau)) / (np.exp(beta * eps) + 1.0)
