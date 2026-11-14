"""Two-particle response functions on a 2D lattice with a (local or k-dependent) self-energy.

* bubble(eps_k, sigma, grid, mu)       chi0(q, i nu_m) = -(T/N) sum_k,n G(k,iw_n) G(k+q, iw_n + i nu_m)
                                        (per spin, >= 0 at nu = 0), computed with FFTs in (r, tau)
* rpa(chi0, U)                          spin / charge RPA:  chi_s = chi0/(1 - U chi0),  chi_c = chi0/(1 + U chi0)
* compressibility check                 chi0(q=0, nu=0) = -dn/dmu  for U = 0 (Ward identity, tested)
* DMFT susceptibility                   with the DMFT local self-energy in G this is the 'bubble
                                        approximation'; the full DMFT chi needs the local irreducible
                                        vertex (DGammaA / dual fermions, roadmap).
"""
from __future__ import annotations

import numpy as np

from .core.matsubara import MatsubaraGrid


def lattice_g(eps_k, sigma, grid: MatsubaraGrid, mu):
    sig = np.asarray(sigma)
    if sig.ndim == 1:
        sig = sig[None, None, :]
    return 1.0 / (grid.iw[None, None, :] + mu - eps_k[..., None] - sig)


def bubble(eps_k, sigma, grid: MatsubaraGrid, mu, hartree=0.0):
    G = lattice_g(eps_k, sigma, grid, mu)
    shift = eps_k - mu + hartree
    Gk_tau = grid.iw_to_tau(G, 1.0, shift[..., None])
    Gr_tau = np.fft.ifft2(Gk_tau, axes=(0, 1)).real
    nk = eps_k.shape[0]
    idx = (-np.arange(nk)) % nk
    Gm = Gr_tau[idx][:, idx][..., ::-1]
    chi_r = Gr_tau * Gm
    chi_q = np.fft.fft2(chi_r, axes=(0, 1)).real
    return grid.tau_to_iv(chi_q)


def rpa(chi0, U):
    spin = chi0 / (1 - U * chi0)
    charge = chi0 / (1 + U * chi0)
    return dict(spin=spin, charge=charge, stoner=float(np.max(U * chi0[..., 0].real)))
