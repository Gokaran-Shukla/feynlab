"""Maximum-entropy analytic continuation of fermionic Matsubara data
(Jarrell & Gubernatis, Phys. Rep. 269, 133 (1996)).

    G(iw_n) = int dw A(w)/(iw_n - w),   minimise  Q = chi^2/2 - alpha S[A],
    S = -int dw [A ln(A/m) - A + m].
A = m e^{x} keeps positivity; Q is minimised with L-BFGS (analytic gradient) on a
descending alpha sequence (warm starts). alpha is chosen by the classic criterion
-2 alpha S = N_eff or, if requested, by the 'historic' chi^2 = N_data.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize


def maxent(iw, G, sigma, omega, model=None, alphas=None, criterion: str = "classic"):
    iw = np.asarray(iw)
    w = np.asarray(omega)
    dw = np.gradient(w)
    m = np.ones_like(w) / (w[-1] - w[0]) if model is None else np.asarray(model, float) / np.sum(np.asarray(model) * dw)
    K = 1.0 / (iw[:, None] - w[None, :]) * dw[None, :]
    Kr = np.vstack([K.real, K.imag]) / sigma
    Gr = np.concatenate([G.real, G.imag]) / sigma
    alphas = np.logspace(3, -4, 36) if alphas is None else alphas
    x = np.zeros_like(w)
    out = []

    for a in alphas:
        def f(x):
            A = m * np.exp(x)
            r = Kr @ A - Gr
            S = -np.sum(dw * (A * x - A + m))
            Q = 0.5 * r @ r - a * S
            dQdA = Kr.T @ r + a * dw * x
            return Q, dQdA * A
        res = minimize(f, x, jac=True, method="L-BFGS-B", options=dict(maxiter=5000, gtol=1e-10, ftol=1e-14))
        x = res.x
        A = m * np.exp(x)
        chi2 = float(np.sum((Kr @ A - Gr) ** 2))
        S = float(-np.sum(dw * (A * x - A + m)))
        B = (np.sqrt(A)[:, None] * (Kr.T @ Kr) * np.sqrt(A)[None, :]) * np.sqrt(dw)[:, None] * np.sqrt(dw)[None, :] / dw.mean()
        lam = np.clip(np.linalg.eigvalsh(0.5 * (B + B.T)), 0, None)
        Neff = float(np.sum(lam / (a * dw.mean() + lam)))
        out.append((a, A.copy(), chi2, S, Neff))
    if criterion == "historic":
        N = len(Gr)
        i = int(np.argmin([abs(c - N) for _, _, c, _, _ in out]))
    else:
        i = int(np.argmin([abs(-2 * a * S - Ne) for a, _, _, S, Ne in out]))
    a, A, chi2, S, Ne = out[i]
    return dict(A=A, alpha=a, chi2=chi2, chi2_per_point=chi2 / len(Gr), omega=w, all=out)
