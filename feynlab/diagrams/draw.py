"""Draw Feynman diagrams (matplotlib): the open fermion line runs left->right along
the axis, closed loops are drawn above it, interaction lines are dashed/wavy."""
from __future__ import annotations

import numpy as np

from .enumerate import X, Y, Diagram


def _wavy(ax, p, q, color, amp=0.06, n=6):
    p, q = np.asarray(p, float), np.asarray(q, float)
    t = np.linspace(0, 1, 200)
    d = q - p
    L = np.hypot(*d) + 1e-12
    nrm = np.array([-d[1], d[0]]) / L
    pts = p[None] + t[:, None] * d[None] + amp * np.sin(2 * np.pi * n * t)[:, None] * nrm[None]
    ax.plot(pts[:, 0], pts[:, 1], color=color, lw=1.2)


def _arrow_arc(ax, p, q, color, bend=0.0):
    from matplotlib.patches import FancyArrowPatch
    ax.add_patch(FancyArrowPatch(p, q, connectionstyle=f"arc3,rad={bend}", arrowstyle="-|>",
                                 mutation_scale=10, color=color, lw=1.4, shrinkA=4, shrinkB=4))


def draw_diagram(d: Diagram, ax=None, title: str | None = None):
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots(figsize=(4, 2.4))
    path = d.open_path()
    pos = {X: (-1.0, 0.0)}
    for i, v in enumerate(path):
        pos[v] = (float(i), 0.0)
    pos[Y] = (float(len(path)), 0.0)
    # place loop vertices on circles above
    for j, loop in enumerate(d.loops()):
        cx = (len(path) - 1) / 2 + 1.4 * (j - (len(d.loops()) - 1) / 2)
        r = 0.45
        for k, v in enumerate(loop):
            a = np.pi / 2 + 2 * np.pi * k / len(loop)
            pos[v] = (cx + r * np.cos(a), 1.3 + 0.6 * j + r * np.sin(a))
    fcol, icol = "#2a78d6", "#eb6834"
    for a, b in d.succ.items():
        if a == b:
            ax.add_patch(plt.Circle((pos[a][0], pos[a][1] + 0.25), 0.25, fill=False, color=fcol, lw=1.4))
            continue
        bend = 0.0 if (pos[a][1] == 0 and pos[b][1] == 0 and abs(pos[a][0] - pos[b][0]) <= 1.01) else 0.35
        _arrow_arc(ax, pos[a], pos[b], fcol, bend)
    for a, b in d.pairs:
        _wavy(ax, pos[a], pos[b], icol)
    for v, (x, y) in pos.items():
        if v in (X, Y):
            continue
        ax.plot(x, y, "o", ms=4, color="#1f1f1e")
        lab = str(v) if d.spin is None else f"{v}{'↑' if d.spin[v] == 0 else '↓'}"
        ax.annotate(lab, (x, y), xytext=(0, -11), textcoords="offset points", ha="center", fontsize=7,
                    color="#6b6a64")
    ax.set_aspect("equal")
    ax.axis("off")
    info = d.describe()
    ax.set_title(title or f"order {d.order}, loops F={info['loops']}, sign {'+' if d.sign > 0 else '-'}",
                 fontsize=8, loc="left")
    return ax


def draw_all(diagrams, ncol: int = 4, fname: str | None = None):
    import matplotlib.pyplot as plt
    n = len(diagrams)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 2.4 * nrow), dpi=130)
    axes = np.atleast_1d(axes).ravel()
    for ax, d in zip(axes, diagrams):
        draw_diagram(d, ax)
    for ax in axes[n:]:
        ax.axis("off")
    fig.tight_layout()
    if fname:
        fig.savefig(fname)
    return fig
