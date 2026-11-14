"""Generation, classification and analysis of Feynman diagrams.

Two diagram families are supported.

1. GenericDiagram  — fermions with a spin-independent, instantaneous two-body
   interaction v (Mattuck Chs. 4, 9, 10; electron gas, nuclear matter).
   Order-n Green's-function diagrams are bijections ('who connects to whom')
   between 2n+1 line *emitters*  {X, 0..2n-1} and 2n+1 *receivers* {Y, 0..2n-1};
   interaction line i joins vertices 2i and 2i+1.  X / Y are the external
   creation / annihilation points.

2. HubbardDiagram  — contact interaction U n_up n_dn (Hubbard / Anderson).
   Interaction line i has an up-end (2i) and a down-end (2i+1); fermion lines
   conserve spin, so diagrams are a pair of permutations.

For both families we provide: connectedness, proper (1PI), skeleton, tadpole
(Hartree-insertion) detection, number of closed fermion loops F, the Feynman
sign (-1)^F, canonical forms for removing duplicates, counting and drawing.
Because connected Green's-function diagrams with fixed external legs have no
non-trivial automorphisms, the number of topologically distinct diagrams equals
(#labelled diagrams)/(n! 2^n)  [generic]  or  (#labelled)/n!  [Hubbard].
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field

import numpy as np

X = -1  # external emitter  (creation point, incoming leg)
Y = -2  # external receiver (annihilation point, outgoing leg)


# ----------------------------------------------------------------------------
# small union-find
# ----------------------------------------------------------------------------
class _DSU:
    def __init__(self, items):
        self.p = {i: i for i in items}

    def find(self, a):
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb

    def n_components(self):
        return len({self.find(a) for a in self.p})


# ----------------------------------------------------------------------------
@dataclass
class Diagram:
    """A Green's-function / self-energy diagram.

    succ : dict  emitter -> receiver for every fermion line.  Internal vertex v
           both emits (outgoing line) and receives (incoming line).
    pairs: list of (a, b) internal vertices joined by interaction lines.
    spin : optional dict vertex -> 0/1 (Hubbard family).
    """
    order: int
    succ: dict
    pairs: list
    spin: dict | None = None
    family: str = "generic"
    _cache: dict = field(default_factory=dict, repr=False)

    # ---------------------------------------------------------------- basics
    @property
    def vertices(self):
        return list(range(2 * self.order))

    def fermion_lines(self):
        return [(a, b) for a, b in self.succ.items()]

    def open_path(self):
        """Vertices along the open fermion line X -> v1 -> ... -> Y."""
        path = []
        cur = self.succ[X]
        while cur != Y:
            path.append(cur)
            cur = self.succ[cur]
        return path

    def loops(self):
        """Closed fermion loops as lists of vertices."""
        on_path = set(self.open_path())
        seen = set(on_path)
        out = []
        for v in self.vertices:
            if v in seen:
                continue
            cyc = [v]
            seen.add(v)
            cur = self.succ[v]
            while cur != v:
                cyc.append(cur)
                seen.add(cur)
                cur = self.succ[cur]
            out.append(cyc)
        return out

    @property
    def n_loops(self):
        return len(self.loops())

    @property
    def sign(self):
        """Feynman sign (-1)^F, F = number of closed fermion loops."""
        return (-1) ** self.n_loops

    # ------------------------------------------------------- connectivity
    def _components(self, drop_fermion=(), drop_interaction=(), include_external=False):
        verts = list(self.vertices)
        if include_external:
            verts += [X, Y]
        d = _DSU(verts)
        for a, b in self.succ.items():
            if (a, b) in drop_fermion:
                continue
            if not include_external and (a == X or b == Y):
                continue
            d.union(a, b)
        for i, (a, b) in enumerate(self.pairs):
            if i in drop_interaction:
                continue
            d.union(a, b)
        return d

    def is_connected(self):
        return self._components(include_external=True).n_components() == 1

    def amputated_connected(self):
        """Connected after removing the two external legs (self-energy part)."""
        if self.succ[X] == Y:
            return False
        return self._components().n_components() == 1

    def is_proper(self):
        """One-particle irreducible: amputated diagram stays connected when any
        single internal fermion line is cut (Mattuck Sec. 10.2)."""
        if not self.amputated_connected():
            return False
        path = self.open_path()
        for a, b in zip(path[:-1], path[1:]):
            if self._components(drop_fermion={(a, b)}).n_components() > 1:
                return False
        return True

    def has_self_loop_tadpole(self):
        return any(self.succ[v] == v for v in self.vertices)

    def has_tadpole(self):
        """True if some interaction line is a bridge whose removal splits off a part
        not containing the external vertices (Hartree-type insertion)."""
        ext = {self.succ[X]}
        for v in self.vertices:
            if self.succ[v] == Y:
                ext.add(v)
        for i in range(len(self.pairs)):
            d = self._components(drop_interaction={i})
            if d.n_components() > 1:
                roots_ext = {d.find(e) for e in ext}
                a, b = self.pairs[i]
                if d.find(a) not in roots_ext or d.find(b) not in roots_ext:
                    return True
        return False

    def has_self_energy_insertion(self):
        """True if an internal fermion line carries a self-energy insertion, i.e.
        the diagram is NOT a skeleton (Mattuck Sec. 11.1)."""
        lines = [(a, b) for a, b in self.succ.items() if a != X and b != Y]
        first = self.succ[X]
        last = [v for v in self.vertices if self.succ[v] == Y][0]
        for (u1, v1), (u2, v2) in itertools.combinations(lines, 2):
            for (a1, b1), (a2, b2) in (((u1, v1), (u2, v2)), ((u2, v2), (u1, v1))):
                # insertion entered at b1 and left from a2
                d = self._components(drop_fermion={(a1, b1), (a2, b2)})
                if d.n_components() < 2:
                    continue
                c = d.find(b1)
                if d.find(a2) != c:
                    continue
                if d.find(a1) == c or d.find(b2) == c:
                    continue
                if d.find(first) == c or d.find(last) == c:
                    continue
                return True
        return False

    def is_skeleton(self):
        return self.is_proper() and not self.has_self_energy_insertion()

    # ---------------------------------------------------- classifications
    def is_hartree_fock_type(self):
        """Order-1 Hartree (tadpole) or Fock (exchange) diagram."""
        return self.order == 1 and self.amputated_connected()

    def describe(self):
        return dict(order=self.order, loops=self.n_loops, sign=self.sign,
                    proper=self.is_proper(), skeleton=self.is_skeleton() if self.is_proper() else False,
                    tadpole=self.has_tadpole(), path=self.open_path(), closed_loops=self.loops())

    # ------------------------------------------------------------ canonical
    def canonical_key(self):
        """Invariant key under relabelling of interaction lines (and end-swaps
        for the generic family)."""
        n = self.order
        best = None
        swaps = [(0,)] if self.family == "hubbard" else list(itertools.product((0, 1), repeat=n))
        for perm in itertools.permutations(range(n)):
            for sw in (itertools.product((0, 1), repeat=n) if self.family == "generic" else [tuple([0] * n)]):
                m = {X: X, Y: Y}
                for old, new in enumerate(perm):
                    a, b = 2 * old, 2 * old + 1
                    if sw[old]:
                        m[a], m[b] = 2 * new + 1, 2 * new
                    else:
                        m[a], m[b] = 2 * new, 2 * new + 1
                key = tuple(sorted((m[a], m[b]) for a, b in self.succ.items()))
                if best is None or key < best:
                    best = key
        return best

    # ------------------------------------------------------------ text form
    def to_string(self):
        """Compact text form: open path then loops, e.g. 'X-0-3-Y | (1 2)'."""
        p = "X-" + "-".join(map(str, self.open_path())) + "-Y"
        loops = " ".join("(" + " ".join(map(str, l)) + ")" for l in self.loops())
        inter = " ".join(f"[{a}~{b}]" for a, b in self.pairs)
        return f"{p} | {loops} | {inter}"


# ----------------------------------------------------------------------------
# generators
# ----------------------------------------------------------------------------
def _generic_labelled(n):
    emit = [X] + list(range(2 * n))
    recv = [Y] + list(range(2 * n))
    pairs = [(2 * i, 2 * i + 1) for i in range(n)]
    for perm in itertools.permutations(recv):
        succ = dict(zip(emit, perm))
        yield Diagram(n, succ, pairs, family="generic")


def _hubbard_labelled(n, ext_spin=0):
    """Spin-resolved labelled diagrams for G_up (ext_spin=0)."""
    pairs = [(2 * i, 2 * i + 1) for i in range(n)]
    up = [2 * i for i in range(n)]
    dn = [2 * i + 1 for i in range(n)]
    if ext_spin == 1:
        up, dn = dn, up
    for pu in itertools.permutations([Y] + up):
        for pd in itertools.permutations(dn):
            succ = dict(zip([X] + up, pu))
            succ.update(dict(zip(dn, pd)))
            spin = {v: (0 if v % 2 == 0 else 1) for v in range(2 * n)}
            yield Diagram(n, succ, pairs, spin=spin, family="hubbard")


def count_diagrams(n: int, family: str = "generic", kind: str = "G", exclude_tadpoles=False,
                   by_loops: bool = False):
    """Number of topologically distinct diagrams.

    kind: 'G' (connected Green's function), 'proper' (self-energy, 1PI),
          'improper' (connected, amputated-connected but reducible), 'skeleton'.
    """
    gen = _generic_labelled(n) if family == "generic" else _hubbard_labelled(n)
    norm = math.factorial(n) * (2**n if family == "generic" else 1)
    counts = {}
    total = 0
    for d in gen:
        if not d.is_connected() or not d.amputated_connected():
            continue
        if kind in ("proper", "skeleton") and not d.is_proper():
            continue
        if kind == "improper" and d.is_proper():
            continue
        if kind == "skeleton" and d.has_self_energy_insertion():
            continue
        if exclude_tadpoles and d.has_tadpole():
            continue
        total += 1
        if by_loops:
            counts[d.n_loops] = counts.get(d.n_loops, 0) + 1
    if by_loops:
        return {k: v // norm for k, v in sorted(counts.items())}
    assert total % norm == 0, "non-free group action detected"
    return total // norm


def unique_diagrams(n: int, family: str = "hubbard", kind: str = "proper", exclude_tadpoles=True):
    """List of topologically distinct diagrams (one representative each)."""
    gen = _generic_labelled(n) if family == "generic" else _hubbard_labelled(n)
    seen = {}
    for d in gen:
        if not d.is_connected() or not d.amputated_connected():
            continue
        if kind in ("proper", "skeleton") and not d.is_proper():
            continue
        if kind == "skeleton" and d.has_self_energy_insertion():
            continue
        if exclude_tadpoles and d.has_tadpole():
            continue
        k = d.canonical_key()
        if k not in seen:
            seen[k] = d
    return list(seen.values())
