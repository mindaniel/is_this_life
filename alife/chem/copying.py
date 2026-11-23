"""RNA-like copying chemistry. A monomer is two beads, like a real nucleotide: a backbone bead (nbb,
sugar-phosphate) and a base bead (nx or ny) bonded to it. Strands join backbone to backbone; bases
hang off the side.

Nothing here says "copy". There are three plain rules plus ordinary physics:

1. Pairing (base pairing). An X base and a Y base close together can stick (energy eps_match, hydrogen
   bonds). Mismatched X-X / Y-Y stick more weakly (eps_mis): that is where copying errors come from.
   Each base holds at most ONE partner: real bases pair by directional hydrogen bonds, so one base faces
   one other. Beads are round, so this one-partner limit is put in by hand (flagged shortcut).
   Pairs form and break by Monte Carlo with Boltzmann odds, so warmer water melts them.
2. Joining (ligation). Two backbone ends that touch, with their bases side by side (lined up), can join with a small chance per check, whatever their
   bases are doing. Real monomers must be "activated" to join (e.g. imidazole-activated nucleotides);
   here that energy is assumed present (flagged; later it should come from a fuel).
3. Breaking (hydrolysis). Any backbone bond can break with a small chance.
Geometry (ordinary molecular stiffness, like the lipids' bends): within a strand, neighbouring bases stack
(a spring base-base, length 1), so the bases line up on one side; a base pair is flat and rigid, so each
backbone sits in line behind its paired base (spring backbone to partner base, length r_pair + 1).

Copying, if it happens, comes out of these: a strand's bases hold incoming monomers side by side, their
backbones touch long enough to join, and the new strand reads the template's order.

Simplifications (flagged): strands have no direction (real RNA is 5'->3' and pairs antiparallel); a strand
can't pair with itself (no hairpins); backbone charge (phosphate, -1) left out.
"""
import numpy as np
from scipy.spatial import cKDTree

from .elements import INDEX

NX, NY, BB = INDEX["nx"], INDEX["ny"], INDEX["nbb"]
BASES = (NX, NY)


def base_of(world):
    """base_of[i] = base bead hanging off backbone bead i (-1 if none); backbone_of the other way."""
    b, t = world.bonds, world.type
    bo, bk = np.full(world.n, -1), np.full(world.n, -1)
    for i, j in ((0, 1), (1, 0)):
        m = (t[b[:, i]] == BB) & np.isin(t[b[:, j]], BASES)
        bo[b[m, i]], bk[b[m, j]] = b[m, j], b[m, i]
    return bo, bk


def bb_bonds(world):
    t = world.type
    return np.flatnonzero((t[world.bonds[:, 0]] == BB) & (t[world.bonds[:, 1]] == BB))


def chain_order(world, m):
    """Backbone beads of molecule m in order (one end to the other)."""
    beads = np.flatnonzero((world.mol == m) & (world.type == BB))
    if len(beads) == 1:
        return list(beads)
    b = world.bonds[bb_bonds(world)]
    b = b[world.mol[b[:, 0]] == m]
    nbr = {i: [] for i in beads}
    for i, j in b:
        nbr[i].append(j)
        nbr[j].append(i)
    ends = [i for i in beads if len(nbr[i]) == 1]
    if not ends:
        return list(beads)
    out, prev, cur = [ends[0]], -1, ends[0]
    while True:
        nxt = [k for k in nbr[cur] if k != prev]
        if not nxt:
            return out
        prev, cur = cur, nxt[0]
        out.append(cur)


def sequence(world, m, bo=None):
    bo = base_of(world)[0] if bo is None else bo
    return "".join("X" if world.type[bo[i]] == NX else "Y" for i in chain_order(world, m))


def add_monomer(world, x, d, base="nx"):
    """Backbone at x, base one bond length along unit vector d."""
    world.add_molecule([x, x + d], ["nbb", base], bonds=[(0, 1)])


def add_strand(world, start, along, side, seq):
    """A straight strand: backbone beads 1 apart along `along`, bases one step to `side` (unit vectors)."""
    k = len(seq)
    bb = start + np.outer(np.arange(k), along)
    pos = np.vstack([bb, bb + side])
    types = ["nbb"] * k + ["nx" if c == "X" else "ny" for c in seq]
    bonds = [(i, i + 1) for i in range(k - 1)] + [(i, k + i) for i in range(k)]
    world.add_molecule(pos, types, bonds=bonds, bends=[(i, i + 2) for i in range(k - 2)])
    n0 = world.n - 2 * k
    for i in range(k - 1):   # stacking between neighbouring bases
        world.add_bond(n0 + k + i, n0 + k + i + 1, 1.0)


class Copying:
    def __init__(self, world, every=50, eps_match=4.0, eps_mis=1.5, reach=1.4, k_pair=30.0, r_pair=1.1,
                 p_pair=0.5, p_join=0.01, join_reach=1.3, stack_reach=1.6, p_break=0.0):
        self.every, self.eps_match, self.eps_mis, self.reach = every, eps_match, eps_mis, reach
        self.k_pair, self.r_pair, self.p_pair = k_pair, r_pair, p_pair
        self.p_join, self.join_reach, self.stack_reach, self.p_break = p_join, join_reach, stack_reach, p_break
        self.partner = np.full(world.n, -1)
        self.joins = self.templated = self.templated_ok = self.breaks = 0
        self.world = world
        self._bk = None   # backbone bead of each base (bonds base-backbone never change)
        world.repel_inside = tuple(set(world.repel_inside) | {NX, NY, BB})   # floppy strands must not fold through
        world.external = self.force

    # ---- pairing: a spring between partner bases, made and broken by Monte Carlo ----------
    def _pairs(self):
        i = np.flatnonzero(self.partner >= 0)
        i = i[i < self.partner[i]]
        return i, self.partner[i]

    def _eps(self, i, j):
        t = self.world.type
        return np.where(t[i] != t[j], self.eps_match, self.eps_mis)

    def _spring(self, i, j):
        """Stretch energy of the pair's springs (base-base and each backbone to the other base)."""
        w = self.world
        if self._bk is None or len(self._bk) != w.n:
            self._bk = base_of(w)[1]
        e = np.zeros(len(i))
        for a, b, r0 in ((i, j, self.r_pair), (self._bk[i], j, self.r_pair + 1.0), (self._bk[j], i, self.r_pair + 1.0)):
            ok = a >= 0
            d = w.delta(a[ok], b[ok])
            e[ok] += 0.5 * self.k_pair * (np.sqrt(np.einsum("ij,ij->i", d, d)) - r0) ** 2
        return e

    def force(self, pos):
        F = np.zeros_like(pos)
        i, j = self._pairs()
        if not len(i):
            return F, 0.0
        w = self.world
        if self._bk is None or len(self._bk) != w.n:
            self._bk = base_of(w)[1]
        bi, bj = self._bk[i], self._bk[j]
        e = -np.sum(self._eps(i, j))
        for a, b, r0 in ((i, j, self.r_pair), (bi, j, self.r_pair + 1.0), (bj, i, self.r_pair + 1.0)):
            ok = a >= 0
            a, b = a[ok], b[ok]
            d = w.delta(a, b)
            r = np.sqrt(np.einsum("ij,ij->i", d, d))
            w._scatter(F, a, b, (-self.k_pair * (r - r0) / np.maximum(r, 1e-9))[:, None] * d)
            e += float(np.sum(0.5 * self.k_pair * (r - r0) ** 2))
        return F, e

    def pair_step(self):
        w = self.world
        kT = float(np.mean(w.kT))
        if len(self.partner) < w.n:
            self.partner = np.r_[self.partner, np.full(w.n - len(self.partner), -1)]
        i, j = self._pairs()   # break: Metropolis on losing the pair energy (minus the spring's stretch)
        if len(i):
            dE = self._eps(i, j) - self._spring(i, j)
            go = w.rng.uniform(size=len(i)) < self.p_pair * np.minimum(1.0, np.exp(-dE / kT))
            self.partner[i[go]] = self.partner[j[go]] = -1
        free = np.flatnonzero(np.isin(w.type, BASES) & (self.partner < 0) & ~w.frozen)
        if len(free) < 2:
            return
        p = w._wrapped()
        pr = cKDTree(p[free], boxsize=w.L).query_pairs(self.reach, output_type="ndarray")
        if not len(pr):
            return
        a, b = free[pr[:, 0]], free[pr[:, 1]]
        keep = w.mol[a] != w.mol[b]
        a, b = a[keep], b[keep]
        order = w.rng.permutation(len(a))
        a, b = a[order], b[order]
        dE = -self._eps(a, b) + self._spring(a, b)
        go = w.rng.uniform(size=len(a)) < self.p_pair * np.minimum(1.0, np.exp(-dE / kT))
        for x, y in zip(a[go], b[go]):
            if self.partner[x] < 0 and self.partner[y] < 0:
                self.partner[x], self.partner[y] = y, x

    # ---- joining and breaking the backbone ----------------------------------------------
    def join_step(self):
        w = self.world
        nb = w.bonds[bb_bonds(w)]
        deg = np.bincount(nb.ravel(), minlength=w.n)
        ends = np.flatnonzero((w.type == BB) & (deg < 2) & ~w.frozen)
        if len(ends) < 2:
            return
        p = w._wrapped()
        pr = cKDTree(p[ends], boxsize=w.L).query_pairs(self.join_reach, output_type="ndarray")
        if not len(pr):
            return
        a, b = ends[pr[:, 0]], ends[pr[:, 1]]
        keep = w.mol[a] != w.mol[b]
        a, b = a[keep], b[keep]
        hit = w.rng.uniform(size=len(a)) < self.p_join
        if not hit.any():
            return
        bo, bk = base_of(w)
        # only monomers lined up side by side can join (their bases must be close enough to stack)
        d = w.delta(bo[a], bo[b])
        hit &= np.einsum("ij,ij->i", d, d) < self.stack_reach ** 2
        bonded = set(map(tuple, np.sort(nb, axis=1)))
        used = set()
        t = w.type
        for x, y in zip(a[hit], b[hit]):
            if x in used or y in used or w.mol[x] == w.mol[y]:
                continue
            used.update((x, y))
            nx_ = [k for k in nb[(nb == x).any(axis=1)].ravel() if k != x]
            ny_ = [k for k in nb[(nb == y).any(axis=1)].ravel() if k != y]
            px, py = self.partner[bo[x]], self.partner[bo[y]]
            if px >= 0 and py >= 0:
                u, v = bk[px], bk[py]
                if (min(u, v), max(u, v)) in bonded:   # both bases held by neighbouring bases of one strand
                    self.templated += 1
                    self.templated_ok += int(t[bo[x]] != t[px] and t[bo[y]] != t[py])
            w.add_bond(x, y)
            w.add_bond(bo[x], bo[y], 1.0)   # the two bases now stack
            new = [(k, y) for k in nx_] + [(x, k) for k in ny_]   # bends keep the backbone straight
            if new:
                w.bends = np.vstack([w.bends, new])
            self.joins += 1

    def break_step(self):
        if self.p_break <= 0:
            return
        w = self.world
        idx = bb_bonds(w)
        hit = idx[w.rng.uniform(size=len(idx)) < self.p_break]
        if not len(hit):
            return
        bo = base_of(w)[0]
        pairs = [tuple(w.bonds[k]) for k in hit]
        for x, y in pairs:   # backbone bond and the stacking bond of its two bases (find by beads: indices shift)
            for a, b in ((x, y), (bo[x], bo[y])):
                k = np.flatnonzero(((w.bonds[:, 0] == a) & (w.bonds[:, 1] == b)) | ((w.bonds[:, 0] == b) & (w.bonds[:, 1] == a)))
                if len(k):
                    w.break_bond(k[0])
            self.breaks += 1
        if len(hit):   # bends across a broken bond go too
            w.bends = w.bends[w.mol[w.bends[:, 0]] == w.mol[w.bends[:, 1]]]

    def react(self):
        self.pair_step()
        self.join_step()
        self.break_step()
