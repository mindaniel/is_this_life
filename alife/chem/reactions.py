"""Chemical reactions: bonds that form and break while the world runs.

Formose (sugar from formaldehyde), simplified. Every carbon is one "carb"
bead (CH2O). A lone carb is formaldehyde; a chain of 2 is glycolaldehyde,
5 is ribose-sized.

  join    two chain ends that touch can bond. Faster when hot (Arrhenius),
          much faster with a Ca2+ ion nearby (catalyst).
          C1 + C1 is very slow (the hard first step of real formose);
          C1 + a chain of 2 or more is fast.
  split   a chain of 4 or more can break into two pieces of at least 2
          (retro-aldol: C4 -> C2 + C2). Faster when hot.

Nothing here codes the formose cycle itself. Grow (C2 -> C4) + split
(C4 -> 2 C2) together should make sugar multiply (autocatalysis), with a
slow start until the first C2 forms by chance. Experiment e06 checks that.

Rates are chances per check (every `every` steps) at temperature kT0,
scaled by exp(-Ea (1/kT - 1/kT0)) at the bead's local temperature.
"""
import numpy as np
from scipy.spatial import cKDTree

from .elements import INDEX

CARB, CA = INDEX["carb"], INDEX["Ca"]


class Formose:
    def __init__(self, every=50, reach=1.3, cat_reach=2.0, kT0=0.6,
                 p_join=0.5, no_cat=0.01, seed_factor=0.0001, Ea_join=4.0,
                 p_split=0.02, Ea_split=8.0, min_split=4):
        self.every, self.reach, self.cat_reach, self.kT0 = every, reach, cat_reach, kT0
        self.p_join, self.no_cat, self.seed_factor, self.Ea_join = p_join, no_cat, seed_factor, Ea_join
        self.p_split, self.Ea_split, self.min_split = p_split, Ea_split, min_split
        self.joins = self.splits = 0

    def _kT(self, world, idx):
        if world.kT_field is None:
            return np.full(len(idx), float(world.kT))
        return np.maximum(world.kT_field(world.pos[idx]), 1e-3)

    def _heat(self, kT, Ea):
        return np.exp(-Ea * (1.0 / kT - 1.0 / self.kT0))

    def react(self, world):
        """One round of joins and splits."""
        carb = np.flatnonzero(world.type == CARB)
        if len(carb) < 2:
            return
        size = np.bincount(world.mol[carb], minlength=world.n_mol)   # chain length per molecule
        cb = world.bonds[(world.type[world.bonds[:, 0]] == CARB) & (world.type[world.bonds[:, 1]] == CARB)]
        degree = np.bincount(cb.ravel(), minlength=world.n)
        self._join(world, carb, size, degree)
        self._split(world)

    def _join(self, world, carb, size, degree):
        ends = carb[degree[carb] < 2]
        if len(ends) < 2:
            return
        p = world._wrapped()
        pairs = cKDTree(p[ends], boxsize=world.L).query_pairs(self.reach, output_type="ndarray")
        if not len(pairs):
            return
        i, j = ends[pairs[:, 0]], ends[pairs[:, 1]]
        keep = world.mol[i] != world.mol[j]   # no rings
        i, j = i[keep], j[keep]
        if not len(i):
            return
        mid = (p[i] + 0.5 * world.delta(i, j)) % world.L
        chance = np.full(len(i), self.p_join)
        ca = np.flatnonzero(world.type == CA)
        has_cat = np.zeros(len(i), dtype=bool)
        if len(ca):
            has_cat = cKDTree(p[ca], boxsize=world.L).query(mid, distance_upper_bound=self.cat_reach)[0] < np.inf
        chance[~has_cat] *= self.no_cat
        both_single = (size[world.mol[i]] == 1) & (size[world.mol[j]] == 1)
        chance[both_single] *= self.seed_factor
        chance *= self._heat(self._kT(world, i), self.Ea_join)
        hit = world.rng.uniform(size=len(i)) < chance
        used = set()
        for a, b in zip(i[hit], j[hit]):
            if a in used or b in used or world.mol[a] == world.mol[b]:
                continue
            used.update((a, b))
            world.add_bond(a, b)
            self.joins += 1

    def _split(self, world):
        b = world.bonds
        cb = np.flatnonzero((world.type[b[:, 0]] == CARB) & (world.type[b[:, 1]] == CARB))
        if not len(cb):
            return
        size = np.bincount(world.mol, minlength=world.n_mol)
        long_ = cb[size[world.mol[b[cb, 0]]] >= self.min_split]
        if not len(long_):
            return
        chance = self.p_split * self._heat(self._kT(world, b[long_, 0]), self.Ea_split)
        hit = long_[world.rng.uniform(size=len(long_)) < chance]
        done = set()
        for k in sorted(hit, reverse=True):   # delete from the end so indices stay valid
            i, j = world.bonds[k]
            m = world.mol[i]
            if m in done:
                continue
            # both pieces must keep at least 2 carbons (no C1 made by splitting)
            mine = world.mol == m
            piece = world._reachable(j, mine & ~np.isin(np.arange(world.n), [i]))
            n_j = int(piece.sum())
            if n_j < 2 or mine.sum() - n_j < 2:   # count now: earlier splits change sizes
                continue
            done.add(m)
            world.break_bond(k)
            self.splits += 1


def chain_lengths(world):
    """Number of carb molecules of each length: counts[n] = chains with n carbons."""
    carb = world.type == CARB
    if not carb.any():
        return np.zeros(1, dtype=int)
    size = np.bincount(world.mol[carb])
    return np.bincount(size[size > 0])
