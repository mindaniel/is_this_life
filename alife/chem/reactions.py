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

CARB, CA, CAT, TAIL, STAIL = INDEX["carb"], INDEX["Ca"], INDEX["cat"], INDEX["tail"], INDEX["stail"]


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
        """One round of joins and splits. Frozen carbons (dry crust) don't react,
        nor carbons bonded to other atoms (glycine, glycolonitrile: not sugar)."""
        carb = np.flatnonzero((world.type == CARB) & ~world.frozen & pure_carb(world))
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
        cb = np.flatnonzero((world.type[b[:, 0]] == CARB) & (world.type[b[:, 1]] == CARB) & ~world.frozen[b[:, 0]]
                            & pure_carb(world)[b[:, 0]])
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


def pure_carb(world):
    """True for beads in molecules made only of carb (sugars and formaldehyde)."""
    other = np.bincount(world.mol[world.type != CARB], minlength=world.n_mol)
    return other[world.mol] == 0


def chain_lengths(world):
    """Number of sugar molecules of each length: counts[n] = chains with n carbons."""
    carb = (world.type == CARB) & pure_carb(world)
    if not carb.any():
        return np.zeros(1, dtype=int)
    size = np.bincount(world.mol[carb])
    return np.bincount(size[size > 0])


class LipidJoin:
    """A catalyst joins two fatty acids into one two-chain lipid (phospholipid-like).

    Real cells make phospholipids with enzymes: two fatty-acid chains on one
    glycerol-phosphate head. Lab protocells (Budin & Szostak 2011; Adamala &
    Szostak 2013 with a catalyst trapped inside) show why it matters: two-chain
    lipids almost never leave a membrane, so a vesicle that makes them keeps its
    fatty acids and grows at its neighbours' expense. Nothing here says "grow".

      join   two single fatty acids (head + one tail chain) whose heads touch
             (within `reach`) bond head to head, with chance p_join if a "cat"
             bead is within `cat_reach` of the pair, p_join * no_cat if not.
      break  each head-head bond breaks with chance p_break (slow hydrolysis).

    `chains`: which tail types count as a joinable chain (STAIL alone = only the
    fatty acids; the long-tail lipids already stand for finished membrane lipids).
    `length`: rest length of the new head-head bond. The default bond (1.0) pulls
    the two heads to ~0.89, closer than neighbours in a leaflet (~1.15): each join
    shrinks the leaflet a little (a model artefact; a real two-chain lipid takes
    about the room of its two fatty acids).

    Shortcut (flagged): the join uses no fuel. In real chemistry it needs an
    activated building block (glycerol phosphate, ATP-like energy); the catalyst
    only speeds it up. Next step: make a food bead that is used up.
    """

    def __init__(self, every=50, reach=1.3, cat_reach=2.5, p_join=0.2, no_cat=0.0, p_break=0.0,
                 chains=(TAIL, STAIL), length=None):
        self.every, self.reach, self.cat_reach = every, reach, cat_reach
        self.chains, self.length = chains, length
        self.p_join, self.no_cat, self.p_break = p_join, no_cat, p_break
        self.joins = self.breaks = 0

    def single_heads(self, world):
        """Head beads of lone fatty acids: molecules of one head + one tail chain of a `chains` type."""
        heads = np.isin(world.type, (INDEX["head"], INDEX["fa_head"]))
        n_head = np.bincount(world.mol[heads], minlength=world.n_mol)
        n_tail = np.bincount(world.mol[np.isin(world.type, self.chains)], minlength=world.n_mol)
        size = np.bincount(world.mol, minlength=world.n_mol)
        lone = (n_head == 1) & (n_tail >= 1) & (size == n_head + n_tail)
        return np.flatnonzero(heads & lone[world.mol])

    def head_bonds(self, world):
        """Indices of bonds that join two heads (the two-chain lipids)."""
        heads = np.isin(world.type, (INDEX["head"], INDEX["fa_head"]))
        b = world.bonds
        return np.flatnonzero(heads[b[:, 0]] & heads[b[:, 1]]) if len(b) else np.zeros(0, int)

    def react(self, world):
        self._break(world)
        self._join(world)

    def _join(self, world):
        heads = self.single_heads(world)
        if len(heads) < 2:
            return
        p = world._wrapped()
        pairs = cKDTree(p[heads], boxsize=world.L).query_pairs(self.reach, output_type="ndarray")
        if not len(pairs):
            return
        i, j = heads[pairs[:, 0]], heads[pairs[:, 1]]
        mid = (p[i] + 0.5 * world.delta(i, j)) % world.L
        chance = np.full(len(i), self.p_join * self.no_cat)
        cat = np.flatnonzero(world.type == CAT)
        if len(cat):
            near = cKDTree(p[cat], boxsize=world.L).query(mid, distance_upper_bound=self.cat_reach)[0] < np.inf
            chance[near] = self.p_join
        hit = world.rng.uniform(size=len(i)) < chance
        used = set()
        for a, b in zip(i[hit], j[hit]):
            if a in used or b in used:
                continue
            used.update((a, b))
            world.add_bond(a, b, self.length)
            self.joins += 1

    def _break(self, world):
        if not self.p_break:
            return
        k = self.head_bonds(world)
        hit = k[world.rng.uniform(size=len(k)) < self.p_break]
        for b in sorted(hit, reverse=True):
            world.break_bond(b)
            self.breaks += 1
