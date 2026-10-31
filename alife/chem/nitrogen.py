"""Nitrogen chemistry: HCN and ammonia, the first step towards pigments.

Beads: "hcn" (HCN, rains in from lightning / impacts), "nh3" (ammonia),
"carb" (formaldehyde, as in reactions.py). Water is implicit, so the H2O
a reaction uses or releases is not tracked.

  strecker   formaldehyde + HCN touching, with ammonia touching too:
             they bond into N-C-C = aminonitrile, which water turns into
             glycine, the simplest amino acid. (The NH3 the second step
             gives back is not tracked.)
  cyano      formaldehyde + HCN with no ammonia: glycolonitrile, a side
             product that goes nowhere here. Ammonia decides which wins.
  polymer    two HCN chain ends touching can bond. Much faster with
             ammonia nearby (a base, real HCN polymer forms at pH ~9).
             HCN + HCN (the dimer) is the slow first step; growing a
             chain is fast. Long HCN chains are the dark "HCN polymer",
             which absorbs light: a crude first pigment.
  hydrolysis a lone HCN reacts with water (-> formamide -> formic acid +
             NH3). The bead becomes NH3; the carbon leaves as formic acid,
             not tracked.

Hydrolysis needs one molecule, polymer needs two to meet. So dilute HCN
should mostly be lost to water, concentrated HCN should polymerise
(why chemists think HCN had to be concentrated: evaporation, freezing).
Experiment e14 checks that.

Rates are chances per check at kT0, scaled by exp(-Ea (1/kT - 1/kT0)).
"""
import numpy as np
from scipy.spatial import cKDTree

from .elements import INDEX

CARB, HCN, NH3 = INDEX["carb"], INDEX["hcn"], INDEX["nh3"]


class Nitrogen:
    def __init__(self, every=50, reach=1.3, base_reach=2.0, kT0=0.6,
                 p_strecker=0.05, p_cyano=0.008, p_poly=0.05, no_base=0.05, dimer_factor=0.01,
                 p_hydro=2e-4, Ea=4.0, Ea_hydro=6.0):
        self.every, self.reach, self.base_reach, self.kT0 = every, reach, base_reach, kT0
        self.p_strecker, self.p_cyano = p_strecker, p_cyano
        self.p_poly, self.no_base, self.dimer_factor = p_poly, no_base, dimer_factor
        self.p_hydro, self.Ea, self.Ea_hydro = p_hydro, Ea, Ea_hydro
        self.glycine = self.cyano = self.poly = self.hydrolysed = 0

    def _heat(self, world, idx, Ea):
        kT = np.full(len(idx), float(world.kT)) if world.kT_field is None else \
            np.maximum(world.kT_field(world.pos[idx]), 1e-3)
        return np.exp(-Ea * (1.0 / kT - 1.0 / self.kT0))

    def react(self, world):
        """One round. Frozen beads (dry crust) don't react."""
        size = np.bincount(world.mol, minlength=world.n_mol)
        live = ~world.frozen
        lone = live & (size[world.mol] == 1)
        p = world._wrapped()
        self._strecker(world, p, lone)
        size = np.bincount(world.mol, minlength=world.n_mol)
        self._polymer(world, p, size, live)
        self._hydrolysis(world, live & (size[world.mol] == 1))

    def _near(self, p, world, idx, points, r):
        """True for each point that has a bead of `idx` within r."""
        if not len(idx):
            return np.zeros(len(points), dtype=bool)
        return cKDTree(p[idx], boxsize=world.L).query(points, distance_upper_bound=r)[0] < np.inf

    def _strecker(self, world, p, lone):
        carb = np.flatnonzero(lone & (world.type == CARB))
        hcn = np.flatnonzero(lone & (world.type == HCN))
        if not len(carb) or not len(hcn):
            return
        hits = cKDTree(p[hcn], boxsize=world.L).query_ball_point(p[carb], self.reach)
        nh3 = np.flatnonzero(lone & (world.type == NH3))
        ntree = cKDTree(p[nh3], boxsize=world.L) if len(nh3) else None
        used = set()
        for c, near in zip(carb, hits):
            if not near:
                continue
            h = hcn[world.rng.choice(near)]
            if h in used:
                continue
            amm = [nh3[k] for k in (ntree.query_ball_point(p[c], self.reach) if ntree else []) if nh3[k] not in used]
            heat = self._heat(world, [c], self.Ea)[0]
            if amm:
                if world.rng.uniform() < self.p_strecker * heat:
                    n = amm[world.rng.integers(len(amm))]
                    world.add_bond(n, c)
                    world.add_bond(c, h)
                    used.update((c, h, n))
                    self.glycine += 1
            elif world.rng.uniform() < self.p_cyano * heat:
                world.add_bond(c, h)
                used.update((c, h))
                self.cyano += 1

    def _polymer(self, world, p, size, live):
        hcn_mol = np.bincount(world.mol[world.type == HCN], minlength=world.n_mol)
        pure = live & (world.type == HCN) & (hcn_mol[world.mol] == size[world.mol])   # HCN-only molecules
        hb = world.bonds[(world.type[world.bonds[:, 0]] == HCN) & (world.type[world.bonds[:, 1]] == HCN)]
        degree = np.bincount(hb.ravel(), minlength=world.n)
        ends = np.flatnonzero(pure & (degree < 2))
        if len(ends) < 2:
            return
        pairs = cKDTree(p[ends], boxsize=world.L).query_pairs(self.reach, output_type="ndarray")
        if not len(pairs):
            return
        i, j = ends[pairs[:, 0]], ends[pairs[:, 1]]
        keep = world.mol[i] != world.mol[j]   # no rings
        i, j = i[keep], j[keep]
        if not len(i):
            return
        mid = (p[i] + 0.5 * world.delta(i, j)) % world.L
        chance = np.full(len(i), self.p_poly)
        chance[~self._near(p, world, np.flatnonzero(world.type == NH3), mid, self.base_reach)] *= self.no_base
        chance[(size[world.mol[i]] == 1) & (size[world.mol[j]] == 1)] *= self.dimer_factor
        chance *= self._heat(world, i, self.Ea)
        hit = world.rng.uniform(size=len(i)) < chance
        used = set()
        for a, b in zip(i[hit], j[hit]):
            if a in used or b in used or world.mol[a] == world.mol[b]:
                continue
            used.update((a, b))
            world.add_bond(a, b)
            self.poly += 1

    def _hydrolysis(self, world, lone):
        hcn = np.flatnonzero(lone & (world.type == HCN))
        if not len(hcn):
            return
        hit = hcn[world.rng.uniform(size=len(hcn)) < self.p_hydro * self._heat(world, hcn, self.Ea_hydro)]
        world.type[hit] = NH3   # same charge (0), so neighbour lists stay valid
        self.hydrolysed += len(hit)


def census(world):
    """What the nitrogen ended up in: lone HCN / NH3 / formaldehyde, glycine,
    glycolonitrile, and HCN chains by length (hcn_chains[n] = chains of n HCN)."""
    t, size = world.type, np.bincount(world.mol, minlength=world.n_mol)
    per = lambda k: np.bincount(world.mol[t == k], minlength=world.n_mol)
    n_c, n_h, n_n = per(CARB), per(HCN), per(NH3)
    mols = np.unique(world.mol)
    s, c, h, n = size[mols], n_c[mols], n_h[mols], n_n[mols]
    pure_h = (h == s)
    return dict(hcn=int(((s == 1) & (h == 1)).sum()), nh3=int(((s == 1) & (n == 1)).sum()),
                carb=int(((s == 1) & (c == 1)).sum()),
                glycine=int(((s == 3) & (c == 1) & (h == 1) & (n == 1)).sum()),
                cyano=int(((s == 2) & (c == 1) & (h == 1)).sum()),
                hcn_chains=np.bincount(h[pure_h & (s > 1)], minlength=2))
