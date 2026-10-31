"""Leaky membranes: small molecules slip through lipid walls, big ones don't.

In the bead model a solute can't squeeze between packed lipid tails, so
nothing ever crosses a membrane. Real membranes leak, and how fast depends
on the molecule:

  formaldehyde (C1)     very fast   small, dissolves in the oily middle
  glycolaldehyde (C2)   fast
  C3-C4 sugars          slow
  C5+ sugars (ribose)   ~never      too big and too polar
  ions (Na, Ca, Cl...)  very slow   charged: the oily middle hates them

Mechanism (Monte Carlo hop): every `every` steps each small molecule gets a
chance, by type and size, to jump 1..`jump` units (up to more than a bilayer
thick) in a random direction, moving as a whole. The jump is kept only if
  - its path passes through lipid tails (it is a crossing; in open water the
    molecule just diffuses as usual), and
  - it lands in clear water (not wedged in the membrane, not on a bead).
"""
import numpy as np
from scipy.spatial import cKDTree

from .elements import INDEX

CARB, TAIL = INDEX["carb"], INDEX["tail"]
TAILS = (TAIL, INDEX["stail"])

# chance per check (every=50) to try a hop; only ~1 in 70 tries is a real crossing
SUGAR_RATE = {1: 0.5, 2: 0.1, 3: 0.02, 4: 0.02}   # carbons in the chain -> rate; 5+ = 0
ION_RATE = 0.005
SMALL_RATE = 0.5   # lone HCN, NH3: small and neutral, cross as fast as formaldehyde


class Leak:
    def __init__(self, every=50, jump=9.0, gap=0.9, sugar_rate=None, ion_rate=ION_RATE, ions=("Na", "K", "Ca", "Mg", "Cl")):
        self.every, self.jump, self.gap = every, jump, gap
        self.sugar_rate = SUGAR_RATE if sugar_rate is None else sugar_rate
        self.ion_rate = ion_rate
        self.ion_types = [INDEX[n] for n in ions if n in INDEX]
        self.hops = 0
        self.wet = None   # optional pos -> bool: hops must land in water (a hop onto dry land got flung back and crashed)

    def rates(self, world):
        """Hop chance for each molecule (0 = never moves this way)."""
        rate = np.zeros(world.n_mol)
        size = np.bincount(world.mol, minlength=world.n_mol)
        carb = world.type == CARB
        carb_mol = np.unique(world.mol[carb])
        for m in carb_mol:
            rate[m] = self.sugar_rate.get(int(size[m]), 0.0)
        ion = np.isin(world.type, self.ion_types) & ~world.frozen
        rate[world.mol[ion]] = self.ion_rate
        small = np.isin(world.type, [INDEX["hcn"], INDEX["nh3"]]) & (size[world.mol] == 1)
        rate[world.mol[small]] = SMALL_RATE
        rate[world.mol[world.frozen]] = 0.0   # rock and dry crust stay put
        return rate

    def leak(self, world):
        rate = self.rates(world)
        movers = np.flatnonzero(world.rng.uniform(size=world.n_mol) < rate)
        if not len(movers):
            return
        L = world.L
        p = world._wrapped()
        tails = np.flatnonzero(np.isin(world.type, TAILS))
        if not len(tails):
            return
        tail_tree = cKDTree(p[tails], boxsize=L)
        along = np.linspace(0.1, 0.9, 9)[:, None]
        moved = False
        for m in movers:
            beads = np.flatnonzero(world.mol == m)
            a = world.rng.uniform(0, 2 * np.pi)
            step = world.rng.uniform(1.0, self.jump) * np.array([np.cos(a), np.sin(a)])
            path = (p[beads[0]] + along * step) % L
            if tail_tree.query(path)[0].min() > 1.0:
                continue   # no membrane in the way: ordinary diffusion handles it
            new = (p[beads] + step) % L
            others = np.ones(world.n, dtype=bool)
            others[beads] = False
            tree = cKDTree(p[others], boxsize=L)
            if tree.query(new)[0].min() < self.gap or tail_tree.query(new)[0].min() < 1.5:
                continue   # must land in clear water, not in a crowd or wedged in the membrane
            if self.wet is not None and not self.wet(new).all():
                continue
            world.pos[beads] += step
            p[beads] = new
            self.hops += 1
            moved = True
        if moved:
            world._pairs = None   # neighbour lists are out of date
