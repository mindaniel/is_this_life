"""2D particle chemistry engine.

Coarse-grained, implicit-solvent model (Cooke-Kremer-Deserno style):
water is not simulated as particles. Its effect is folded into the pair rules:
oily tails attract each other because water pushes them together
(hydrophobic effect). Heads only repel, so they stay at the water surface.

Charges: + and - attract, like charges repel. Salt water weakens the pull
with distance (screening): U = coulomb * qi * qj * exp(-r / debye) / r.
More salt = shorter debye length = charges feel each other less far.

Units: bead size sigma = 1, energy epsilon = 1, mass = 1.
"""
import numpy as np
from scipy.spatial import cKDTree

from . import fast

from . import elements as el

HEAD, TAIL = el.INDEX["head"], el.INDEX["tail"]
TYPE_NAMES = el.NAMES
N_TYPES = len(TYPE_NAMES)


class Rules:
    """Interaction table: the 'chemistry'. Data, not code.

    Built from the element table (elements.py).
    sigma[a, b]  size of the repulsive core between type a and b
    eps[a, b]    depth of oily attraction (0 = no attraction)
    charge[a]    electric charge of type a
    wc           attraction range beyond contact
    coulomb      charge strength. Bjerrum length (where two +1 charges feel
                 1 kT) is ~0.7 nm in water, about one bead, so ~0.6 at kT 0.6
    debye        screening length in sigma. ~2 = pond (~30 mM salt),
                 ~0.5 = seawater
    """

    def __init__(self):
        n = N_TYPES
        size = np.array(el.SIZE)
        oily = np.array(el.OILY)
        self.sigma = (size[:, None] + size[None, :]) / 2
        for (a, b), s in el.PAIR_SIZE.items():
            self.sigma[el.INDEX[a], el.INDEX[b]] = self.sigma[el.INDEX[b], el.INDEX[a]] = s
        self.eps = np.sqrt(oily[:, None] * oily[None, :])
        self.charge = np.array(el.CHARGE, dtype=float)
        assert self.sigma.shape == (n, n)
        self.wc = 1.6
        self.coulomb = 0.6
        self.debye = 2.0
        self.bond_k = 100.0   # spring between neighbours in a molecule
        self.bond_r0 = 1.0
        self.bend_k = 10.0    # spring between beads i and i+2: keeps lipids straight
        self.bend_r0 = 2.0

    @property
    def rc(self):
        return 2 ** (1 / 6) * self.sigma

    @property
    def cutoff(self):
        return float(self.rc.max() + self.wc)

    @property
    def charge_cutoff(self):
        """Beyond 3 screening lengths the charge force is under 5%: ignore it."""
        return 3.0 * self.debye


class World:
    """Periodic 2D box of beads, bonds, and a Langevin thermostat."""

    def __init__(self, box=40.0, rules=None, kT=1.0, gamma=1.0, dt=0.01, seed=0):
        self.L = float(box)
        self.rules = rules or Rules()
        self.kT, self.gamma, self.dt = kT, gamma, dt
        self.rng = np.random.default_rng(seed)
        self.pos = np.zeros((0, 2))
        self.vel = np.zeros((0, 2))
        self.type = np.zeros(0, dtype=int)
        self.mol = np.zeros(0, dtype=int)
        self.bonds = np.zeros((0, 2), dtype=int)
        self.bends = np.zeros((0, 2), dtype=int)
        self.frozen = np.zeros(0, dtype=bool)   # beads that never move (mineral)
        self.n_mol = 0
        self.t = 0
        self.energy = 0.0
        self.flow = None   # optional function pos -> water velocity (currents)
        self.kT_field = None   # optional function pos -> temperature per bead (zones)
        self.external = None   # optional function pos -> (force, energy), e.g. a pond shore
        self._pairs = None
        self._charged_pairs = None
        self._pos_at_build = None
        self.skin = 0.4
        self.fast = True   # compiled force loops (alife/chem/fast.py) when numba is installed
        self.field_every = 10   # steps between re-reading flow and kT_field at bead positions
        # bead types that also push on beads of their own molecule. Lipids are stiff
        # (bends), so they skip that; floppy chains made by reactions (sugars) need it
        self.repel_inside = ()

    # ---- building -------------------------------------------------------
    def add_molecule(self, positions, types, bonds=(), bends=(), frozen=False):
        """Add beads as one molecule. bonds/bends use local bead indices.

        types can be type numbers or element names ("Na", "tail", ...).
        frozen beads stay where they are put (a mineral surface).
        """
        start = len(self.pos)
        positions = np.asarray(positions, dtype=float) % self.L
        types = [TYPE_NAMES.index(t) if isinstance(t, str) else t for t in types]
        self.pos = np.vstack([self.pos, positions])
        v = np.zeros(positions.shape) if frozen else self.rng.normal(0, np.sqrt(self.kT), positions.shape)
        self.vel = np.vstack([self.vel, v])
        self.frozen = np.concatenate([self.frozen, np.full(len(types), frozen)])
        self.type = np.concatenate([self.type, types]).astype(int)
        self.mol = np.concatenate([self.mol, np.full(len(types), self.n_mol)])
        if len(bonds):
            self.bonds = np.vstack([self.bonds, np.asarray(bonds) + start])
        if len(bends):
            self.bends = np.vstack([self.bends, np.asarray(bends) + start])
        self.n_mol += 1
        self._pairs = None

    # ---- reactions: bonds that form and break while running ---------------
    def add_bond(self, i, j):
        """Bond beads i and j (a reaction joins them). Their molecules become one."""
        self.bonds = np.vstack([self.bonds, [i, j]])
        a, b = self.mol[i], self.mol[j]
        if a != b:
            self.mol[self.mol == b] = a
        self._pairs = None

    def break_bond(self, k):
        """Remove bond number k. If that cuts a molecule in two, one piece gets a new id."""
        i, j = self.bonds[k]
        self.bonds = np.delete(self.bonds, k, axis=0)
        mine = self.mol == self.mol[j]
        piece = self._reachable(j, mine)
        if not piece[i]:
            self.mol[piece] = self.n_mol
            self.n_mol += 1
        self._pairs = None

    def _reachable(self, start, within):
        """Beads connected to `start` by bonds, searching only beads in `within`."""
        b = self.bonds[within[self.bonds[:, 0]] & within[self.bonds[:, 1]]]
        seen = np.zeros(self.n, dtype=bool)
        seen[start] = True
        todo = [start]
        while todo:
            x = todo.pop()
            for y in np.r_[b[b[:, 0] == x, 1], b[b[:, 1] == x, 0]]:
                if not seen[y]:
                    seen[y] = True
                    todo.append(y)
        return seen

    @property
    def n(self):
        return len(self.pos)

    # ---- geometry -------------------------------------------------------
    def delta(self, i, j):
        """Shortest vector from bead i to bead j through periodic walls."""
        d = self.pos[j] - self.pos[i]
        return d - self.L * np.round(d / self.L)

    def _wrapped(self):
        p = np.mod(self.pos, self.L)
        p[p >= self.L] = 0.0
        return p

    def _group(self):
        """Beads in the same group skip non-bonded forces: same molecule, unless
        their type is in repel_inside (then each bead is its own group)."""
        if not len(self.repel_inside):
            return self.mol
        group = self.mol.copy()
        alone = np.isin(self.type, self.repel_inside)
        group[alone] = -1 - np.flatnonzero(alone)
        return group

    def _neighbour_pairs(self):
        """Verlet list: rebuild only when something moved more than skin/2."""
        if self._pairs is not None:
            moved = self.pos - self._pos_at_build
            moved -= self.L * np.round(moved / self.L)
            if np.max(np.einsum("ij,ij->i", moved, moved)) < (self.skin / 2) ** 2:
                return self._pairs
        p = self._wrapped()
        if self.fast and fast.HAVE_NUMBA and self.L >= 3 * (self.rules.charge_cutoff + self.skin):
            everyone = np.arange(self.n)
            self._pairs = fast.cell_pairs(p, self.L, self.rules.cutoff + self.skin,
                                          self._group(), self.frozen, everyone)
            charged = np.flatnonzero(self.rules.charge[self.type] != 0)
            no_rock = np.zeros(self.n, dtype=bool)   # charged list never skipped rock pairs
            self._charged_pairs = fast.cell_pairs(p, self.L, self.rules.charge_cutoff + self.skin,
                                                  self._group(), no_rock, charged)
            self._pos_at_build = self.pos.copy()
            return self._pairs
        tree = cKDTree(p, boxsize=self.L)
        pairs = tree.query_pairs(self.rules.cutoff + self.skin, output_type="ndarray")
        # beads in the same molecule interact only through bonds
        group = self._group()
        pairs = pairs[group[pairs[:, 0]] != group[pairs[:, 1]]]
        # frozen beads (rock) never move, so rock-rock forces do nothing
        pairs = pairs[~(self.frozen[pairs[:, 0]] & self.frozen[pairs[:, 1]])]
        self._pairs = pairs

        # charges reach farther, but only charged beads need the long list
        charged = np.flatnonzero(self.rules.charge[self.type] != 0)
        cp = np.zeros((0, 2), dtype=int)
        if len(charged) > 1:
            ctree = cKDTree(p[charged], boxsize=self.L)
            cp = charged[ctree.query_pairs(self.rules.charge_cutoff + self.skin, output_type="ndarray")]
            cp = cp[group[cp[:, 0]] != group[cp[:, 1]]]
        self._charged_pairs = cp
        self._pos_at_build = self.pos.copy()
        return pairs

    # ---- forces ---------------------------------------------------------
    def forces(self):
        if not (self.fast and fast.HAVE_NUMBA):
            return self.forces_numpy()
        r = self.rules
        F = np.zeros_like(self.pos)
        pos = np.ascontiguousarray(self.pos)
        pairs = self._neighbour_pairs()
        energy = fast.pair_forces(pos, self.L, pairs, self.type, r.sigma, r.rc, r.eps, r.wc, F)
        cp = self._charged_pairs
        if cp is not None and len(cp):
            energy += fast.charge_forces(pos, self.L, cp, self.type, r.charge, r.coulomb, r.debye,
                                         r.charge_cutoff, F)
        for idx, k, r0 in ((self.bonds, r.bond_k, r.bond_r0),
                           (self.bends, r.bend_k, r.bend_r0)):
            if len(idx):
                energy += fast.spring_forces(pos, self.L, idx, k, r0, F)
        return F, energy

    def forces_numpy(self):
        """Reference version of forces(), plain numpy (slower, easy to read)."""
        r = self.rules
        F = np.zeros_like(self.pos)
        energy = 0.0

        pairs = self._neighbour_pairs()
        if len(pairs):
            i, j = pairs[:, 0], pairs[:, 1]
            d = self.delta(i, j)
            dist = np.sqrt(np.einsum("ij,ij->i", d, d))
            ti, tj = self.type[i], self.type[j]
            b = r.sigma[ti, tj]
            rc = r.rc[ti, tj]
            eps = r.eps[ti, tj]
            fmag = np.zeros_like(dist)   # >0 pushes apart, <0 pulls together

            rep = dist < rc
            sr6 = (b[rep] / dist[rep]) ** 6
            fmag[rep] = 24 * (2 * sr6 ** 2 - sr6) / dist[rep]
            energy += np.sum(4 * (sr6 ** 2 - sr6) + 1 - eps[rep])

            att = (~rep) & (dist < rc + r.wc) & (eps > 0)
            x = np.pi * (dist[att] - rc[att]) / (2 * r.wc)
            fmag[att] = -eps[att] * np.pi / (2 * r.wc) * np.sin(2 * x)
            energy += np.sum(-eps[att] * np.cos(x) ** 2)

            f = (fmag / np.maximum(dist, 1e-9))[:, None] * d
            self._scatter(F, i, j, f)

        cp = self._charged_pairs
        if cp is not None and len(cp):
            i, j = cp[:, 0], cp[:, 1]
            d = self.delta(i, j)
            dist = np.maximum(np.sqrt(np.einsum("ij,ij->i", d, d)), 0.5)
            near = dist < r.charge_cutoff
            i, j, d, dist = i[near], j[near], d[near], dist[near]
            qq = r.coulomb * r.charge[self.type[i]] * r.charge[self.type[j]]
            screen = np.exp(-dist / r.debye)
            fmag = qq * screen * (1 / dist ** 2 + 1 / (r.debye * dist))   # >0 pushes apart
            cut = r.charge_cutoff
            energy += np.sum(qq * (screen / dist - np.exp(-cut / r.debye) / cut))
            self._scatter(F, i, j, (fmag / dist)[:, None] * d)

        for idx, k, r0 in ((self.bonds, r.bond_k, r.bond_r0),
                           (self.bends, r.bend_k, r.bend_r0)):
            if not len(idx):
                continue
            i, j = idx[:, 0], idx[:, 1]
            d = self.delta(i, j)
            dist = np.sqrt(np.einsum("ij,ij->i", d, d))
            stretch = dist - r0
            f = (-k * stretch / np.maximum(dist, 1e-9))[:, None] * d
            self._scatter(F, i, j, f)
            energy += np.sum(0.5 * k * stretch ** 2)

        return F, energy

    def _scatter(self, F, i, j, f):
        """Add force f to bead j and -f to bead i (Newton's third law)."""
        for ax in range(2):
            F[:, ax] += np.bincount(j, f[:, ax], minlength=self.n)
            F[:, ax] -= np.bincount(i, f[:, ax], minlength=self.n)

    # ---- time -----------------------------------------------------------
    def step(self, n=1, max_force=None):
        """Langevin dynamics: forces + friction + random kicks (heat).

        Friction is drag against the water. If the water flows (self.flow),
        beads are dragged toward the local water velocity instead of rest.
        If self.kT_field is set, each bead feels its local temperature.
        Both fields vary over tens of units while a bead moves ~0.03 per step,
        so they are re-read every `field_every` steps, not every step.
        """
        c1 = np.exp(-self.gamma * self.dt)
        c2 = np.sqrt((1 - c1 ** 2) * self.kT)
        u = 0.0
        for s in range(n):
            if s % self.field_every == 0:
                if self.kT_field is not None:
                    c2 = np.sqrt((1 - c1 ** 2) * np.maximum(self.kT_field(self.pos), 0.0))[:, None]
                u = 0.0 if self.flow is None else self.flow(self.pos)
            F, self.energy = self.forces()
            if self.external is not None:
                f_ext, e_ext = self.external(self.pos)
                F += f_ext
                self.energy += e_ext
            if max_force is not None:
                np.clip(F, -max_force, max_force, out=F)
            self.vel += F * self.dt
            self.vel = u + c1 * (self.vel - u) + c2 * self.rng.standard_normal(self.vel.shape)
            self.vel[self.frozen] = 0.0
            self.pos += self.vel * self.dt
            self.t += 1
