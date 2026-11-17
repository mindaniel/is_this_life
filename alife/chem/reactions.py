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

from . import fast
from .elements import INDEX

CARB, CA, CAT, TAIL, STAIL = INDEX["carb"], INDEX["Ca"], INDEX["cat"], INDEX["tail"], INDEX["stail"]
FA_HEAD, FA_ACID = INDEX["fa_head"], INDEX["fa_acid"]
FUEL, WASTE = INDEX["fuel"], INDEX["waste"]
PIGMENT = INDEX["pigment"]
HEADS = (INDEX["head"], FA_HEAD, FA_ACID)


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
        heads = np.isin(world.type, HEADS)
        n_head = np.bincount(world.mol[heads], minlength=world.n_mol)
        n_tail = np.bincount(world.mol[np.isin(world.type, self.chains)], minlength=world.n_mol)
        size = np.bincount(world.mol, minlength=world.n_mol)
        lone = (n_head == 1) & (n_tail >= 1) & (size == n_head + n_tail)
        return np.flatnonzero(heads & lone[world.mol])

    def head_bonds(self, world):
        """Indices of bonds that join two heads (the two-chain lipids)."""
        heads = np.isin(world.type, HEADS)
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


class Protonation:
    """Fatty-acid heads take up and give off H+ with the water: COOH (fa_acid, neutral) <-> COO- (fa_head).

    Real carboxylic acids: alone in water, half are charged at pH = pKa (~4.8). In a wall,
    charged neighbours push each other, so taking off one more H+ costs more: the wall's
    half-charged point moves up (lab: pH ~7-9 for fatty-acid vesicles). A COOH next to a COO-
    shares the H+ (rules.pair_eps fa_acid-fa_head), which holds such walls together.
    Nothing here sets the wall's charge: it comes out of pH, pKa and the neighbours.

    Every `every` steps, `sweeps` passes of constant-pH Monte Carlo (fast.titrate) over all
    fatty-acid heads in random order: a flip is kept with chance min(1, exp(-dE / kT)), dE = change
    in pair + charge energy - kT ln10 (pH - pKa) for giving off the H+. Proton hops take
    nanoseconds, far quicker than lipids move, so a full pass each check is fair.
    Shortcut (flagged): H+ and counter-ions are implicit (pH is one number everywhere; salt
    screening as before).

    Sealed insides (`closed` = bead numbers of heads facing water inside a vesicle, set by the
    caller, e.g. with `facing`): those heads don't swap H+ with the outside water. A bag's
    inside holds ~0 free H+ (pH 5 in 1000 nm^3 is 0.006 ions), so its H+ live on its heads;
    with `hand_over` they hop between neighbouring heads (COOH -> COO-, count kept).
    Defaults (closed None, hand_over False) give the old behaviour bit for bit.
    """

    def __init__(self, ph, pka=4.8, every=50, sweeps=1, hand_over=False, hop_reach=1.5):
        self.ph, self.pka, self.every, self.sweeps = ph, pka, every, sweeps
        self.hand_over, self.hop_reach = hand_over, hop_reach
        self.closed = None
        self.flips = self.hops = 0

    def lists(self, world):
        """Fatty-acid heads, and their neighbours as CSR lists: short (pair forces) and charge."""
        heads = np.flatnonzero(np.isin(world.type, (FA_HEAD, FA_ACID)))
        where = np.full(world.n, -1)
        where[heads] = np.arange(len(heads))
        pairs = world._neighbour_pairs()
        short = []
        for x, y in ((0, 1), (1, 0)):
            ok = where[pairs[:, x]] >= 0
            short.append(np.stack([where[pairs[ok, x]], pairs[ok, y]], 1))
        sel = np.flatnonzero((world.rules.charge[world.type] != 0) | (where >= 0))
        group = world._group()
        cp = sel[cKDTree(world._wrapped()[sel], boxsize=world.L).query_pairs(world.rules.charge_cutoff,
                                                                              output_type="ndarray")]
        cp = cp[group[cp[:, 0]] != group[cp[:, 1]]]
        charge = []
        for x, y in ((0, 1), (1, 0)):
            ok = where[cp[:, x]] >= 0
            charge.append(np.stack([where[cp[ok, x]], cp[ok, y]], 1))

        def csr(parts):
            e = np.concatenate(parts)
            e = e[np.argsort(e[:, 0], kind="stable")]
            ptr = np.concatenate([[0], np.cumsum(np.bincount(e[:, 0], minlength=len(heads)))])
            return ptr, np.ascontiguousarray(e[:, 1])
        return heads, csr(short), csr(charge)

    def run(self, world, apply=True, order=None):
        """One pass; returns dE per head (before its flip). apply=False changes nothing."""
        heads, (sp, sn), (cpt, cn) = self.lists(world)
        if not len(heads):
            return np.zeros(0)
        r = world.rules
        order = world.rng.permutation(len(heads)) if order is None else order
        u = world.rng.uniform(size=len(heads)) if apply else np.zeros(len(heads))
        mu = world.kT * np.log(10) * (self.ph - self.pka)
        dE = np.zeros(len(heads))
        typ = world.type if apply else world.type.copy()
        open_ = np.ones(len(heads), bool) if self.closed is None else ~np.isin(heads, self.closed)
        pos = np.ascontiguousarray(world.pos)
        n = fast.titrate(pos, world.box(), heads, order, sp, sn, cpt, cn, typ,
                         FA_ACID, FA_HEAD, r.sigma, r.rc, r.eps, r.wc, r.charge, r.coulomb, r.debye,
                         r.charge_cutoff, world.kT, mu, u, apply, dE, open_)
        if apply and self.hand_over:
            a = np.repeat(np.arange(len(heads)), np.diff(sp))   # head pairs from the short list
            where = np.full(world.n, -1)
            where[heads] = np.arange(len(heads))
            b = where[sn]
            ok = b >= 0
            a, b = a[ok], b[ok]
            d = world.pos[heads[a]] - world.pos[heads[b]]
            d -= world.box() * np.round(d / world.box())
            ok = (d ** 2).sum(1) < self.hop_reach ** 2
            a, b = a[ok], b[ok]
            k = world.rng.permutation(len(a))
            h = fast.hand_over(pos, world.box(), heads, a[k], b[k], sp, sn, cpt, cn, typ, FA_ACID, FA_HEAD,
                               r.sigma, r.rc, r.eps, r.wc, r.charge, r.coulomb, r.debye, r.charge_cutoff,
                               world.kT, world.rng.uniform(size=len(a)))
            self.hops += h
            n += h
        if n:
            self.flips += n
            world._pairs = None   # the charged-pair list changed
        return dE

    @staticmethod
    def facing(world, step=0.5):
        """Fatty-acid heads (COOH or COO-) facing water sealed inside a vesicle: the water nearest
        to a point `step` beyond the head (away from its first tail bead) is a pocket (water_pockets3)."""
        from .analysis import water_pockets3
        heads = np.flatnonzero(np.isin(world.type, (FA_HEAD, FA_ACID)))
        if not len(heads):
            return heads
        tail = np.full(world.n, -1)
        bonds = world.bonds
        for x, y in ((0, 1), (1, 0)):
            ok = np.isin(world.type[bonds[:, y]], (TAIL, STAIL))
            tail[bonds[ok, x]] = bonds[ok, y]
        d = world.pos[heads] - world.pos[tail[heads]]
        d -= world.box() * np.round(d / world.box())
        pts = world.pos[heads] + step * d / np.linalg.norm(d, axis=1)[:, None]
        inside = water_pockets3(world, points=pts, near_water=True)[2]
        return heads[inside]

    def react(self, world):
        for _ in range(self.sweeps):
            self.run(world)

    def charged(self, world, beads=None):
        """Fraction of fatty-acid heads (among `beads`, default all) that are COO-."""
        t = world.type if beads is None else world.type[beads]
        n = np.isin(t, (FA_HEAD, FA_ACID)).sum()
        return float((t == FA_HEAD).sum() / n) if n else 0.0


class Burn:
    """A catalyst burns fuel: fuel -> waste + H+ (energy goes downhill, never back).

    Real example: an activated molecule (an acid anhydride or a thioester) that a catalyst splits
    with water, freeing a carboxylic acid's H+. Every `every` steps each fuel bead within `cat_reach`
    of a "cat" bead burns with chance p_burn. Its H+ goes onto the nearest COO- head within
    `h_reach` (COO- -> COOH); with none in reach it doesn't burn (nothing to take the H+).
    heat (in kT): energy each burn gives out as heat, a kick to the `heat_beads` nearest beads
    (random directions, so their motion gains `heat` kT on average). The thermostat (water)
    carries it off within ~100 steps, as real heat leaves a nanometre spot in picoseconds.
    Shortcuts (flagged): the H+ jumps straight to that head (real: diffuses there in nanoseconds);
    burn chance doesn't depend on temperature.
    """

    def __init__(self, every=50, cat_reach=2.5, p_burn=0.1, h_reach=4.0, heat=0.0, heat_beads=8):
        self.every, self.cat_reach, self.p_burn, self.h_reach = every, cat_reach, p_burn, h_reach
        self.heat, self.heat_beads = heat, heat_beads
        self.burns = 0

    def react(self, world):
        fuel = np.flatnonzero(world.type == FUEL)
        cat = np.flatnonzero(world.type == CAT)
        if not len(fuel) or not len(cat):
            return 0
        p = world._wrapped()
        near = cKDTree(p[cat], boxsize=world.L).query(p[fuel], distance_upper_bound=self.cat_reach)[0] < np.inf
        fuel = fuel[near]
        if not len(fuel):
            return 0
        fuel = fuel[world.rng.uniform(size=len(fuel)) < self.p_burn]
        base = np.flatnonzero(world.type == FA_HEAD)
        if not len(fuel) or not len(base):
            return 0
        tree = cKDTree(p[base], boxsize=world.L)
        n = 0
        for f in fuel:
            dist, k = tree.query(p[f], k=min(8, len(base)), distance_upper_bound=self.h_reach)
            for d, j in zip(np.atleast_1d(dist), np.atleast_1d(k)):
                if d < np.inf and world.type[base[j]] == FA_HEAD:
                    world.type[base[j]] = FA_ACID
                    world.type[f] = WASTE
                    n += 1
                    if self.heat:
                        self.kick(world, p, f)
                    break
        if n:
            self.burns += n
            world._pairs = None
        return n

    def kick(self, world, p, f):
        """Give `heat` kT of motion to the beads around bead f (heat from the burn)."""
        heat_kick(world, p, f, self.heat, self.heat_beads)


def heat_kick(world, p, f, heat, beads=8, tree=None):
    """Give `heat` kT of motion to the `beads` beads nearest bead f (random directions)."""
    near = (tree or cKDTree(p, boxsize=world.L)).query(p[f], k=beads)[1]
    near = near[~world.frozen[near]]
    v = world.rng.standard_normal((len(near), world.dim))
    v *= np.sqrt(2 * heat * world.kT / len(near)) / np.linalg.norm(v, axis=1)[:, None]
    world.vel[near] += v


class PHShocks:
    """pH of the outside water with sudden events: at each (step, jump) the pH jumps by `jump`
    (volcanic gas or acid rain: down; alkaline spring water: up), then drifts back to `base`
    with time constant `relax` steps (mixing with the rest of the pond). Shocks add up.
    Shortcut (flagged): the whole outside changes at once (one pH number).
    A third number holds the jump for that many steps before it starts to fade (a long dry
    spell of volcanic gas, a spring that keeps flowing).
    """

    def __init__(self, base, shocks=(), relax=5000):
        self.base, self.relax = base, relax
        self.shocks = [tuple(s) + (0,) * (3 - len(s)) for s in shocks]   # (step, jump, hold)

    @staticmethod
    def parse(text):
        """"10000:-3,60000:+2:20000" -> [(10000, -3.0, 0), (60000, 2.0, 20000)]"""
        if not text:
            return []
        out = []
        for item in text.split(","):
            v = item.split(":")
            out.append((int(v[0]), float(v[1]), int(v[2]) if len(v) > 2 else 0))
        return out

    @staticmethod
    def random(rate, steps, seed, size=(1.0, 3.0), hold=0):
        """Random events: on average `rate` per 10k steps (Poisson), each up or down by a size
        drawn from `size`. Own random stream, so the world's random numbers are untouched."""
        rng = np.random.default_rng(seed + 7919)
        out, t = [], 0.0
        while True:
            t += rng.exponential(10000 / rate)
            if t >= steps:
                return out
            out.append((int(t), float(rng.choice((-1, 1)) * rng.uniform(*size)), hold))

    def __call__(self, t):
        ph = self.base
        for s, j, hold in self.shocks:
            if t >= s:
                ph += j if t < s + hold else j * np.exp(-(t - s - hold) / self.relax)
        return ph


class LightPump:
    """Light-driven H+ pump (bacteriorhodopsin-like): a pigment in the wall absorbs light and uses
    the energy to push an H+ from the outside water onto a head facing the sealed inside.

    Every `every` steps while the light is on, each pigment bead absorbs a photon with chance
    p_abs and tries to put an H+ on the nearest COO- facing the inside (`closed`) within `reach`.
    Cost = energy change of that head (charges, pairing) + kT ln10 (pH - pKa) (taking an H+ from
    the outside water) - `gain` kT stored from the photon. Kept with chance min(1, exp(-cost/kT)):
    the more acid the inside already is, the harder the push, so the pump stalls at a gradient.
    Light: on for the first half of each `period` (day), off the second half (night).
    Shortcuts (flagged): the pigment is a fatty acid whose head is a pigment bead (all oily, so it
    sits in the oil); a real pump is a protein spanning the wall, built to point one way. Here the
    direction comes from "inside heads only". gain 10 kT ~ what real pumps keep of a ~80 kT photon.

    Other pigments: `pigment` = function world -> bead indices that absorb (e.g. dark HCN polymer,
    nitrogen.dark_beads); then `any_head=True` lets the H+ go onto any nearby COO- (no built-in
    direction), and gain 0 turns the H+ part off. heat (kT): each absorbed photon also warms the
    `heat_beads` beads around the pigment (heat_kick); absorbed = photons so far.
    sun (sun.Sun): smooth day instead of on/off, absorb chance p_abs x brightness; then gain and heat
    are in kT at Earth's temperature (sun.kT_earth), since a photon's energy doesn't change with
    the pond's temperature.
    """

    def __init__(self, every=50, reach=3.5, p_abs=0.2, gain=10.0, period=20000, pka=4.8,
                 pigment=None, any_head=False, heat=0.0, heat_beads=8, sun=None):
        self.every, self.reach, self.p_abs, self.gain, self.period, self.pka = every, reach, p_abs, gain, period, pka
        self.pigment, self.any_head, self.heat, self.heat_beads = pigment, any_head, heat, heat_beads
        self.sun = sun
        self.pumps = self.tries = self.absorbed = 0

    def light(self, t):
        if self.sun:
            return self.sun.brightness(t) > 0
        return self.period <= 0 or (t % self.period) < self.period / 2

    def react(self, world, ph, closed, prot):
        if not self.light(world.t):
            return 0
        pig = np.flatnonzero(world.type == PIGMENT) if self.pigment is None else self.pigment(world)
        p_abs = self.p_abs * self.sun.brightness(world.t) if self.sun else self.p_abs
        pig = pig[world.rng.uniform(size=len(pig)) < p_abs]
        unit = self.sun.kT_earth / world.kT if self.sun else 1.0   # gain, heat: in kT at Earth's temperature
        if not len(pig):
            return 0
        self.absorbed += len(pig)
        if self.heat:
            p = world._wrapped()
            tree = cKDTree(p, boxsize=world.L)
            for g in pig:
                heat_kick(world, p, g, self.heat * unit, self.heat_beads, tree)
        if self.any_head:
            closed = np.flatnonzero(world.type == FA_HEAD)
        if not self.gain or closed is None or not len(closed):
            return 0
        heads, (sp, sn), (cpt, cn) = prot.lists(world)
        where = np.full(world.n, -1)
        where[heads] = np.arange(len(heads))
        p = world._wrapped()
        tree = cKDTree(p[closed], boxsize=world.L)
        r = world.rules
        pos = np.ascontiguousarray(world.pos)
        mu = world.kT * np.log(10) * (ph - self.pka)
        n = 0
        for g in pig:
            dist, k = tree.query(p[g], k=min(12, len(closed)), distance_upper_bound=self.reach)
            for d, j in zip(np.atleast_1d(dist), np.atleast_1d(k)):
                if d == np.inf:
                    break
                i = closed[j]
                if world.type[i] != FA_HEAD:
                    continue
                self.tries += 1
                cost = fast._flip_energy(pos, world.box(), where[i], i, FA_ACID, sp, sn, cpt, cn, world.type,
                                         r.sigma, r.rc, r.eps, r.wc, r.charge, r.coulomb, r.debye,
                                         r.charge_cutoff, mu - self.gain * unit * world.kT)
                if cost <= 0 or world.rng.uniform() < np.exp(-cost / world.kT):
                    world.type[i] = FA_ACID
                    n += 1
                break
        if n:
            self.pumps += n
            world._pairs = None
        return n
