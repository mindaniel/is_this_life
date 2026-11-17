"""e16: a catalyst inside: does the vesicle that makes two-chain lipid win?

Two fatty-acid vesicles side by side (A at x = Lx/4, B at 3 Lx/4), same size,
same sugar inside and out. A also holds --cats catalyst beads ("cat": a short
peptide / ribozyme, too big to leave); B holds the same number of extra sugar
beads instead (same osmotic push, no chemistry). The catalyst joins two fatty
acids of A's inner wall into one two-chain lipid (LipidJoin in reactions.py).
Lab result to compare (Budin & Szostak 2011; Adamala & Szostak 2013): two-chain
lipids stay put, so the vesicle that makes them holds its fatty acids and takes
up loose ones; its neighbour shrinks. Nothing is coded to make A grow.

Logged every --every steps, for A and B (the lipid cluster lining each one's
inner beads: the sugar and catalyst/extra sugar put inside at the start): fatty-acid chains in the wall,
two-chain lipids in it, share of its inner beads still sealed in a pocket (A_in, B_in),
sealed pockets; loose lipids; joins so far. A_fa/B_fa: single-chain fatty acids in the wall;
A_len/B_len: their mean tail length (beads; with --fa-lengths).
A_left/B_left: how many left it since the last row (Hops: checked every 50 steps, so a fatty
acid that leaves and comes back between two checks is missed).
With --ph: fatty-acid heads switch between COOH and COO- (reactions.Protonation; --pka, --hbond
= COOH-COO- pairing); A_q/B_q = charged share of the fatty-acid heads in each wall.
--alcohol f: that share of fatty-acid molecules are fatty alcohols (fa_oh, OH head, H-bonds to COO-).
--sealed-h: heads facing a sealed inside keep their H+ (no swap with the outside water, only
hops between neighbouring heads); A_qin/B_qin = charged share of those inner heads.
--fuel N: N fuel beads inside each bag; A's catalyst burns them (reactions.Burn: fuel -> waste,
its H+ onto a nearby COO-), so A's inside turns acid. Turns on --sealed-h. burns = so far.
--shocks "t:dpH,...": sudden pH jumps of the outside water, drifting back (--relax steps);
ph column = outside pH. --acid-oily x: COOH heads a little oily (real COOH is far less
water-loving than COO-; lets neutral fatty acids slip into the oil, as at low pH in the lab).
A third number per shock holds it that many steps before it fades ("5000:-3:20000");
--shock-rate r adds random shocks (r per 10k steps, size --shock-size "1:3", up or down).
--heat h: each burn gives out h kT of heat to the beads around it (Burn.kick).
--pump N: N of A's fatty acids carry a pigment head (light-driven H+ pump, reactions.LightPump:
pushes H+ from outside onto A's inner heads while the light is on); --light-period P (day = first
half), 0 = always on; --pump-gain (kT kept per photon). pumps = so far, light = 1/0.
--rain "t:frac,...": rain adds water (frac of the volume) over --rain-over steps, evaporation takes
it back over --evap steps (0 = never). Each bag moves as a whole with its contents; the outside
gets more dilute. water = volume now / at start.
--save-every K: save the whole world (World.save) to <name>_state.pkl every K steps and at the end.
--start state.pkl: carry on from a saved world (no build) for --steps more steps, with whatever
conditions the flags give now (rules and bead set come from the file; --reseed s: new random
stream). Saved worlds can be joined first: experiments/merge_worlds.py.
--hcn N / --nh3 M: add N HCN and M ammonia beads in the water outside the bags (at the start, or
when carrying on). With any HCN in the world, the nitrogen chemistry runs (nitrogen.Nitrogen,
rates calibrated at this kT): HCN chains grow where HCN is concentrated (experiments/hcn_pond.py).
--glow: pure-HCN chains of --dark-len or more (dark polymer, real HCN polymer is brown-black)
absorb light (--p-abs, --light-period): each photon gives --glow-heat kT of heat around it, and
with --glow-gain g > 0 tries to put an H+ onto a nearby COO- (photo-acid; speculative, flagged).
Logs hcn (lone), dark (beads in dark chains), longest, dark_wall (dark beads touching a lipid),
absorbed (photons so far), glow_acid (H+ put on by light so far).
--tar-len n: HCN chains of n or more turn into "tar" beads (long HCN polymer is insoluble), oiliness
--tar-oily (a guess). Logs dark_tail (dark beads touching a lipid tail: sunk into a wall).
--sun-d d: a smooth sun (alife/chem/sun.py) for a planet at distance d (Earth = 1): light rises and
falls through the day (brightness 1/d^2), the temperature follows (mean --kT-earth x d^-0.5, small
afternoon-peaked swing), set every 50 steps. Reaction rates keep their Earth calibration
(kT0 = --kT-earth), so hotter = faster. Logs kT and sun (brightness).

    python experiments/e16_catalyst.py --steps 200000 --video runs/video/e16_s0.mp4
    python experiments/e16_catalyst.py --cats 0      # control: no catalyst anywhere
"""
import argparse
import csv
import os
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from alife.chem.analysis import cluster_labels, water_pockets3   # noqa: E402
from alife.chem.elements import INDEX   # noqa: E402
from alife.chem.nitrogen import HCN, NH3, TAR, Nitrogen, census, dark_beads   # noqa: E402
from alife.chem.reactions import HEADS, PIGMENT, Burn, LightPump, LipidJoin, PHShocks, Protonation   # noqa: E402
from alife.chem.sun import Sun   # noqa: E402
from alife.chem.world import STAIL, TAIL, TAILS, World   # noqa: E402
from e15_split3d import feed, mixed_vesicle, offset, parse_lengths, scatter_fast, shear_flow   # noqa: E402

FIELDS = ["step", "A_chains", "A_double", "B_chains", "B_double", "A_in", "B_in", "same_cluster", "pockets", "volumes",
          "loose", "joins", "added", "A_fa", "B_fa", "A_len", "B_len", "A_left", "B_left", "A_q", "B_q",
          "A_oh", "B_oh", "A_qin", "B_qin", "ph", "burns", "pumps", "light", "water",
          "hcn", "dark", "longest", "dark_wall", "dark_tail", "absorbed", "glow_acid", "kT", "sun"]
FA_OH = INDEX["fa_oh"]
CHAIN_HEADS = HEADS + (FA_OH,)   # every single-tail head (fatty alcohols too) counts as a chain


def build(args):
    box = (args.box_x, args.box, args.box)
    w = World(box=box, kT=args.kT, seed=args.seed, dim=3)
    w.rules.set_oily("stail", args.short_oily)
    if args.ph is not None:
        w.rules.set_pair_eps("fa_acid", "fa_head", args.hbond)
        if args.acid_oily:
            w.rules.set_oily("fa_acid", args.acid_oily)
    L = np.array(box)
    centres = [np.array([L[0] / 4, L[1] / 2, L[2] / 2]), np.array([3 * L[0] / 4, L[1] / 2, L[2] / 2])]
    n = sum(mixed_vesicle(w, c, args.radius, args.mix, args.charged, lengths=args.lengths, acid=args.acid,
                          alcohol=args.alcohol) for c in centres)
    r_in, outer = args.radius - 3.5, args.radius + 2.5
    inside = [lambda p, c=c: np.linalg.norm(offset(w, p, c), axis=1) < r_in for c in centres]
    n_in = int(round(args.osm * 4 / 3 * np.pi * r_in ** 3))
    w.markers = []   # bead numbers inside A and B (beads are only ever appended, so these stay valid)
    for keep, extra in ((inside[0], "cat"), (inside[1], "sugar")):
        start = w.n
        scatter_fast(w, n_in, keep)
        scatter_fast(w, args.cats, keep, name=extra)   # A: catalyst; B: as many plain sugar beads
        if args.fuel:
            scatter_fast(w, args.fuel, keep, name="fuel")   # both bags: same fuel, only A can burn it
        w.markers.append(np.arange(start, w.n))
    if args.pump:   # pigment heads on N random fatty acids of A's wall
        fah = np.flatnonzero(np.isin(w.type, (INDEX["fa_head"], INDEX["fa_acid"], INDEX["head"]))
                             & (np.linalg.norm(offset(w, w.pos, centres[0]), axis=1) < args.radius + 3)
                             & np.isin(w.mol, w.mol[w.type == STAIL]))
        w.type[w.rng.choice(fah, args.pump, replace=False)] = PIGMENT
    n_out = int(round(args.osm * (np.prod(L) - 2 * 4 / 3 * np.pi * (outer + 1.0) ** 3)))
    scatter_fast(w, n_out, lambda p: np.all([np.linalg.norm(offset(w, p, c), axis=1) > outer + 1.0
                                             for c in centres], axis=0))
    w.step(300, max_force=50)
    w.t = 0
    w.n_start = w.n_mol
    if args.shear:
        w.flow = shear_flow(args.box, args.shear)
    return w, n


def bag_centre(w, markers):
    """Middle of a bag's contents (marker beads), taking the periodic box into account."""
    d = w.pos[markers] - w.pos[markers[0]]
    d -= w.L * np.round(d / w.L)
    return (w.pos[markers[0]] + d.mean(0)) % w.L


def add_nitrogen(w, args):
    """--hcn / --nh3 beads into the water outside both bags (radius + 4 from each bag's middle)."""
    centres = [bag_centre(w, m) for m in w.markers]
    away = args.radius + 4.0
    outside = lambda p: np.all([np.linalg.norm(offset(w, p, c), axis=1) > away for c in centres], axis=0)
    for name, count in (("hcn", args.hcn), ("nh3", args.nh3)):
        if count:
            scatter_fast(w, count, outside, name=name)


def nitrogen_row(w, glow, dark_len):
    c = census(w)
    h = c["hcn_chains"]
    dark = dark_beads(w, dark_len)
    lipid = np.flatnonzero(np.isin(w.type, TAILS + HEADS))
    tails = np.flatnonzero(np.isin(w.type, TAILS))
    touch = inside = 0
    if len(dark) and len(lipid):
        p = w._wrapped()
        touch = int((cKDTree(p[lipid], boxsize=w.L).query(p[dark], distance_upper_bound=1.5)[0] < np.inf).sum())
        inside = int((cKDTree(p[tails], boxsize=w.L).query(p[dark], distance_upper_bound=1.5)[0] < np.inf).sum())
    row = {"hcn": c["hcn"], "dark": len(dark), "longest": int(np.flatnonzero(h)[-1]) if h[2:].any() else 1,
           "dark_wall": touch, "dark_tail": inside}
    if glow:
        row["absorbed"], row["glow_acid"], row["light"] = glow.absorbed, glow.pumps, int(glow.light(w.t))
    return row


def wall_of(w, markers, lab, lipid_beads):
    """Lipid cluster lining the marker beads: majority vote of each marker's nearest lipid bead."""
    tree = cKDTree(w._wrapped()[lipid_beads], boxsize=w.L)
    _, k = tree.query(w._wrapped()[markers])
    return int(np.bincount(lab[w.mol[lipid_beads][k]]).argmax())


def single_fa(w, rx):
    """Per molecule: a fatty acid (short tail) not joined into a two-chain lipid."""
    fa = np.zeros(w.n_mol, bool)
    fa[w.mol[w.type == STAIL]] = True
    fa[w.mol[w.bonds[rx.head_bonds(w)][:, 0]]] = False
    fa[w.mol[w.type == PIGMENT]] = False   # pigment carriers are not counted as fatty acids
    return fa


class Hops:
    """Counts fatty acids leaving A's and B's wall: in the wall at one check, still a single fatty
    acid but out of that wall at the next. Analysis only (draws no random numbers)."""

    def __init__(self):
        self.prev, self.left = None, [0, 0]

    def check(self, w, rx):
        lab = cluster_labels(w)
        lipid_beads = np.flatnonzero(np.isin(w.mol, np.unique(w.mol[np.isin(w.type, TAILS)])))
        fa = single_fa(w, rx)
        now = [fa & (lab == wall_of(w, m, lab, lipid_beads)) for m in w.markers]
        if self.prev is not None:
            for k in range(2):
                was = np.zeros(w.n_mol, bool)
                was[:len(self.prev[k])] = self.prev[k]
                self.left[k] += int((was & fa & ~now[k]).sum())
        self.prev = now

    def take(self):
        out, self.left = self.left, [0, 0]
        return out


class Rain:
    """Rain adds water, evaporation takes it back (implicit water: the box grows / shrinks).
    Each event (step, frac): the volume grows by frac x start volume over `over` steps, then
    shrinks back over `evap` steps (0 = stays). Lipid clumps move as a whole, and beads sealed
    in a bag move with that bag (World.rescale group), so only the outside gets more dilute."""

    def __init__(self, events, over=2000, evap=20000):
        self.events, self.over, self.evap = events, over, evap
        self.now = 1.0

    @staticmethod
    def parse(text):
        return [(int(a), float(b)) for a, b in (x.split(":") for x in text.split(","))] if text else []

    def target(self, t):
        v = 1.0
        for s, frac in self.events:
            if t >= s:
                up = min(1.0, (t - s) / self.over)
                down = min(1.0, max(0.0, t - s - self.over) / self.evap) if self.evap else 0.0
                v += frac * (up - down)
        return v

    def apply(self, w):
        goal = self.target(w.t)
        if abs(goal / self.now - 1) < 0.002:
            return
        group = cluster_labels(w)
        lipid_beads = np.flatnonzero(np.isin(w.mol, np.unique(w.mol[np.isin(w.type, TAILS)])))
        allm = np.concatenate(w.markers)
        sealed = water_pockets3(w, beads=allm)[2]
        lab = group.copy()
        for m in w.markers:   # a bag's sealed contents ride with its wall
            c = wall_of(w, m, lab, lipid_beads)
            inside = m[sealed[np.isin(allm, m)]]
            group[w.mol[inside]] = c
        w.rescale((goal / self.now) ** (1 / 3), group, unwrap=True)
        self.now = goal


def measure(w, rx, prot=None):
    lab = cluster_labels(w)
    tails = np.isin(w.type, TAILS)
    lipids = np.unique(w.mol[tails])
    lipid_beads = np.flatnonzero(np.isin(w.mol, lipids))
    size = np.bincount(lab[lipids], minlength=lab.max() + 1)
    heads = np.isin(w.type, CHAIN_HEADS)
    double = np.zeros(w.n_mol, bool)
    double[w.mol[w.bonds[rx.head_bonds(w)][:, 0]]] = True
    row = {"step": w.t}
    walls = {"A": wall_of(w, w.markers[0], lab, lipid_beads), "B": wall_of(w, w.markers[1], lab, lipid_beads)}
    for k, c in walls.items():
        mine = heads & (lab[w.mol] == c)
        row[f"{k}_chains"] = int(mine.sum())                         # one head per fatty-acid chain
        row[f"{k}_double"] = int(np.unique(w.mol[mine & double[w.mol]]).size)
    row["same_cluster"] = int(walls["A"] == walls["B"])
    vol, _, sealed = water_pockets3(w, beads=np.concatenate(w.markers))
    row["A_in"] = round(float(sealed[:len(w.markers[0])].mean()), 3)
    row["B_in"] = round(float(sealed[len(w.markers[0]):].mean()), 3)
    row["pockets"] = len(vol)
    row["volumes"] = " ".join(str(int(v)) for v in sorted(vol, reverse=True))
    row["loose"] = int((size[lab[lipids]] == 1).sum())
    row["joins"] = rx.joins - rx.breaks
    row["added"] = w.n_mol - w.n_start
    fa = single_fa(w, rx)
    n_tail = np.bincount(w.mol[w.type == STAIL], minlength=w.n_mol)
    for k, c in walls.items():
        mine = fa & (lab == c)
        row[f"{k}_fa"] = int(mine.sum())
        row[f"{k}_len"] = round(float(n_tail[mine].mean()), 3) if mine.any() else 0
        fah = np.isin(w.type, (INDEX["fa_head"], INDEX["fa_acid"])) & (lab[w.mol] == c)
        row[f"{k}_q"] = round(float((w.type[fah] == INDEX["fa_head"]).mean()), 3) if fah.any() else ""
        row[f"{k}_oh"] = int(((w.type == FA_OH) & (lab[w.mol] == c)).sum())
        if prot is not None and prot.closed is not None:
            inn = np.zeros(w.n, bool)
            inn[prot.closed] = True
            inn &= fah
            row[f"{k}_qin"] = round(float((w.type[inn] == INDEX["fa_head"]).mean()), 3) if inn.any() else ""
    return row


class Film:
    """Left: a thin slice through both vesicles (z middle), whole box in x-y. Right: wall size over time."""

    def __init__(self, w, path, fps=6):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from alife.viz import Movie

        class V:
            pass
        self.fig, (self.a, self.b) = plt.subplots(1, 2, figsize=(13, 5), gridspec_kw={"width_ratios": [2, 1]})
        v = V()
        v.fig = self.fig
        self.movie = Movie(v, path, fps=fps)
        self.hist = []

    def grab(self, w, row, rx):
        p = w._wrapped()
        L = np.broadcast_to(w.L, (3,))
        cut = np.abs(p[:, 2] - L[2] / 2) < 1.5
        tails, heads = np.isin(w.type, TAILS), np.isin(w.type, CHAIN_HEADS)
        double = np.zeros(w.n_mol, bool)
        double[w.mol[w.bonds[rx.head_bonds(w)][:, 0]]] = True
        two = double[w.mol]
        self.a.clear()
        for m, col, sz in ((w.type == INDEX["sugar"], "#DDDDDD", 4), (tails & ~two, "#EF9F27", 6),
                           (heads & ~two, "#D85A30", 8), (two & (tails | heads), "#7B2FBE", 10),
                           (w.type == INDEX["cat"], "#2E8B57", 22)):
            self.a.scatter(p[cut & m, 0], p[cut & m, 1], s=sz, c=col)
        cat = w.type == INDEX["cat"]   # catalysts out of the slice too (faint), so you see where A is
        self.a.scatter(p[cat & ~cut, 0], p[cat & ~cut, 1], s=10, c="#2E8B57", alpha=0.3)
        self.a.set_xlim(0, L[0]), self.a.set_ylim(0, L[1]), self.a.set_aspect("equal")
        self.a.set_title("slice through both vesicles: A (left) has catalyst (green); "
                         "purple = two-chain lipid", fontsize=9)
        self.hist.append((row["step"] / 1000, row["A_chains"], row["B_chains"], row["A_double"]))
        h = np.array(self.hist)
        self.b.clear()
        self.b.plot(h[:, 0], h[:, 1], c="#2E8B57", label="A wall (catalyst)")
        self.b.plot(h[:, 0], h[:, 2], c="#888888", label="B wall (no catalyst)")
        self.b.plot(h[:, 0], 2 * h[:, 3], c="#7B2FBE", ls="--", label="chains in two-chain lipids (A)")
        self.b.set_xlabel("step (thousands)"), self.b.set_ylabel("fatty-acid chains")
        self.b.legend(fontsize=8)
        self.fig.suptitle(f"step {w.t}", y=0.02, va="bottom", fontsize=10)
        self.movie.grab()

    def close(self):
        self.movie.close()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--radius", type=float, default=10.0)
    p.add_argument("--box", type=float, default=48.0, help="y and z size")
    p.add_argument("--box-x", type=float, default=80.0, help="x size (the two vesicles sit along x)")
    p.add_argument("--mix", type=float, default=0.5, help="fraction of fatty acids (short tail) in the walls")
    p.add_argument("--short-oily", type=float, default=0.85)
    p.add_argument("--charged", type=float, default=0.5, help="fraction of fatty-acid heads charged COO-")
    p.add_argument("--osm", type=float, default=0.15)
    p.add_argument("--cats", type=int, default=20, help="catalyst beads in A (B gets as many sugar beads)")
    p.add_argument("--p-join", type=float, default=0.2, help="join chance per check with a catalyst near")
    p.add_argument("--no-cat", type=float, default=0.0, help="join chance factor without a catalyst")
    p.add_argument("--p-break", type=float, default=0.0, help="chance per check that a two-chain lipid splits")
    p.add_argument("--fa-only", action="store_true", help="only fatty acids (short tails) join, not long-tail lipids")
    p.add_argument("--length", type=float, default=None,
                   help="rest length of the head-head bond (default bond 1.0 squeezes the pair to ~0.89)")
    p.add_argument("--fa-lengths", default=None,
                   help='random fatty-acid tail lengths, beads:weight, e.g. "1:0.3,2:0.4,3:0.3" (default: all 2)')
    p.add_argument("--ph", type=float, default=None,
                   help="heads switch COOH <-> COO- at this pH (default: fixed, --charged of them COO-)")
    p.add_argument("--pka", type=float, default=4.8, help="pKa of a lone fatty acid in water")
    p.add_argument("--hbond", type=float, default=0.5, help="COOH-COO- pairing strength (tail-tail = 1)")
    p.add_argument("--alcohol", type=float, default=0.0, help="share of fatty-acid molecules that are fatty alcohols")
    p.add_argument("--sealed-h", action="store_true", help="heads facing a sealed inside keep their H+ (needs --ph)")
    p.add_argument("--fuel", type=int, default=0, help="fuel beads inside each bag (A's catalyst burns them)")
    p.add_argument("--p-burn", type=float, default=0.1, help="burn chance per check for fuel next to a catalyst")
    p.add_argument("--shocks", default="", help='outside pH jumps "step:dpH,..." e.g. "5000:-3"')
    p.add_argument("--relax", type=float, default=5000, help="steps for a pH shock to fade (time constant)")
    p.add_argument("--acid-oily", type=float, default=0.0, help="oiliness of COOH heads (fa_acid)")
    p.add_argument("--shock-rate", type=float, default=0.0, help="random pH shocks per 10k steps")
    p.add_argument("--shock-size", default="1:3", help="random shock size range (pH units)")
    p.add_argument("--heat", type=float, default=0.0, help="heat per burn (kT) given to nearby beads")
    p.add_argument("--pump", type=int, default=0, help="pigment pumps in A's wall (light-driven H+ in)")
    p.add_argument("--pump-gain", type=float, default=10.0, help="kT kept per absorbed photon")
    p.add_argument("--p-abs", type=float, default=0.2, help="photon absorbed per pigment per check (light on)")
    p.add_argument("--light-period", type=int, default=20000, help="day+night length in steps (0 = always on)")
    p.add_argument("--rain", default="", help='rain events "step:frac,..." (frac of volume added)')
    p.add_argument("--rain-over", type=int, default=2000)
    p.add_argument("--evap", type=int, default=20000, help="steps for the rain water to evaporate (0 = never)")
    p.add_argument("--save-every", type=int, default=0, help="save the world every N steps (and at the end)")
    p.add_argument("--start", default=None, help="carry on from a saved world (.pkl) for --steps more steps")
    p.add_argument("--reseed", type=int, default=None, help="with --start: new random stream from this seed")
    p.add_argument("--hcn", type=int, default=0, help="HCN beads added outside the bags")
    p.add_argument("--nh3", type=int, default=0, help="ammonia beads added outside the bags")
    p.add_argument("--glow", action="store_true", help="dark HCN polymer absorbs light")
    p.add_argument("--dark-len", type=int, default=5, help="HCN chain length that counts as dark pigment")
    p.add_argument("--glow-heat", type=float, default=10.0, help="heat (kT) per absorbed photon")
    p.add_argument("--glow-gain", type=float, default=0.0, help="photo-acid: kT per photon to put H+ on a COO- (0 = off)")
    p.add_argument("--tar-len", type=int, default=0, help="HCN chains of this length or more turn into oily tar (0 = off)")
    p.add_argument("--tar-oily", type=float, default=0.5, help="oiliness of tar beads (a guess)")
    p.add_argument("--sun-d", type=float, default=0.0, help="smooth sun, planet at this distance (Earth = 1; 0 = off)")
    p.add_argument("--kT-earth", type=float, default=1.1, help="temperature at distance 1")
    p.add_argument("--feed", type=int, default=0, help="add --batch loose fatty acids every N steps (0 = never)")
    p.add_argument("--batch", type=int, default=1)
    p.add_argument("--shear", type=float, default=0.0)
    p.add_argument("--kT", type=float, default=1.1)
    p.add_argument("--steps", type=int, default=200000)
    p.add_argument("--every", type=int, default=1000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", default="runs/e16")
    p.add_argument("--video", default=None)
    p.add_argument("--fps", type=int, default=6)
    args = p.parse_args()
    args.lengths = parse_lengths(args.fa_lengths)
    args.acid = "head" if args.ph is None else "fa_acid"
    if args.fuel or args.pump:
        args.sealed_h = True
    if (args.sealed_h or args.shocks or args.shock_rate or args.acid_oily) and args.ph is None:
        p.error("--sealed-h, --fuel, --pump, --shocks and --acid-oily need --ph")

    name = (f"cats{args.cats}_p{args.p_join:g}_nc{args.no_cat:g}_br{args.p_break:g}_mix{args.mix:g}"
            f"_q{args.charged:g}_f{args.feed}x{args.batch}_sh{args.shear:g}" + ("_fa" if args.fa_only else "")
            + (f"_len{args.length:g}" if args.length else "")
            + ("_mix" + "-".join(map(str, args.lengths[0])) if args.lengths else "")
            + (f"_ph{args.ph:g}_pka{args.pka:g}_hb{args.hbond:g}" if args.ph is not None else "")
            + (f"_oh{args.alcohol:g}" if args.alcohol else "") + ("_sealH" if args.sealed_h else "")
            + (f"_fuel{args.fuel}_pb{args.p_burn:g}" if args.fuel else "")
            + ("_shock" + args.shocks.replace(":", "x").replace(",", "_") + f"_r{args.relax:g}" if args.shocks else "")
            + (f"_ao{args.acid_oily:g}" if args.acid_oily else "")
            + (f"_srate{args.shock_rate:g}x{args.shock_size.replace(':', '-')}_r{args.relax:g}" if args.shock_rate else "")
            + (f"_heat{args.heat:g}" if args.heat else "")
            + (f"_pump{args.pump}_g{args.pump_gain:g}_pa{args.p_abs:g}_P{args.light_period}" if args.pump else "")
            + ("_rain" + args.rain.replace(":", "x").replace(",", "_") + f"_o{args.rain_over}_ev{args.evap}"
               if args.rain else "")
            + (f"_hcn{args.hcn}_nh{args.nh3}" if args.hcn or args.nh3 else "")
            + (f"_tar{args.tar_len}o{args.tar_oily:g}" if args.tar_len else "")
            + (f"_sun{args.sun_d:g}" if args.sun_d else "")
            + (f"_glow{args.dark_len}_h{args.glow_heat:g}_g{args.glow_gain:g}_pa{args.p_abs:g}_P{args.light_period}"
               if args.glow else "") + f"_s{args.seed}")
    if args.start:
        name = os.path.splitext(os.path.basename(args.start))[0].replace("_state", "") + "_then_" + name
        if args.reseed is not None:
            name += f"_rs{args.reseed}"
    os.makedirs(args.out, exist_ok=True)
    if args.start:
        w, _ = World.load(args.start)
        if args.reseed is not None:
            w.rng = np.random.default_rng(args.reseed)
        n = len(np.unique(w.mol[np.isin(w.type, TAILS)]))
        if args.shear:
            w.flow = shear_flow(args.box, args.shear)
    else:
        w, n = build(args)
    if args.hcn or args.nh3:
        add_nitrogen(w, args)
        w.step(200, max_force=50)   # new beads may sit close to old ones
        w.t -= 200
    nit = None
    if (w.type == HCN).any() or (w.type == NH3).any():
        w.repel_inside = tuple(sorted(set(w.repel_inside) | {HCN, NH3, TAR}))   # floppy chains: no folding through
        nit = Nitrogen(every=50, kT0=args.kT_earth if args.sun_d else w.kT,   # rates as calibrated, at this temperature
                       tar_len=args.tar_len or None)
        if args.tar_len:
            w.rules.set_oily("tar", args.tar_oily)
    end = w.t + args.steps
    rx = LipidJoin(p_join=args.p_join, no_cat=args.no_cat, p_break=args.p_break, length=args.length,
                   chains=(STAIL,) if args.fa_only else (TAIL, STAIL))
    prot = Protonation(args.ph, args.pka, every=rx.every, hand_over=args.sealed_h) if args.ph is not None else None
    shocks = PHShocks.parse(args.shocks)
    if args.shock_rate:
        lo, hi = map(float, args.shock_size.split(":"))
        shocks += [(s + w.t, j, hd) for s, j, hd in PHShocks.random(args.shock_rate, args.steps, args.seed, (lo, hi))]
        print("random shocks:", shocks, flush=True)
    ph_of = PHShocks(args.ph, shocks, args.relax) if args.ph is not None else None
    burn = Burn(every=rx.every, p_burn=args.p_burn, heat=args.heat) if args.fuel else None
    sun = Sun(args.sun_d, args.light_period, args.kT_earth) if args.sun_d else None
    pump = LightPump(every=rx.every, p_abs=args.p_abs, gain=args.pump_gain, period=args.light_period,
                     pka=args.pka, sun=sun) if args.pump else None
    rain = Rain(Rain.parse(args.rain), args.rain_over, args.evap) if args.rain else None
    glow = LightPump(every=rx.every, p_abs=args.p_abs, gain=args.glow_gain, period=args.light_period, pka=args.pka,
                     pigment=lambda w: dark_beads(w, args.dark_len), any_head=True,
                     heat=args.glow_heat, sun=sun) if args.glow else None
    if glow and args.glow_gain and prot is None:
        p.error("--glow-gain needs --ph")
    if args.sealed_h:
        prot.closed = Protonation.facing(w)
    print(f"{name}: {n} lipids, {w.n} beads", flush=True)
    film = None
    if args.video:
        film = Film(w, args.video, fps=args.fps)
    log = open(os.path.join(args.out, name + ".csv"), "w", newline="")
    writer = csv.DictWriter(log, FIELDS)
    writer.writeheader()
    t0 = time.time()
    hops = Hops()
    hops.check(w, rx)
    state = os.path.join(args.out, name + "_state.pkl")
    while w.t < end:
        for _ in range(args.every // rx.every):
            w.step(rx.every)
            if sun:
                w.kT = sun.kT(w.t)
            if prot:
                prot.ph = ph_of(w.t)
                prot.react(w)
            if pump:
                pump.react(w, prot.ph, prot.closed, prot)
            if glow:
                glow.react(w, prot.ph if prot else 7.0, prot.closed if prot else None, prot)
            if nit:
                nit.react(w)
            if burn:
                burn.react(w)
            rx.react(w)
            if rain:
                rain.apply(w)
            if args.feed and w.t % args.feed == 0:
                feed(w, args.batch, 0, clear=1.0, charged=args.charged, lengths=args.lengths, acid=args.acid,
                     alcohol=args.alcohol)
            hops.check(w, rx)
        row = measure(w, rx, prot)
        row["A_left"], row["B_left"] = hops.take()
        if prot:
            row["ph"] = round(prot.ph, 3)
        if burn:
            row["burns"] = burn.burns
        if pump:
            row["pumps"], row["light"] = pump.pumps, int(pump.light(w.t))
        if rain:
            row["water"] = round(rain.now, 4)
        if sun:
            row["kT"], row["sun"] = round(w.kT, 4), round(sun.brightness(w.t), 3)
        if nit:
            row.update(nitrogen_row(w, glow, args.dark_len))
        if args.sealed_h:
            prot.closed = Protonation.facing(w)   # who faces a sealed inside, refreshed every row
        if args.save_every and w.t % args.save_every == 0:
            w.save(state, args=vars(args))
        writer.writerow(row)
        log.flush()
        if film:
            film.grab(w, row, rx)
        if w.t % (10 * args.every) == 0:
            print(f"step {w.t} {row} ({time.time() - t0:.0f}s)", flush=True)
    log.close()
    w.save(state, args=vars(args))
    np.savez_compressed(os.path.join(args.out, name + "_final.npz"), pos=w.pos, type=w.type, mol=w.mol, L=w.L,
                        bonds=w.bonds, marker_A=w.markers[0], marker_B=w.markers[1])
    if film:
        film.close()


if __name__ == "__main__":
    main()
