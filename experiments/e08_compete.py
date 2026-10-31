"""Experiment 8: do protocells compete for lipids? (growth by stealing)

Two ring vesicles side by side in water, same size. Lipids now and then pop
out of a membrane and rejoin one (a few are always loose). Ring A holds
--sugar-a sugar beads, ring B --sugar-b. Sugar can't leave and pulls water
in, so a sugar-filled vesicle is swollen and its membrane stretched. In the
lab (Chen, Roberts & Szostak 2004), swollen fatty-acid vesicles grew by
taking lipids from relaxed neighbours, which shrank. Same here?

Membrane food: --feed N loose lipids at the start, and/or one new lipid
every --drip steps at a random free spot away from everything (a slow
steady supply, like lipids made elsewhere washing in).

Logged every --log-every steps: lipids in each ring's membrane, loose
lipids, and each ring's water area.

    python experiments/e08_compete.py --sugar-a 40 --drip 2000 --seed 0   # runs/e08_compete/a40_b0_f0_d2000_seed0
    python experiments/e08_compete.py --sugar-a 0 --sugar-b 0    # control: both empty
    python experiments/e08_compete.py --video runs/video/compete.mp4
    python experiments/e08_compete.py --mix 0.5 --hop-oily 1.8 --drip 500   # walls that trade lipids

--mix: fraction of each wall made of 2-tail lipids with short tails ("stail",
oiliness --hop-oily), like fatty acids: they leave and join membranes often
(e10 showed 3-tail walls hardly ever trade). Then the drip adds 2-tail lipids too.
"""
import argparse
import csv
import os
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from alife.chem import TAILS, Rules, World, add_lipid
from alife.chem.analysis import cluster_labels, water_pockets
from alife.chem.elements import INDEX

FIELDS = ["step", "lipids_a", "lipids_b", "loose", "merged", "area_a", "area_b", "added", "two_a", "two_b"]


def ring(world, center, R, tails=3, mix=None):
    """Bilayer ring: inner heads point in, outer heads out. Returns its molecule ids.
    mix: each lipid is 2-tail (short "stail" tails) with this chance; spacing stays 3-tail.

    Each layer gets one lipid per unit of length at its middle radius: counted
    at the head radius, the outer layer's tails were crammed (0.75 apart) and
    some starts blew up.
    """
    first = world.n_mol
    length = tails * world.rules.bond_r0
    gap = 2 * length + 0.6
    for r, out, mid in ((R, 0.0, R + length / 2), (R + gap, np.pi, R + gap - length / 2)):
        n = int(2 * np.pi * mid / 1.0)
        for k in range(n):
            a = 2 * np.pi * k / n
            xy = center + r * np.array([np.cos(a), np.sin(a)])
            if mix is not None and world.rng.uniform() < mix:
                add_lipid(world, *xy, a + out, n_tails=2, tail="stail")
            else:
                add_lipid(world, *xy, a + out, n_tails=tails)
    return np.arange(first, world.n_mol)


def fill(world, center, radius, count, name="sugar", gap=1.1):
    """Put count beads at free spots within radius of center."""
    L = world.L
    for k in range(count):
        tree = cKDTree(world._wrapped(), boxsize=L)
        for _ in range(5000):
            r, a = radius * np.sqrt(world.rng.uniform()), world.rng.uniform(0, 2 * np.pi)
            p = (center + r * np.array([np.cos(a), np.sin(a)])) % L
            if tree.query(p)[0] > gap:
                world.add_molecule([p], [name])
                break
        else:
            raise RuntimeError(f"no room for {name} {k + 1}")


def feed(world, count, clear=3.0, tries=5000, tails=3):
    """Loose lipids at random spots with nothing within `clear` of any of their beads.
    tails=2: short-tail ("stail") lipids, like fatty acids."""
    L = world.L
    for k in range(count):
        tree = cKDTree(world._wrapped(), boxsize=L)
        for _ in range(tries):
            p, a = world.rng.uniform(0, L, 2), world.rng.uniform(0, 2 * np.pi)
            beads = (p + np.outer(np.arange(4) * world.rules.bond_r0, [np.cos(a), np.sin(a)])) % L
            if tree.query(beads)[0].min() > clear:
                add_lipid(world, *(p % L), a, n_tails=tails, tail="stail" if tails == 2 else INDEX["tail"])
                break
        else:
            return k   # no room left
    return count


def build(args):
    rules = Rules()
    rules.bend_k = 50
    rules.bend_r0 = 2 * rules.bond_r0
    if args.mix:
        rules.set_oily("stail", args.hop_oily)
    L = args.box
    world = World(box=L, rules=rules, kT=args.kT, seed=args.seed)
    ca = np.array([L / 4, L / 2])
    cb = np.array([3 * L / 4, L / 2])
    mix = args.mix or None
    world.ring_a, world.ring_b = ring(world, ca, args.radius, mix=mix), ring(world, cb, args.radius, mix=mix)
    world.centers = (ca, cb)
    fill(world, ca, args.radius - 1.5, args.sugar_a)
    fill(world, cb, args.radius - 1.5, args.sugar_b)
    feed(world, args.feed)
    world.step(1000, max_force=50.0)
    world.t = 0
    world.radius, world.n_start = args.radius, world.n_mol
    return world


def centre(world, beads):
    """Mean position of beads in the periodic box (circular mean)."""
    th = world._wrapped()[beads] / world.L * 2 * np.pi
    return (np.arctan2(np.sin(th).mean(0), np.cos(th).mean(0)) % (2 * np.pi)) / (2 * np.pi) * world.L


def measure(world):
    """Lipids in each ring's cluster (the cluster holding most of its first lipids), loose lipids, areas."""
    lab = cluster_labels(world)
    tail = np.isin(world.type, TAILS)
    lipids = np.unique(world.mol[tail])
    two = np.bincount(world.mol[tail], minlength=world.n_mol)[lipids] == 2
    size = np.bincount(lab[lipids], minlength=lab.max() + 1)
    home_a = np.bincount(lab[world.ring_a]).argmax()
    home_b = np.bincount(lab[world.ring_b]).argmax()
    s = {"step": world.t, "lipids_a": int(size[home_a]), "lipids_b": int(size[home_b]),
         "loose": int((size[lab[lipids]] == 1).sum()), "merged": int(home_a == home_b),
         "added": world.n_mol - world.n_start,
         "two_a": int((two & (lab[lipids] == home_a)).sum()), "two_b": int((two & (lab[lipids] == home_b)).sum())}
    # water area of the biggest pocket near the middle of each ring's lipids
    areas, where = water_pockets(world, centers=True)
    for key, home in (("area_a", home_a), ("area_b", home_b)):
        c = centre(world, np.isin(world.mol, lipids[lab[lipids] == home]))
        best = 0
        for a, xy in zip(areas, where):
            d = np.array(xy) - c
            d -= world.L * np.round(d / world.L)
            if np.hypot(*d) < world.radius and a > best:
                best = a
        s[key] = int(best)
    return s


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--radius", type=float, default=12.0)
    p.add_argument("--box", type=float, default=120.0)
    p.add_argument("--sugar-a", type=int, default=40)
    p.add_argument("--sugar-b", type=int, default=0)
    p.add_argument("--feed", type=int, default=0, help="loose lipids in the water at the start")
    p.add_argument("--drip", type=int, default=0, help="add one loose lipid every N steps (0 = never)")
    p.add_argument("--mix", type=float, default=0.0, help="fraction of 2-tail lipids in each wall (and drip 2-tail)")
    p.add_argument("--hop-oily", type=float, default=1.8, help="oiliness of the 2-tail lipids' tails")
    p.add_argument("--kT", type=float, default=0.6)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--steps", type=int, default=200000)
    p.add_argument("--log-every", type=int, default=1000)
    p.add_argument("--video", default=None, help="save a movie, one frame per 400 steps")
    p.add_argument("--fps", type=int, default=10)
    p.add_argument("--out", default="runs/e08_compete")
    args = p.parse_args()

    mix = f"_mix{args.mix:g}_o{args.hop_oily:g}" if args.mix else ""
    out = os.path.join(args.out, f"a{args.sugar_a}_b{args.sugar_b}_f{args.feed}_d{args.drip}{mix}_seed{args.seed}")
    os.makedirs(out, exist_ok=True)
    world = build(args)
    movie = None
    if args.video:
        import matplotlib
        matplotlib.use("Agg")
        from alife.viz import Movie, View
        view = View(world, title=f"A: {args.sugar_a} sugar   B: {args.sugar_b} sugar")
        view.legend()
        caption = view.ax.text(0.01, 0.99, "", transform=view.ax.transAxes, fontsize=10, va="top", zorder=5)
        movie = Movie(view, args.video, fps=args.fps)
    log = open(os.path.join(out, "log.csv"), "w", newline="")
    writer = csv.DictWriter(log, FIELDS)
    writer.writeheader()
    start = time.time()
    while world.t < args.steps:
        world.step(200)
        if args.drip and world.t % args.drip == 0:
            feed(world, 1, tails=2 if args.mix else 3)
        if movie and world.t % 400 == 0:
            s = measure(world)
            caption.set_text(f"lipids in A (left) {s['lipids_a']}   B (right) {s['lipids_b']}   loose {s['loose']}"
                             + ("   MERGED" if s["merged"] else ""))
            view.update()
            movie.grab()
        if world.t % args.log_every == 0:
            s = measure(world)
            writer.writerow(s)
            log.flush()
            if world.t % (20 * args.log_every) == 0:
                print(f"step {world.t}  A {s['lipids_a']}  B {s['lipids_b']}  loose {s['loose']}  "
                      f"merged {s['merged']}  area A {s['area_a']} B {s['area_b']}  ({time.time() - start:.0f}s)",
                      flush=True)
    log.close()
    if movie:
        movie.close()
        print(f"video: {args.video} ({movie.frames} frames)")


if __name__ == "__main__":
    main()
