"""Experiment 10: can a vesicle grow its wall and split? (mixed lipids + shear)

In the lab (Zhu & Szostak 2009) fatty-acid vesicles fed with more fatty acid
grow wall faster than they take in water: they stretch into long threads,
and gentle shear breaks the threads into small closed vesicles that keep
most of their contents. Our 3-tail lipids hold a ring but almost never move
between clusters; 2-tail ones move all the time but can't hold a ring alone.

Here one ring has a wall mixing both (--mix = fraction of 2-tail). The
2-tail lipids have short-chain tails ("stail") whose oiliness --hop-oily
sets (1 = same as a normal tail; higher = stickier: real fatty acids leave
a membrane less often the longer their chain). Loose
2-tail lipids drip in (one every --drip steps, like lipids washing in), and
the water slides in bands (shear flow u_x = shear sin(2 pi y / L)) that
stretch whatever floats in it. Sugar (--sugar) sits inside.

Logged every --log-every steps: lipids in the biggest cluster (3- and 2-tail),
loose lipids, how many closed water pockets (vesicles) there are and their
areas, sugar in each.

    python experiments/e10_grow.py --mix 0.5 --drip 500 --shear 0.05 --video runs/video/grow.mp4
    python experiments/e10_grow.py --mix 0 --drip 0 --shear 0     # control: plain 3-tail ring
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

import e08_compete as e08

FIELDS = ["step", "biggest", "biggest_3", "biggest_2", "loose", "clusters20", "pockets", "areas", "added"]


def ring(world, center, R, mix):
    """Bilayer ring like e08's (3-tail spacing); each lipid is 2-tail with chance mix."""
    first = world.n_mol
    length = 3 * world.rules.bond_r0
    gap = 2 * length + 0.6
    for r, out, mid in ((R, 0.0, R + length / 2), (R + gap, np.pi, R + gap - length / 2)):
        n = int(2 * np.pi * mid / 1.0)
        for k in range(n):
            a = 2 * np.pi * k / n
            if world.rng.uniform() < mix:
                add_lipid(world, *(center + r * np.array([np.cos(a), np.sin(a)])), a + out, n_tails=2, tail="stail")
            else:
                add_lipid(world, *(center + r * np.array([np.cos(a), np.sin(a)])), a + out, n_tails=3)
    return np.arange(first, world.n_mol)


def bag(world):
    """Centre and outer radius of the biggest lipid cluster (the bag)."""
    lab = cluster_labels(world)
    lipids = np.unique(world.mol[np.isin(world.type, TAILS)])
    big = np.bincount(lab[lipids]).argmax()
    beads = np.isin(world.mol, lipids[lab[lipids] == big])
    c = e08.centre(world, beads)
    d = world._wrapped()[beads] - c
    d -= world.L * np.round(d / world.L)
    return c, np.percentile(np.hypot(d[:, 0], d[:, 1]), 90)


def drip(world, near=0.0, clear=3.0, tries=5000):
    """One loose 2-tail lipid at a random spot clear of everything.

    near > 0: within `near` of the bag's outer edge (the bag sits in a lipid-rich spot)."""
    L = world.L
    tree = cKDTree(world._wrapped(), boxsize=L)
    if near:
        c, R = bag(world)
    for _ in range(tries):
        if near:
            r, a0 = R + clear + near * world.rng.uniform(), world.rng.uniform(0, 2 * np.pi)
            p = (c + r * np.array([np.cos(a0), np.sin(a0)])) % L
        else:
            p = world.rng.uniform(0, L, 2)
        a = world.rng.uniform(0, 2 * np.pi)
        beads = (p + np.outer(np.arange(3) * world.rules.bond_r0, [np.cos(a), np.sin(a)])) % L
        if tree.query(beads)[0].min() > clear:
            add_lipid(world, *p, a, n_tails=2, tail="stail")
            return 1
    return 0


def shear_flow(L, shear):
    k = 2 * np.pi / L

    def flow(pos):
        v = np.zeros_like(pos)
        v[:, 0] = shear * np.sin(k * pos[:, 1])
        return v

    return flow


def build(args):
    rules = Rules()
    rules.bend_k = 50
    rules.bend_r0 = 2 * rules.bond_r0
    rules.set_oily("stail", args.hop_oily)
    L = args.box
    world = World(box=L, rules=rules, kT=args.kT, seed=args.seed)
    c = np.array([L / 2, L / 4])   # where the shear is strongest (sin peaks at y = L/4)
    world.ring = ring(world, c, args.radius, args.mix)
    e08.fill(world, c, args.radius - 1.5, args.sugar)
    world.step(1000, max_force=50.0)
    world.t = 0
    world.n_start = world.n_mol
    if args.shear:
        world.flow = shear_flow(L, args.shear)
    return world


def measure(world):
    lab = cluster_labels(world)
    tail = np.isin(world.type, TAILS)
    lipids = np.unique(world.mol[tail])
    ntails = np.bincount(world.mol[tail], minlength=world.n_mol)[lipids]
    size = np.bincount(lab[lipids], minlength=lab.max() + 1)
    big = size.argmax()
    inbig = lab[lipids] == big
    areas, _ = water_pockets(world, centers=True)
    areas = sorted((int(a) for a in areas), reverse=True)
    return {"step": world.t, "biggest": int(size[big]), "biggest_3": int((inbig & (ntails == 3)).sum()),
            "biggest_2": int((inbig & (ntails == 2)).sum()), "loose": int((size[lab[lipids]] == 1).sum()),
            "clusters20": int((size >= 20).sum()), "pockets": len(areas), "areas": " ".join(map(str, areas)),
            "added": world.n_mol - world.n_start}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--radius", type=float, default=12.0)
    p.add_argument("--box", type=float, default=120.0)
    p.add_argument("--mix", type=float, default=0.5, help="fraction of 2-tail lipids in the ring")
    p.add_argument("--hop-oily", type=float, default=1.0, help="oiliness of the 2-tail lipids' tails")
    p.add_argument("--sugar", type=int, default=20)
    p.add_argument("--drip", type=int, default=500, help="add one loose 2-tail lipid every N steps (0 = never)")
    p.add_argument("--near", type=float, default=0.0,
                   help="drip food within this distance outside the bag (0 = anywhere)")
    p.add_argument("--shear", type=float, default=0.05, help="water speed of the shear bands (0 = still)")
    p.add_argument("--kT", type=float, default=0.6)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--steps", type=int, default=100000)
    p.add_argument("--log-every", type=int, default=1000)
    p.add_argument("--video", default=None, help="save a movie, one frame per 400 steps")
    p.add_argument("--fps", type=int, default=10)
    p.add_argument("--out", default="runs/e10_grow")
    args = p.parse_args()

    out = os.path.join(args.out, f"mix{args.mix}_o{args.hop_oily}_d{args.drip}_n{args.near}_sh{args.shear}_s{args.sugar}_seed{args.seed}")
    os.makedirs(out, exist_ok=True)
    world = build(args)
    movie = None
    if args.video:
        import matplotlib
        matplotlib.use("Agg")
        from alife.viz import Movie, View
        view = View(world, title=f"{args.mix:.0%} 2-tail wall (oily {args.hop_oily}), food every {args.drip} steps, shear {args.shear}")
        view.legend()
        caption = view.ax.text(0.01, 0.99, "", transform=view.ax.transAxes, fontsize=10, va="top", zorder=5)
        movie = Movie(view, args.video, fps=args.fps)
    log = open(os.path.join(out, "log.csv"), "w", newline="")
    writer = csv.DictWriter(log, FIELDS)
    writer.writeheader()
    start = time.time()
    while world.t < args.steps:
        world.step(100)
        if args.drip and world.t % args.drip == 0:
            drip(world, near=args.near)
        if movie and world.t % 400 == 0:
            s = measure(world)
            caption.set_text(f"step {world.t}   wall lipids {s['biggest']} (3-tail {s['biggest_3']}, "
                             f"2-tail {s['biggest_2']})   loose {s['loose']}   vesicles {s['pockets']} "
                             f"[{s['areas']}]")
            view.update()
            movie.grab()
        if world.t % args.log_every == 0:
            s = measure(world)
            writer.writerow(s)
            log.flush()
            if world.t % (10 * args.log_every) == 0:
                print(f"step {world.t}  wall {s['biggest']} (3t {s['biggest_3']} 2t {s['biggest_2']})  "
                      f"loose {s['loose']}  big clusters {s['clusters20']}  vesicles {s['pockets']} "
                      f"[{s['areas']}]  added {s['added']}  ({time.time() - start:.0f}s)", flush=True)
    log.close()
    if movie:
        movie.close()
        print(f"video: {args.video} ({movie.frames} frames)")


if __name__ == "__main__":
    main()
