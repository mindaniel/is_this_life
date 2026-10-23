"""Experiment 7: is a vesicle a reactor? (food in, sugar made inside, sugar stays)

One ring vesicle. Inside: a few sugar seeds (C2) and Ca2+ (the catalyst).
Outside: formaldehyde (food), no catalyst. With leaky membranes
(alife/chem/membrane.py), formaldehyde slips in, Ca2+ helps join it onto
the seeds (formose), and the grown sugar is too big to leak back out. If
that works, the inside piles up sugar from outside food: a first, very
simple metabolism.

Conditions:
  reactor      leaky membrane, Ca2+ inside
  --sealed     nothing crosses: food can't get in (control)
  --ca none    leaky, no catalyst anywhere (control)
  --ca outside leaky, Ca2+ in the open water instead (is the bag needed?)

    python experiments/e07_reactor.py --seed 0                       # runs/e07_reactor/reactor_seed0
    python experiments/e07_reactor.py --sealed --seed 0
    python experiments/e07_reactor.py --video runs/video/reactor.mp4
"""
import argparse
import csv
import os
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from alife.chem import Rules, World, add_lipid
from alife.chem.analysis import water_pockets
from alife.chem.elements import INDEX
from alife.chem.membrane import Leak
from alife.chem.reactions import CARB, Formose

FIELDS = ["step", "in_c1", "in_c2_4", "in_c5plus", "out_c1", "out_c2plus", "in_ca", "area", "hops", "joins"]


def build(args):
    R = args.radius
    rules = Rules()
    rules.bend_k = 50
    rules.bend_r0 = 2 * rules.bond_r0
    L = 2 * R + args.margin
    world = World(box=L, rules=rules, kT=args.kT, seed=args.seed)
    world.repel_inside = (CARB,)
    c = np.array([L / 2, L / 2])
    # bilayer ring: inner layer heads point in, outer layer heads point out
    for r, out in ((R, 0.0), (R + 6.6, np.pi)):
        n = int(2 * np.pi * r / 1.0)
        for k in range(n):
            a = 2 * np.pi * k / n
            add_lipid(world, *(c + r * np.array([np.cos(a), np.sin(a)])), a + out, n_tails=3)
    world.n_lipid_beads = world.n
    placed = []

    def spot(lo, hi, gap=1.1):
        """Free point at distance lo..hi from the centre (uniform by area)."""
        tree = cKDTree(world._wrapped(), boxsize=L)
        for _ in range(20000):
            r = np.sqrt(world.rng.uniform(lo ** 2, hi ** 2))
            a = world.rng.uniform(0, 2 * np.pi)
            p = (c + r * np.array([np.cos(a), np.sin(a)])) % L
            d = p - np.array(placed) if placed else np.zeros((0, 2))
            d -= L * np.round(d / L)
            if tree.query(p)[0] > gap and not (np.hypot(*d.T) < gap).any():
                placed.append(p)
                return p
        raise RuntimeError("no room")

    inside, outside = (0.0, R - 1.5), (R + 6.6 + 2.5, 0.5 * L - 0.5)
    for _ in range(args.seeds):   # C2 sugar seeds (skip formose's slow first step)
        p = spot(*inside)
        a = world.rng.uniform(0, 2 * np.pi)
        q = p + np.array([np.cos(a), np.sin(a)])
        if np.hypot(*(q - c)) > R - 1.5 or cKDTree(world._wrapped(), boxsize=L).query(q % L)[0] < 1.1:
            q = p - np.array([np.cos(a), np.sin(a)])
        placed.append(q)
        world.add_molecule([p, q], ["carb", "carb"], bonds=[(0, 1)])
    if args.ca != "none":
        where = inside if args.ca == "inside" else outside
        for name, count in (("Ca", args.n_ca), ("Cl", 2 * args.n_ca)):
            for _ in range(count):
                world.add_molecule([spot(*where)], [name])
    for _ in range(args.food):
        world.add_molecule([spot(*outside)], ["carb"])
    world.step(200, max_force=50.0)
    world.t = 0
    return world


def measure(world, formose, leak):
    """Carbons inside / outside the ring, by chain length."""
    L = world.L
    p = world._wrapped()
    lip = np.zeros(world.n, dtype=bool)
    lip[:world.n_lipid_beads] = True
    # ring centre: circular mean of lipid beads (periodic box); inside = closer than mid-membrane
    th = p[lip] / L * 2 * np.pi
    c = (np.arctan2(np.sin(th).mean(0), np.cos(th).mean(0)) % (2 * np.pi)) / (2 * np.pi) * L
    d = p - c
    d -= L * np.round(d / L)
    r = np.hypot(*d.T)
    inside = r < np.median(r[lip])
    carb = world.type == CARB
    size = np.bincount(world.mol, minlength=world.n_mol)[world.mol]
    areas = water_pockets(world)
    s = {"step": world.t,
         "in_c1": int((inside & carb & (size == 1)).sum()),
         "in_c2_4": int((inside & carb & (size >= 2) & (size <= 4)).sum()),
         "in_c5plus": int((inside & carb & (size >= 5)).sum()),
         "out_c1": int((~inside & carb & (size == 1)).sum()),
         "out_c2plus": int((~inside & carb & (size >= 2)).sum()),
         "in_ca": int((inside & (world.type == INDEX["Ca"])).sum()),
         "area": int(max(areas)) if len(areas) else 0,
         "hops": leak.hops if leak else 0, "joins": formose.joins}
    return s


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--radius", type=float, default=8.0, help="inner radius of the ring")
    p.add_argument("--margin", type=float, default=40.0, help="box = 2 radius + margin")
    p.add_argument("--food", type=int, default=150, help="formaldehyde outside")
    p.add_argument("--seeds", type=int, default=3, help="C2 sugars inside at the start")
    p.add_argument("--ca", choices=["inside", "outside", "none"], default="inside")
    p.add_argument("--n-ca", type=int, default=4)
    p.add_argument("--sealed", action="store_true", help="no leaks: nothing crosses the membrane")
    p.add_argument("--kT", type=float, default=0.6)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--steps", type=int, default=60000)
    p.add_argument("--log-every", type=int, default=1000)
    p.add_argument("--video", default=None, help="save a movie, one frame per 200 steps")
    p.add_argument("--fps", type=int, default=10)
    p.add_argument("--name", default=None)
    p.add_argument("--out", default="runs/e07_reactor")
    args = p.parse_args()

    name = args.name or ("sealed" if args.sealed else "reactor" if args.ca == "inside" else f"ca_{args.ca}")
    out = os.path.join(args.out, f"{name}_seed{args.seed}")
    os.makedirs(out, exist_ok=True)
    world = build(args)
    formose = Formose(every=50)
    leak = None if args.sealed else Leak()
    movie = None
    if args.video:
        import matplotlib
        matplotlib.use("Agg")
        from alife.viz import Movie, View
        view = View(world, title=name)
        view.legend()
        caption = view.ax.text(0.01, 0.99, "", transform=view.ax.transAxes, fontsize=9, va="top", zorder=5)
        movie = Movie(view, args.video, fps=args.fps)
    log = open(os.path.join(out, "log.csv"), "w", newline="")
    writer = csv.DictWriter(log, FIELDS)
    writer.writeheader()
    start = time.time()
    frame_every = 200
    while world.t < args.steps:
        world.step(50)
        formose.react(world)
        if leak:
            leak.leak(world)
        if movie and world.t % frame_every == 0:
            s = measure(world, formose, leak)
            caption.set_text(f"inside: food C1 {s['in_c1']}  small sugar C {s['in_c2_4']}  big sugar C {s['in_c5plus']}\n"
                             f"outside: food C1 {s['out_c1']}  sugar C {s['out_c2plus']}")
            view.update()
            movie.grab()
        if world.t % args.log_every == 0:
            s = measure(world, formose, leak)
            writer.writerow(s)
            log.flush()
            if world.t % (10 * args.log_every) == 0:
                print(f"step {world.t}  inside C1 {s['in_c1']} C2-4 {s['in_c2_4']} C5+ {s['in_c5plus']}  "
                      f"outside C1 {s['out_c1']} C2+ {s['out_c2plus']}  Ca in {s['in_ca']}  area {s['area']}  "
                      f"hops {s['hops']}  ({time.time() - start:.0f}s)", flush=True)
    log.close()
    if movie:
        movie.close()
        print(f"video: {args.video} ({movie.frames} frames)")


if __name__ == "__main__":
    main()
