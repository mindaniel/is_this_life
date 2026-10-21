"""Experiment 6: sugar from formaldehyde (formose), in a small box of water.

Formaldehyde (CH2O, one carb bead) rains into ponds from lightning and UV in
the air. With Ca2+ as catalyst, carbons join into chains: 2 = glycolaldehyde,
5 = ribose-sized. See alife/chem/reactions.py for the rules.

Checks that the rules behave like real formose:
  - no Ca: almost nothing happens          --ca 0
  - hotter is faster                       --kT 0.7
  - slow start, then sugar multiplies      compare --seed-c2 0 with --seed-c2 5
    (a little glycolaldehyde skips the slow start)

    python experiments/e06_formose.py                          # live window
    python experiments/e06_formose.py --headless --name base   # runs/e06_formose/base_seed0
    python experiments/e06_formose.py --plot                   # compare all runs
"""
import argparse
import csv
import glob
import os
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from alife.chem import World, add_ions
from alife.chem.reactions import CARB, Formose, chain_lengths

MAX_LEN = 8   # log chains 1..7 separately, 8 = 8 or longer
FIELDS = ["step", "joins", "splits"] + [f"c{n}" for n in range(1, MAX_LEN + 1)] + ["sugar_carbons"]


def build(args):
    world = World(box=args.box, kT=args.kT, seed=args.seed)
    world.repel_inside = (CARB,)   # sugar chains are floppy: their beads must not overlap
    add_ions(world, "carb", args.formaldehyde)
    for _ in range(args.seed_c2):   # a little glycolaldehyde (C2) at free spots
        tree = cKDTree(world._wrapped(), boxsize=args.box)
        while True:
            two = world.rng.uniform(0, args.box, 2) + np.array([[0, 0], [1.0, 0]])
            if tree.query(two % args.box)[0].min() >= 1.0:
                break
        world.add_molecule(two, ["carb", "carb"], bonds=[(0, 1)])
    if args.ca:
        add_ions(world, "Ca", args.ca)
        add_ions(world, "Cl", 2 * args.ca)
    world.step(200, max_force=50.0)
    world.t = 0
    chem = Formose(every=args.every)
    return world, chem


def advance(world, chem, n):
    done = 0
    while done < n:
        k = min(chem.every, n - done)
        world.step(k)
        chem.react(world)
        done += k


def info(world, chem):
    counts = chain_lengths(world)
    row = {"step": world.t, "joins": chem.joins, "splits": chem.splits}
    for n in range(1, MAX_LEN):
        row[f"c{n}"] = int(counts[n]) if n < len(counts) else 0
    row[f"c{MAX_LEN}"] = int(counts[MAX_LEN:].sum())
    # carbons sitting in sugar-sized chains (5 or more): what osmosis would see
    lengths = np.arange(len(counts))
    row["sugar_carbons"] = int((counts * lengths)[5:].sum())
    return row


def run_headless(world, chem, args):
    out = os.path.join(args.out, f"{args.name}_seed{args.seed}")
    os.makedirs(out, exist_ok=True)
    log = open(os.path.join(out, "log.csv"), "w", newline="")
    writer = csv.DictWriter(log, FIELDS)
    writer.writeheader()
    start = time.time()
    while world.t < args.steps:
        advance(world, chem, args.log_every)
        row = info(world, chem)
        writer.writerow(row)
        log.flush()
        if world.t % (args.steps // 10) == 0:
            print(f"step {world.t}  " + " ".join(f"C{n}:{row[f'c{n}']}" for n in range(1, MAX_LEN + 1))
                  + f"  sugar carbons {row['sugar_carbons']}  ({time.time() - start:.0f}s)", flush=True)
    log.close()


def run_live(world, chem, args):
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    from alife.viz import View

    view = View(world, title="formose")

    def frame(_):
        advance(world, chem, args.steps_per_frame)
        row = info(world, chem)
        s = {"clusters": 0, "largest": 0, "free_lipids": 0, "vesicles": 0}
        view.ax.set_title(f"formose  step {world.t}  " + " ".join(f"C{n}:{row[f'c{n}']}" for n in range(1, 7)),
                          fontsize=9, loc="left")
        return view.update(s)

    _anim = FuncAnimation(view.fig, frame, interval=1, blit=False, cache_frame_data=False)
    plt.show()


def plot(args):
    """Share of carbon in formaldehyde (C1), small sugars (C2-4) and big sugars (C5+)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    runs = {}
    for path in sorted(glob.glob(os.path.join(args.out, "*_seed*", "log.csv"))):
        name = os.path.basename(os.path.dirname(path)).split("_seed")[0]
        rows = list(csv.DictReader(open(path)))
        if rows:
            runs.setdefault(name, []).append(rows)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
    titles = ["formaldehyde C1", "small sugars C2-C4", "big sugars C5+"]
    groups = [range(1, 2), range(2, 5), range(5, MAX_LEN + 1)]
    colors = plt.cm.tab10.colors
    for c, (name, logs) in enumerate(runs.items()):
        for k, rows in enumerate(logs):
            t = np.array([int(r["step"]) for r in rows])
            carbons = np.array([[int(r[f"c{n}"]) * n for n in range(1, MAX_LEN + 1)] for r in rows], dtype=float)
            total = carbons.sum(axis=1, keepdims=True)
            for ax, g in zip(axes, groups):
                ax.plot(t, carbons[:, [n - 1 for n in g]].sum(axis=1) / total[:, 0], color=colors[c % 10],
                        lw=0.8, alpha=0.7, label=name if k == 0 else None)
        last = np.array([[int(rows[-1][f"c{n}"]) * n for n in range(1, MAX_LEN + 1)] for rows in logs], dtype=float)
        share = last / last.sum(axis=1, keepdims=True)
        print(f"{name:12s} {len(logs)} runs, carbon at end: C1 {share[:, 0].mean():.0%}  "
              f"C2-4 {share[:, 1:4].sum(1).mean():.0%}  C5+ {share[:, 4:].sum(1).mean():.0%}")
    for ax, t in zip(axes, titles):
        ax.set_title(t)
        ax.set_xlabel("step")
    axes[0].set_ylabel("share of all carbon")
    axes[0].legend(fontsize=8)
    path = os.path.join(args.out, "compare.png")
    fig.savefig(path, dpi=100, bbox_inches="tight")
    print("saved", path)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--formaldehyde", type=int, default=400)
    p.add_argument("--ca", type=int, default=20, help="Ca2+ ions (catalyst), with 2 Cl- each")
    p.add_argument("--seed-c2", type=int, default=0, help="glycolaldehyde (C2) added at the start")
    p.add_argument("--box", type=float, default=40.0)
    p.add_argument("--kT", type=float, default=0.6)
    p.add_argument("--every", type=int, default=50, help="steps between reaction checks")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--steps", type=int, default=100000)
    p.add_argument("--log-every", type=int, default=500)
    p.add_argument("--steps-per-frame", type=int, default=100)
    p.add_argument("--name", default="base")
    p.add_argument("--headless", action="store_true")
    p.add_argument("--plot", action="store_true")
    p.add_argument("--out", default="runs/e06_formose")
    args = p.parse_args()

    if args.plot:
        plot(args)
        return
    world, chem = build(args)
    print(f"{args.name}: {args.formaldehyde} formaldehyde, {args.ca} Ca, {args.seed_c2} C2, kT {args.kT}")
    if args.headless:
        run_headless(world, chem, args)
    else:
        run_live(world, chem, args)


if __name__ == "__main__":
    main()
