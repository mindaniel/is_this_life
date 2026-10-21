"""Experiment 3: what do other elements do to membranes?

Early-Earth lipids were likely fatty acids: a head with charge -1 (COO-).
Every negative head needs a positive ion nearby to keep the water neutral.
Worlds:

  neutral   uncharged heads (same as e01, baseline)
  na        fatty acids + Na+ (sodium, +1)
  ca        fatty acids + Na+ + extra Ca2+ / Cl-  (calcium, +2)
  mg        fatty acids + Na+ + extra Mg2+ / Cl-  (magnesium, +2)

Real labs: fatty-acid membranes are wrecked by too much Ca2+ (soap scum).
Mg2+ is gentler because it keeps its water shell. Do our rules show that?

    python experiments/e03_ions.py --world ca            # live window
    python experiments/e03_ions.py --world ca --headless
    python experiments/e03_ions.py --plot
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

from alife.chem import TAIL, TYPE_NAMES, Rules, World, add_ions, fill_with_lipids, net_charge
from alife.chem.analysis import summary

FIELDS = ["step", "clusters", "largest", "free_lipids", "vesicles", "bridges", "energy_per_lipid"]
FA = TYPE_NAMES.index("fa_head")


def build(args):
    rules = Rules()
    rules.bend_k = args.bend
    rules.bend_r0 = 2 * rules.bond_r0
    rules.eps[TAIL, TAIL] = args.stick
    rules.debye = args.debye
    world = World(box=args.box, rules=rules, kT=args.kT, seed=args.seed)
    head = "head" if args.world == "neutral" else "fa_head"
    fill_with_lipids(world, args.lipids, n_tails=args.tails, head=head)
    if args.world != "neutral":
        add_ions(world, "Na", args.lipids)            # one Na+ per COO-
    if args.world in ("ca", "mg"):
        metal = "Ca" if args.world == "ca" else "Mg"
        add_ions(world, metal, args.metal)
        add_ions(world, "Cl", 2 * args.metal)         # keep total charge zero
    assert abs(net_charge(world)) < 1e-9, net_charge(world)
    world.step(200, max_force=50.0)
    world.t = 0
    return world


def bridges(world, contact=1.6):
    """+2 ions touching heads of two or more different lipids: a salt bridge."""
    metal = np.flatnonzero(np.isin(world.type, [TYPE_NAMES.index("Ca"), TYPE_NAMES.index("Mg")]))
    heads = np.flatnonzero(world.type == FA)
    if not len(metal) or not len(heads):
        return 0
    tree = cKDTree(world._wrapped()[heads], boxsize=world.L)
    near = tree.query_ball_point(world._wrapped()[metal], contact)
    return int(sum(len(set(world.mol[heads[h]] for h in hs)) >= 2 for hs in near))


def info(world):
    s = summary(world)
    s["bridges"] = bridges(world)
    return s


def run_headless(world, args):
    import matplotlib
    matplotlib.use("Agg")
    from alife.viz import View

    out = os.path.join(args.out, f"{args.world}_seed{args.seed}")
    os.makedirs(out, exist_ok=True)
    view = View(world, title=args.world)
    log = open(os.path.join(out, "log.csv"), "w", newline="")
    writer = csv.DictWriter(log, FIELDS, extrasaction="ignore")
    writer.writeheader()
    start = time.time()
    snap_every = args.steps // args.snapshots
    view.save(os.path.join(out, "step_000000.png"), info(world))
    while world.t < args.steps:
        world.step(args.log_every)
        s = info(world)
        writer.writerow(s)
        log.flush()
        if world.t % snap_every == 0:
            view.save(os.path.join(out, f"step_{world.t:06d}.png"), s)
            print(f"{s}  ({time.time() - start:.0f}s)", flush=True)
    log.close()


def run_live(world, args):
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    from alife.viz import View

    view = View(world, title=args.world)

    def frame(_):
        world.step(args.steps_per_frame)
        return view.update(info(world))

    _anim = FuncAnimation(view.fig, frame, interval=1, blit=False, cache_frame_data=False)
    plt.show()


def plot(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"neutral": "#888780", "na": "#378ADD", "ca": "#444441", "mg": "#639922"}
    panels = ["clusters", "largest", "free_lipids", "bridges", "energy_per_lipid"]
    fig, axes = plt.subplots(len(panels), 1, figsize=(9, 11), sharex=True)
    for path in sorted(glob.glob(os.path.join(args.out, "*", "log.csv"))):
        name = os.path.basename(os.path.dirname(path))
        rows = list(csv.DictReader(open(path)))
        t = [int(r["step"]) for r in rows]
        for ax, key in zip(axes, panels):
            ax.plot(t, [float(r[key]) for r in rows], color=colors.get(name.split("_seed")[0], "k"),
                    lw=1, alpha=0.8, label=name)
        last = rows[-1]
        print(f"{name:14s} clusters {last['clusters']:>3}  largest {last['largest']:>3}  free {last['free_lipids']:>3}  "
              f"vesicles {last['vesicles']}  bridges {last['bridges']:>3}")
    for ax, key in zip(axes, panels):
        ax.set_ylabel(key)
    axes[0].legend(fontsize=8, ncol=4)
    axes[-1].set_xlabel("step")
    path = os.path.join(args.out, "compare.png")
    fig.savefig(path, dpi=100, bbox_inches="tight")
    print("saved", path)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--world", choices=["neutral", "na", "ca", "mg"], default="ca")
    p.add_argument("--lipids", type=int, default=200)
    p.add_argument("--metal", type=int, default=40, help="number of Ca2+ or Mg2+ ions")
    p.add_argument("--debye", type=float, default=2.0, help="screening length (less = saltier)")
    p.add_argument("--tails", type=int, default=3)
    p.add_argument("--box", type=float, default=64.0)
    p.add_argument("--kT", type=float, default=0.6)
    p.add_argument("--stick", type=float, default=1.0)
    p.add_argument("--bend", type=float, default=50.0)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--steps", type=int, default=100000)
    p.add_argument("--log-every", type=int, default=1000)
    p.add_argument("--snapshots", type=int, default=5)
    p.add_argument("--steps-per-frame", type=int, default=50)
    p.add_argument("--headless", action="store_true")
    p.add_argument("--plot", action="store_true")
    p.add_argument("--out", default="runs/e03_ions")
    args = p.parse_args()

    if args.plot:
        return plot(args)
    world = build(args)
    print(f"{args.world}: {world.n} beads, net charge {net_charge(world):+.0f}, debye {args.debye}")
    if args.headless:
        run_headless(world, args)
    else:
        run_live(world, args)


if __name__ == "__main__":
    main()
