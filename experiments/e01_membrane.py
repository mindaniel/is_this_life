"""Experiment 1: do membranes self-assemble from oily lipids in water?

Starts with lipids scattered randomly. Only rules: tails stick to tails,
everything has a hard core, bonds hold molecules together, heat jiggles.
Watch for micelles (small clumps), bilayers (double lines), vesicles (rings).
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from alife.chem import TAIL, Rules, World, fill_with_lipids
from alife.chem.analysis import summary


def build(args):
    rules = Rules()
    rules.bend_k = args.bend
    rules.bend_r0 = 2 * rules.bond_r0
    rules.eps[TAIL, TAIL] = args.stick
    world = World(box=args.box, rules=rules, kT=args.kT, seed=args.seed)
    fill_with_lipids(world, args.lipids, n_tails=args.tails)
    # random start has some tight contacts: relax gently before real physics
    world.step(200, max_force=50.0)
    world.t = 0
    return world


def run_headless(world, args):
    import matplotlib
    matplotlib.use("Agg")
    from alife.viz import View

    os.makedirs(args.out, exist_ok=True)
    view = View(world, title="membrane")
    view.save(os.path.join(args.out, "step_000000.png"), summary(world))
    start = time.time()
    chunk = args.steps // args.snapshots
    for _ in range(args.snapshots):
        world.step(chunk)
        info = summary(world)
        view.save(os.path.join(args.out, f"step_{world.t:06d}.png"), info)
        print(f"{info}  ({time.time() - start:.0f}s)", flush=True)


def run_live(world, args):
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    from alife.viz import View

    view = View(world, title="membrane")

    def frame(_):
        world.step(args.steps_per_frame)
        return view.update(summary(world))

    _anim = FuncAnimation(view.fig, frame, interval=1, blit=False, cache_frame_data=False)
    plt.show()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--lipids", type=int, default=200)
    p.add_argument("--tails", type=int, default=3, help="tail beads per lipid")
    p.add_argument("--box", type=float, default=64.0)
    p.add_argument("--kT", type=float, default=0.6, help="temperature (heat noise)")
    p.add_argument("--stick", type=float, default=1.0, help="tail-tail attraction strength")
    p.add_argument("--bend", type=float, default=50.0, help="lipid stiffness")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--headless", action="store_true", help="save PNGs instead of a window")
    p.add_argument("--steps", type=int, default=150000)
    p.add_argument("--snapshots", type=int, default=8)
    p.add_argument("--steps-per-frame", type=int, default=50)
    p.add_argument("--out", default="runs/e01_membrane")
    args = p.parse_args()

    world = build(args)
    print(f"{args.lipids} lipids, {world.n} beads, box {world.L}, kT {args.kT}")
    if args.headless:
        run_headless(world, args)
    else:
        run_live(world, args)


if __name__ == "__main__":
    main()
