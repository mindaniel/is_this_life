"""Experiment 2: does a changing planet help membranes organise?

Same lipids as e01, three worlds:
  steady   constant temperature
  cycle    temperature swings up and down (seasons)
  shocks   seasons + rare asteroid impacts (heat spike + blast)

The last part of every run is calm (steady base temperature), so all three
are judged at the same temperature. Measured: clusters, largest cluster,
vesicles (sealed water pockets), energy per lipid (lower = better packed).

    python experiments/e02_environment.py --env shocks          # live window
    python experiments/e02_environment.py --env shocks --headless
    python experiments/e02_environment.py --plot                # compare logs
"""
import argparse
import csv
import glob
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from alife.chem import TAIL, Rules, World, fill_with_lipids
from alife.chem.analysis import summary
from alife.env import Environment

FIELDS = ["step", "kT", "light", "clusters", "largest", "free_lipids", "vesicles", "energy_per_lipid"]


def build(args):
    rules = Rules()
    rules.bend_k = args.bend
    rules.bend_r0 = 2 * rules.bond_r0
    rules.eps[TAIL, TAIL] = args.stick
    world = World(box=args.box, rules=rules, kT=args.kT, seed=args.seed)
    fill_with_lipids(world, args.lipids, n_tails=args.tails)
    world.step(200, max_force=50.0)
    world.t = 0

    calm_after = int(args.steps * (1 - args.calm))
    env = Environment(base_kT=args.kT, calm_after=calm_after, seed=args.seed + 1000)
    if args.env in ("cycle", "shocks"):
        env.amp, env.period = args.amp, args.period
        env.day_period, env.day_heat = args.day_period, args.day_heat
        env.current = args.current
    if args.env == "shocks":
        env.shock_every = args.shock_every
    return world, env


def info(world, env):
    s = summary(world)
    s["kT"] = round(float(world.kT), 3)
    if env.day_heat and not env.calm(world.t):
        s["light"] = round(float(env.light), 3)
    return s


def run_headless(world, env, args):
    import matplotlib
    matplotlib.use("Agg")
    from alife.viz import View

    tag = args.env
    if args.env != "steady":
        tag += "_sun" if args.day_heat else ""
        tag += "_flow" if args.current else ""
    out = os.path.join(args.out, f"{tag}_seed{args.seed}")
    os.makedirs(out, exist_ok=True)
    view = View(world, title=args.env)
    log = open(os.path.join(out, "log.csv"), "w", newline="")
    writer = csv.DictWriter(log, FIELDS, extrasaction="ignore")
    writer.writeheader()
    start = time.time()
    snap_every = args.steps // args.snapshots
    n_events = 0
    while world.t < args.steps:
        env.run(world, args.log_every)
        s = info(world, env)
        writer.writerow(s)
        log.flush()
        if len(env.events) > n_events:
            n_events = len(env.events)
            print(f"  asteroid at step {env.events[-1][0]}", flush=True)
        if world.t % snap_every == 0:
            view.save(os.path.join(out, f"step_{world.t:06d}.png"), s)
            print(f"{s}  ({time.time() - start:.0f}s)", flush=True)
    log.close()


def run_live(world, env, args):
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    from alife.viz import View

    view = View(world, title=args.env)

    def frame(_):
        env.run(world, args.steps_per_frame)
        return view.update(info(world, env))

    _anim = FuncAnimation(view.fig, frame, interval=1, blit=False, cache_frame_data=False)
    plt.show()


def plot(args):
    """Overlay the logs of every run in args.out."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"steady": "#888780", "cycle": "#378ADD", "shocks": "#D85A30",
              "cycle_sun": "#1D9E75", "shocks_sun": "#7F77DD", "shocks_sun_flow": "#D4537E"}
    panels = ["kT", "clusters", "largest", "vesicles", "energy_per_lipid"]
    fig, axes = plt.subplots(len(panels), 1, figsize=(9, 11), sharex=True)
    for path in sorted(glob.glob(os.path.join(args.out, "*", "log.csv"))):
        name = os.path.basename(os.path.dirname(path))
        env = name.split("_seed")[0]
        rows = list(csv.DictReader(open(path)))
        t = [int(r["step"]) for r in rows]
        for ax, key in zip(axes, panels):
            ax.plot(t, [float(r[key]) for r in rows], color=colors.get(env, "k"),
                    lw=1, alpha=0.8, label=name)
        last = rows[-1]
        print(f"{name:16s} clusters {last['clusters']:>3}  largest {last['largest']:>3}  "
              f"vesicles {last['vesicles']}  energy {float(last['energy_per_lipid']):.3f}")
    for ax, key in zip(axes, panels):
        ax.set_ylabel(key)
    axes[0].legend(fontsize=8, ncol=3)
    axes[-1].set_xlabel("step")
    path = os.path.join(args.out, "compare.png")
    fig.savefig(path, dpi=100, bbox_inches="tight")
    print("saved", path)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--env", choices=["steady", "cycle", "shocks"], default="shocks")
    p.add_argument("--lipids", type=int, default=200)
    p.add_argument("--tails", type=int, default=3)
    p.add_argument("--box", type=float, default=64.0)
    p.add_argument("--kT", type=float, default=0.6, help="base temperature")
    p.add_argument("--stick", type=float, default=1.0)
    p.add_argument("--bend", type=float, default=50.0)
    p.add_argument("--amp", type=float, default=0.15, help="season swing in kT")
    p.add_argument("--period", type=int, default=20000, help="steps per season cycle")
    p.add_argument("--day-period", type=int, default=2000, help="steps per day/night (Earth spin)")
    p.add_argument("--day-heat", type=float, default=0.0, help="extra kT at noon from sunlight")
    p.add_argument("--current", type=float, default=0.0, help="water swirl speed (sun-driven)")
    p.add_argument("--shock-every", type=int, default=40000, help="mean steps between asteroids")
    p.add_argument("--calm", type=float, default=0.15, help="fraction of run at the end held steady")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--steps", type=int, default=200000)
    p.add_argument("--log-every", type=int, default=1000)
    p.add_argument("--snapshots", type=int, default=10)
    p.add_argument("--steps-per-frame", type=int, default=100)
    p.add_argument("--headless", action="store_true")
    p.add_argument("--plot", action="store_true", help="compare finished runs")
    p.add_argument("--out", default="runs/e02_environment")
    args = p.parse_args()

    if args.plot:
        return plot(args)
    world, env = build(args)
    print(f"{args.env}: {args.lipids} lipids, box {world.L}, base kT {args.kT}")
    if args.headless:
        run_headless(world, env, args)
    else:
        run_live(world, env, args)


if __name__ == "__main__":
    main()
