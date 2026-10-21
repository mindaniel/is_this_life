"""Experiment 4: one big connected world that is not the same everywhere.

A 3 x 3 grid of zones in one box. Matter moves freely between zones.
  - zones differ: a hot vent, cold corners, shaded and sunny places
  - rock walls make pools: one almost closed (narrow gap), one open (wide gap),
    and a ridge with a strait through it
  - weather: seasons, day/night sun, drifting currents
  - asteroids are LOCAL: heat + blast near the impact, the rest survives

Questions: where do membranes and vesicles form? Do survivors refill a hit zone?
Do semi-closed pools keep lipids concentrated?

    python experiments/e04_world.py              # live window (slow: 9x e01)
    python experiments/e04_world.py --headless   # snapshots + log in runs/e04_world
"""
import argparse
import csv
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from alife.chem import TAIL, Rules, World, add_ions, add_rock_wall, fill_with_lipids
from alife.chem.analysis import summary
from alife.env import Environment, Terrain

FIELDS = ["step", "kT", "clusters", "largest", "free_lipids", "vesicles", "energy_per_lipid", "zones"]

# zone map, row 0 = bottom. heat = extra kT, sun = share of sunlight
HEAT = [[0.00, -0.10, 0.00],
        [0.20, 0.00, -0.10],     # hot vent on the left middle
        [0.00, 0.05, -0.15]]     # cold top right
SUN = [[1.0, 1.0, 0.6],
       [0.3, 1.0, 1.0],          # vent is deep: little sun
       [1.0, 0.6, 1.0]]


def square_pool(world, cx, cy, side, gap):
    """Four rock walls around (cx, cy). gap = open fraction in the right wall."""
    h = side / 2
    a, b = 0.5 - gap / 2, 0.5 + gap / 2
    add_rock_wall(world, (cx - h, cy - h), (cx + h, cy - h))
    add_rock_wall(world, (cx - h, cy + h), (cx + h, cy + h))
    add_rock_wall(world, (cx - h, cy - h), (cx - h, cy + h))
    add_rock_wall(world, (cx + h, cy - h), (cx + h, cy + h), gaps=[(a, b)])


def build(args):
    rules = Rules()
    rules.bend_k = args.bend
    rules.bend_r0 = 2 * rules.bond_r0
    rules.eps[TAIL, TAIL] = args.stick
    world = World(box=args.box, rules=rules, kT=args.kT, seed=args.seed)
    z = args.box / 3
    square_pool(world, 0.5 * z, 2.5 * z, 0.65 * z, gap=0.08)   # almost closed pool, top left
    square_pool(world, 2.5 * z, 0.5 * z, 0.65 * z, gap=0.40)   # open pool, bottom right
    add_rock_wall(world, (1.5 * z, 1.1 * z), (1.5 * z, 2.3 * z), gaps=[(0.4, 0.6)])   # ridge with strait

    head = "fa_head" if args.fatty else "head"
    fill_with_lipids(world, args.lipids, n_tails=args.tails, head=head)
    if args.fatty:
        add_ions(world, "Na", args.lipids)
    world.step(200, max_force=50.0)
    world.t = 0

    env = Environment(base_kT=args.kT, amp=args.amp, period=args.period,
                      day_period=args.day_period, day_heat=args.day_heat,
                      current=args.current, eddies=3,
                      shock_every=args.shock_every, shock_heat=0.5, blast_radius=args.blast,
                      terrain=Terrain(args.box, HEAT, SUN), shock_local=True, seed=args.seed + 1000)
    return world, env


def zone_counts(world):
    """Lipids (counted by head bead) in each zone, bottom-left to top-right."""
    heads = world._wrapped()[np.isin(world.mol, np.unique(world.mol[world.type == TAIL]))
                             & (world.type != TAIL)]
    ij = np.minimum((heads / (world.L / 3)).astype(int), 2)
    counts = np.zeros((3, 3), dtype=int)
    np.add.at(counts, (ij[:, 1], ij[:, 0]), 1)
    return counts


def info(world, env):
    s = summary(world)
    s["kT"] = round(float(world.kT), 3)
    s["zones"] = " ".join(str(c) for c in zone_counts(world).ravel())
    return s


def draw(view, world, env, s):
    view.show_field(env.kT_field_at(world.t, world.L), vmin=env.base_kT - 0.3, vmax=env.base_kT + 0.6)
    return view.update(s)


def run_headless(world, env, args):
    import matplotlib
    matplotlib.use("Agg")
    from alife.viz import View

    out = os.path.join(args.out, f"seed{args.seed}")
    os.makedirs(out, exist_ok=True)
    view = View(world, title="world", size=10)
    log = open(os.path.join(out, "log.csv"), "w", newline="")
    writer = csv.DictWriter(log, FIELDS, extrasaction="ignore")
    writer.writeheader()
    start = time.time()
    snap_every = args.steps // args.snapshots
    s = info(world, env)
    draw(view, world, env, s)
    view.save(os.path.join(out, "step_000000.png"), s)
    n_events = 0
    while world.t < args.steps:
        env.run(world, args.log_every)
        s = info(world, env)
        writer.writerow(s)
        log.flush()
        if len(env.events) > n_events:
            n_events = len(env.events)
            t, x, y = env.events[-1]
            print(f"  asteroid at step {t}, ({x:.0f}, {y:.0f})", flush=True)
        if world.t % snap_every == 0:
            draw(view, world, env, s)
            view.save(os.path.join(out, f"step_{world.t:06d}.png"), s)
            print(f"{s}  ({time.time() - start:.0f}s)", flush=True)
    log.close()


def run_live(world, env, args):
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    from alife.viz import View

    view = View(world, title="world", size=9)

    def frame(_):
        env.run(world, args.steps_per_frame)
        return draw(view, world, env, info(world, env))

    _anim = FuncAnimation(view.fig, frame, interval=1, blit=False, cache_frame_data=False)
    plt.show()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lipids", type=int, default=1800)
    p.add_argument("--fatty", action="store_true", help="fatty acid heads (-1) + Na+ ions")
    p.add_argument("--tails", type=int, default=3)
    p.add_argument("--box", type=float, default=192.0)
    p.add_argument("--kT", type=float, default=0.6)
    p.add_argument("--stick", type=float, default=1.0)
    p.add_argument("--bend", type=float, default=50.0)
    p.add_argument("--amp", type=float, default=0.1, help="season swing")
    p.add_argument("--period", type=int, default=20000)
    p.add_argument("--day-period", type=int, default=2000)
    p.add_argument("--day-heat", type=float, default=0.1)
    p.add_argument("--current", type=float, default=0.3)
    p.add_argument("--shock-every", type=int, default=30000, help="mean steps between asteroids")
    p.add_argument("--blast", type=float, default=20.0, help="asteroid blast and heat radius")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--steps", type=int, default=100000)
    p.add_argument("--log-every", type=int, default=1000)
    p.add_argument("--snapshots", type=int, default=10)
    p.add_argument("--steps-per-frame", type=int, default=50)
    p.add_argument("--headless", action="store_true")
    p.add_argument("--out", default="runs/e04_world")
    args = p.parse_args()

    world, env = build(args)
    print(f"world: {world.n} beads ({args.lipids} lipids), box {world.L}")
    if args.headless:
        run_headless(world, env, args)
    else:
        run_live(world, env, args)


if __name__ == "__main__":
    main()
