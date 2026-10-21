"""Experiment 5: wet-dry cycles. Does a pond that dries and refills make vesicles?

On early Earth, ponds on volcanic land dried in the sun and refilled with
rain, again and again. Drying crowds lipids together into stacks; rain
spreads them out again. In the lab, this is how lipids are made into
vesicles (Deamer's wet-dry cycling).

Three ponds, same lipids, same temperature:
  wet      pond always full
  wetdry   pond dries to r_dry and refills, once per --period steps
  rain     wetdry, and refilling stirs the water (rain splash, inflow)

--sugar N dissolves N sugar beads in the pond: too big to cross a membrane,
so sugar trapped in a vesicle pushes its wall out (osmotic pressure).

--chem makes the sugar in the pond instead (formose, alife/chem/reactions.py):
rain brings --rain-carb formaldehyde per refill, Ca2+ (--ca) is the catalyst,
and a meteorite every ~--meteor-every steps drops --meteor-sugar ready-made
C5 sugars (0 = none). --leak lets small molecules cross membranes
(alife/chem/membrane.py): formaldehyde fast, ions very slowly, C5+ never.

Logged every --log-every steps: pond radius, clusters, every water pocket
with its lining score (heads facing in) and coverage, so real vesicles can
be counted and timed afterwards.

    python experiments/e05_wetdry.py --pond wetdry                # live window
    python experiments/e05_wetdry.py --pond wetdry --headless     # runs/e05_wetdry
    python experiments/e05_wetdry.py --plot                       # compare all logs
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

from alife.chem import TAIL, Rules, World, add_lipid
from alife.chem.analysis import lining, summary, water_pockets
from alife.chem.elements import INDEX
from alife.chem.membrane import Leak
from alife.chem.reactions import CARB, Formose, chain_lengths
from alife.env import Pond

FIELDS = ["step", "R", "clusters", "largest", "free_lipids", "vesicles", "real", "energy_per_lipid",
          "c1", "c2_4", "c5plus", "meteors", "pockets"]


def build(args):
    rules = Rules()
    rules.bend_k = args.bend
    rules.bend_r0 = 2 * rules.bond_r0
    rules.eps[TAIL, TAIL] = args.stick
    world = World(box=args.box, rules=rules, kT=args.kT, seed=args.seed)
    world.repel_inside = (CARB,)   # floppy sugar chains must not fold through themselves
    # lipids on grid cells inside the full pond, random direction (grid = no overlaps)
    L, length = args.box, args.tails * rules.bond_r0
    spacing = length + 1.2
    n = int(L / spacing)
    centers = (np.stack(np.meshgrid(np.arange(n), np.arange(n)), -1).reshape(-1, 2) + 0.5) * spacing
    centers += (L - n * spacing) / 2
    centers = centers[np.hypot(*(centers - L / 2).T) < args.r_wet - length]
    if len(centers) < args.lipids:
        raise ValueError(f"pond of radius {args.r_wet} fits only {len(centers)} lipids")
    rng = world.rng
    for c in rng.permutation(centers)[:args.lipids]:
        angle = rng.uniform(0, 2 * np.pi)
        start = c - 0.5 * length * np.array([np.cos(angle), np.sin(angle)])
        add_lipid(world, *start, angle, n_tails=args.tails)
    add_sugar(world, args.sugar, args.r_wet - 2)
    if args.chem:
        for name, count in (("Ca", args.ca), ("Cl", 2 * args.ca)):
            for p in free_spots(world, count, L / 2, args.r_wet - 2):
                world.add_molecule([p], [name])
    world.step(200, max_force=50.0)
    world.t = 0
    pond = Pond(L, args.r_wet, args.r_dry, 0 if args.pond == "wet" else args.period,
                rain=args.rain if args.pond == "rain" else 0.0, seed=args.seed + 1000)
    return world, pond


def free_spots(world, count, center, radius, min_gap=1.0, tries=500):
    """Up to `count` random points within `radius` of center, min_gap from any bead and each other."""
    L = world.L
    pts = []
    tree = cKDTree(world._wrapped(), boxsize=L) if world.n else None
    for _ in range(count):
        for _ in range(tries):
            r, a = radius * np.sqrt(world.rng.uniform()), world.rng.uniform(0, 2 * np.pi)
            p = (center + r * np.array([np.cos(a), np.sin(a)])) % L
            d = p - np.array(pts) if pts else np.zeros((0, 2))
            d -= L * np.round(d / L)
            if (tree is None or tree.query(p)[0] >= min_gap) and not (np.hypot(*d.T) < min_gap).any():
                pts.append(p)
                break
    return pts


class Chemistry:
    """Formose in the pond: rain brings formaldehyde, meteorites bring sugar."""

    def __init__(self, args):
        self.formose = Formose(every=50)
        self.leak = Leak(every=50) if args.leak else None
        self.rain_carb, self.period = args.rain_carb, args.period
        self.meteor_every, self.meteor_sugar = args.meteor_every, args.meteor_sugar
        self.rained = 0.0
        self.meteors = 0

    def __call__(self, world, pond, k):
        """Every chunk of k steps: react, then add rain and maybe a meteorite."""
        self.formose.react(world)
        if self.leak is not None:
            self.leak.leak(world)
        L = world.L
        R = pond.radius(world.t)
        if self.period and self.rain_carb:
            # formaldehyde arrives with the refill, at the pace the pond refills
            # (refill_speed integrates to period/pi over one cycle)
            self.rained += self.rain_carb * pond.refill_speed(world.t) * k * np.pi / self.period
            n = int(self.rained)
            self.rained -= n
            for p in free_spots(world, n, L / 2, R - 1):
                world.add_molecule([p], ["carb"])
        if self.meteor_sugar and world.rng.uniform() < k / self.meteor_every:
            # a small rock lands in the pond: a burst of ready-made C5 sugars near one spot
            r, a = max(R - 10, 0) * np.sqrt(world.rng.uniform()), world.rng.uniform(0, 2 * np.pi)
            hit = L / 2 + r * np.array([np.cos(a), np.sin(a)])
            added = 0
            for p in free_spots(world, self.meteor_sugar * 20, hit, 10.0):   # crowded pond: try many spots
                ang = world.rng.uniform(0, 2 * np.pi)
                chain = p + np.outer(np.arange(5) - 2, [np.cos(ang), np.sin(ang)])
                if cKDTree(world._wrapped(), boxsize=L).query(chain % L)[0].min() < 0.9:
                    continue   # chain would overlap something: try the next spot
                world.add_molecule(chain, ["carb"] * 5, bonds=[(i, i + 1) for i in range(4)])
                added += 1
                if added == self.meteor_sugar:
                    break
            self.meteors += 1


def add_sugar(world, count, radius, min_gap=1.0, tries=500):
    """Dissolve sugar beads at free spots inside the pond."""
    from scipy.spatial import cKDTree
    L = world.L
    for k in range(count):
        tree = cKDTree(world._wrapped(), boxsize=L)
        for _ in range(tries):
            r, a = radius * np.sqrt(world.rng.uniform()), world.rng.uniform(0, 2 * np.pi)
            p = (L / 2 + r * np.array([np.cos(a), np.sin(a)])) % L
            if tree.query(p)[0] >= min_gap:
                break
        else:
            raise ValueError(f"no room for sugar {k + 1}")
        world.add_molecule([p], ["sugar"])


def advance(world, pond, n, chunk=50, chem=None):
    """Step n, moving the shore (and new rain swirls) every chunk steps."""
    done = 0
    while done < n:
        k = min(chunk, n - done)
        if chem is not None:
            chem(world, pond, k)   # before the shore: new beads need the shore too
        world.external = pond.shore_at(world.t, world.frozen)
        world.flow = pond.flow_at(world.t)
        world.step(k)
        done += k


def info(world, pond, min_lining, min_cover, chem=None):
    s = summary(world)
    counts = chain_lengths(world)
    carbons = counts * np.arange(len(counts))
    s["c1"], s["c2_4"], s["c5plus"] = int(carbons[1:2].sum()), int(carbons[2:5].sum()), int(carbons[5:].sum())
    s["meteors"] = chem.meteors if chem else 0
    s["R"] = round(pond.radius(world.t), 2)
    areas, where, masks = water_pockets(world, centers=True, masks=True)
    m_ = masks[0].shape[0] if masks else 1
    sugar = world._wrapped()[world.type == INDEX["sugar"]].astype(int) % m_
    # carbons in big sugars (chains of 5+, too big to leak out: they make the pressure) and small ones
    carb = np.flatnonzero(world.type == CARB)
    size = np.bincount(world.mol[carb], minlength=world.n_mol)[world.mol[carb]]
    big = world._wrapped()[carb[size >= 5]].astype(int) % m_
    small = world._wrapped()[carb[size < 5]].astype(int) % m_
    ions = world._wrapped()[np.isin(world.type, [INDEX[n] for n in ("Na", "Ca", "Cl")])].astype(int) % m_
    parts, real = [], 0
    for a, (x, y), m in zip(areas, where, masks):
        sc, cov = lining(world, m)
        real += sc >= min_lining and cov >= min_cover
        # inside = on a pocket grid cell (cell = 1 unit). trapped = sugar beads + big-sugar carbons
        trapped = int(m[sugar[:, 0], sugar[:, 1]].sum()) + int(m[big[:, 0], big[:, 1]].sum())
        parts.append(f"{x:.1f},{y:.1f},{a},{sc:.2f},{cov:.2f},{trapped},{int(m[small[:, 0], small[:, 1]].sum())},"
                     f"{int(m[ions[:, 0], ions[:, 1]].sum())}")
    s["pockets"] = ";".join(parts)
    s["real"] = real
    return s


def save_frame(world, path):
    np.savez_compressed(path, pos=world._wrapped(), type=world.type, mol=world.mol,
                        bonds=world.bonds, L=world.L, t=world.t)


def draw(view, world, pond, s):
    import matplotlib.patches as mpatches
    if getattr(view, "shore", None) is None:
        view.shore = mpatches.Circle((world.L / 2, world.L / 2), pond.r_wet, fill=False,
                                     ec="#5F5E5A", lw=1.5, ls="--", zorder=2)
        view.ax.add_patch(view.shore)
    view.shore.set_radius(s["R"])
    return view.update(s)


def run_headless(world, pond, args, chem=None):
    import matplotlib
    matplotlib.use("Agg")
    from alife.viz import View

    out = os.path.join(args.out, f"{args.name or args.pond}_seed{args.seed}")
    os.makedirs(out, exist_ok=True)
    view = View(world, title=args.pond)
    log = open(os.path.join(out, "log.csv"), "w", newline="")
    writer = csv.DictWriter(log, FIELDS, extrasaction="ignore")
    writer.writeheader()
    start = time.time()
    snap_every = args.steps // args.snapshots
    while world.t < args.steps:
        advance(world, pond, args.log_every, chem=chem)
        s = info(world, pond, args.min_lining, args.min_cover, chem)
        writer.writerow(s)
        log.flush()
        if s["vesicles"]:
            save_frame(world, os.path.join(out, f"frame_{world.t:06d}.npz"))
        if world.t % snap_every == 0:
            draw(view, world, pond, s)
            view.save(os.path.join(out, f"step_{world.t:06d}.png"), s)
            print(f"step {world.t} R {s['R']} clusters {s['clusters']} largest {s['largest']} "
                  f"pockets {s['vesicles']} real {s['real']}  carbons C1 {s['c1']} C2-4 {s['c2_4']} "
                  f"C5+ {s['c5plus']} meteors {s['meteors']}  ({time.time() - start:.0f}s)", flush=True)
    log.close()


def run_live(world, pond, args, chem=None):
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    from alife.viz import View

    view = View(world, title=args.pond)

    def frame(_):
        advance(world, pond, args.steps_per_frame, chem=chem)
        return draw(view, world, pond, info(world, pond, args.min_lining, args.min_cover, chem))

    _anim = FuncAnimation(view.fig, frame, interval=1, blit=False, cache_frame_data=False)
    plt.show()


def lifetimes(rows, L, match=6.0, every=1):
    """Follow real vesicles (pockets passing the lining test) between logs.

    rows: [(step, [(x, y), ...]), ...]. Returns list of lifetimes in steps.
    """
    tracks, done = [], []
    for step, ps in rows:
        new, taken = [], set()
        for tr in tracks:
            x0, y0 = tr[-1]
            best = None
            for k, (x, y) in enumerate(ps):
                d = np.array([x - x0, y - y0])
                d -= L * np.round(d / L)
                if k not in taken and np.hypot(*d) < match:
                    best = k
                    break
            if best is None:
                done.append(len(tr) * every)
            else:
                taken.add(best)
                tr.append(ps[best])
                new.append(tr)
        new += [[p] for k, p in enumerate(ps) if k not in taken]
        tracks = new
    return done + [len(tr) * every for tr in tracks]


def plot(args):
    """Compare wet and wetdry logs: pond size, clusters, real vesicles."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"wet": "#378ADD", "wetdry": "#D85A30", "rain": "#1D9E75"}
    panels = ["R", "largest", "free_lipids", "real"]
    fig, axes = plt.subplots(len(panels), 1, figsize=(9, 9), sharex=True)
    totals = {}
    for path in sorted(glob.glob(os.path.join(args.out, "*_seed*", "log.csv"))):
        name = os.path.basename(os.path.dirname(path))
        pond = name.split("_seed")[0]
        rows = list(csv.DictReader(open(path)))
        if len(rows) < 2:
            continue
        t = [int(r["step"]) for r in rows]
        for ax, key in zip(axes, panels):
            ax.plot(t, [float(r[key]) for r in rows], color=colors.get(pond, "k"), lw=0.8, alpha=0.6)
        every = t[1] - t[0]
        real_xy = []
        for r in rows:
            ps = [tuple(float(v) for v in p.split(",")) for p in r["pockets"].split(";") if p]
            real_xy.append((int(r["step"]), [(x, y) for x, y, a, sc, cov, *_ in ps
                                            if sc >= args.min_lining and cov >= args.min_cover]))
        life = lifetimes(real_xy, args.box, every=every)
        tot = totals.setdefault(pond, {"runs": 0, "vesicles": 0, "life": [], "logs_with": 0, "logs": 0})
        tot["runs"] += 1
        tot["vesicles"] += len(life)
        tot["life"] += life
        tot["logs_with"] += sum(1 for _, ps in real_xy if ps)
        tot["logs"] += len(real_xy)
    for pond, tot in totals.items():
        life = np.array(tot["life"]) if tot["life"] else np.array([0])
        print(f"{pond:7s} {tot['runs']} runs: {tot['vesicles']} real vesicles, "
              f"present {100 * tot['logs_with'] / tot['logs']:.1f}% of the time, "
              f"median life {np.median(life):.0f}, longest {life.max()} steps")
    for ax, key in zip(axes, panels):
        ax.set_ylabel(key)
    axes[0].set_title("blue: always wet   orange: wet-dry   green: wet-dry + rain")
    axes[-1].set_xlabel("step")
    path = os.path.join(args.out, "compare.png")
    fig.savefig(path, dpi=100, bbox_inches="tight")
    print("saved", path)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pond", choices=["wet", "wetdry", "rain"], default="wetdry")
    p.add_argument("--lipids", type=int, default=300)
    p.add_argument("--tails", type=int, default=3)
    p.add_argument("--box", type=float, default=96.0)
    p.add_argument("--r-wet", type=float, default=46.0, help="pond radius when full")
    p.add_argument("--r-dry", type=float, default=22.0, help="pond radius at driest (crowded)")
    p.add_argument("--period", type=int, default=20000, help="steps per dry + refill cycle")
    p.add_argument("--rain", type=float, default=1.0, help="top swirl speed while refilling")
    p.add_argument("--sugar", type=int, default=0, help="sugar beads dissolved in the pond")
    p.add_argument("--chem", action="store_true", help="make sugar in the pond (formose)")
    p.add_argument("--rain-carb", type=int, default=40, help="formaldehyde brought by each refill")
    p.add_argument("--ca", type=int, default=20, help="Ca2+ catalyst ions (with 2 Cl- each), --chem only")
    p.add_argument("--meteor-every", type=int, default=50000, help="mean steps between meteorites")
    p.add_argument("--meteor-sugar", type=int, default=0, help="C5 sugars per meteorite (0 = none)")
    p.add_argument("--leak", action="store_true", help="small molecules leak through membranes (--chem only)")
    p.add_argument("--name", default=None, help="run folder name (default: --pond)")
    p.add_argument("--kT", type=float, default=0.6)
    p.add_argument("--stick", type=float, default=1.0)
    p.add_argument("--bend", type=float, default=50.0)
    p.add_argument("--min-lining", type=float, default=0.8)
    p.add_argument("--min-cover", type=float, default=0.9)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--steps", type=int, default=200000)
    p.add_argument("--log-every", type=int, default=200)
    p.add_argument("--snapshots", type=int, default=40)
    p.add_argument("--steps-per-frame", type=int, default=50)
    p.add_argument("--headless", action="store_true")
    p.add_argument("--plot", action="store_true")
    p.add_argument("--out", default="runs/e05_wetdry")
    args = p.parse_args()

    if args.plot:
        plot(args)
        return
    world, pond = build(args)
    chem = Chemistry(args) if args.chem else None
    print(f"{args.pond}: {args.lipids} lipids, {args.sugar} sugar, pond r {args.r_wet}"
          + (f" -> {args.r_dry} every {args.period} steps" if args.pond != "wet" else "")
          + (f", rain {args.rain}" if args.pond == "rain" else ""))
    if args.headless:
        run_headless(world, pond, args, chem)
    else:
        run_live(world, pond, args, chem)


if __name__ == "__main__":
    main()
