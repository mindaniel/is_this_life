"""Experiment 9: a landscape of different ponds under weather.

A patch of volcanic land with several ponds, none alike:
  - sizes and shapes differ (wobbly shorelines: oval, lobed, ragged)
  - some are shallow and dry out further, some deep and barely shrink
  - each has its own bottom:
      plain   nothing special
      clay    charged mineral flakes (clay speeds up vesicle making in the lab:
              Hanczyc, Fujikawa & Szostak 2003); Na+ balances their charge
      rocky   rock islands that split the current and give surfaces
      hot     a warm spring in one spot (kT +0.15 at its centre): everything
              there shakes harder and reactions run faster
The land slopes: each pond sits at its own height, and when a pond is
nearly full it overflows down a stream into its nearest lower neighbour
(--stream-speed), carrying water and whatever floats in it. Streams run in
wet seasons and dry up in summer. Weather (four seasons; wind from one
slowly wandering direction pushing water through each pond; storms whose
middle whirls the ponds round) is on unless --no-weather. --chem: formose with rain,
Ca2+ and meteorites, --leak, --crust as in e05.

The landscape (positions, shapes, kinds) comes from its own random numbers
(--land-seed), so one landscape can be run with many --seed histories.

Logged every --log-every steps: e05's columns plus, per pond, lipids in it
and vesicles in it ("pond_lipids", "pond_vesicles", space separated).

    python experiments/e09_landscape.py --seed 0 --video runs/video/landscape.mp4
    python experiments/e09_landscape.py --chem --leak --crust --seed 0
"""
import argparse
import csv
import os
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))

from alife.chem import TAIL, Rules, World, add_lipid
from alife.chem.elements import INDEX
from alife.chem.reactions import CARB
from alife.env import Crust, Pond, Weather
import e05_wetdry as e05

FIELDS = e05.FIELDS[:-1] + ["pond_lipids", "pond_vesicles", "pockets"]
KINDS = ("plain", "clay", "rocky", "hot")
LAND = "#E9DFC9"
KIND_COLORS = {"plain": "#5F5E5A", "clay": "#A0522D", "rocky": "#2C2C2A", "hot": "#D85A30"}


def landscape(args):
    """Box size, pond centres, sizes, depths, shapes, kinds, heights and streams (own random numbers)."""
    if args.layout == "basin":
        return basin(args)
    rng = np.random.default_rng(args.land_seed)
    L = np.ceil(np.sqrt(args.ponds)) * args.spacing
    side = int(np.ceil(np.sqrt(args.ponds)))
    cells = (np.stack(np.meshgrid(np.arange(side), np.arange(side)), -1).reshape(-1, 2) + 0.5) * args.spacing
    centers = cells[:args.ponds] + rng.uniform(-0.08, 0.08, (args.ponds, 2)) * args.spacing
    scales = rng.uniform(0.8, 1.15, args.ponds)
    depths = rng.uniform(0.6, 1.0, args.ponds)
    shapes = []
    for _ in range(args.ponds):
        waves = [(2, rng.uniform(0.05, 0.2), rng.uniform(0, 2 * np.pi)),
                 (3, rng.uniform(0.0, 0.12), rng.uniform(0, 2 * np.pi)),
                 (int(rng.integers(4, 7)), rng.uniform(0.0, 0.05), rng.uniform(0, 2 * np.pi))]
        shapes.append(waves)
    kinds = [KINDS[i % len(KINDS)] for i in range(args.ponds)] if args.kinds is None else args.kinds.split(",")
    heights = rng.uniform(size=args.ponds)   # drawn last: the rest of the landscape stays as before
    streams = []
    for i in range(args.ponds):   # overflow runs to the nearest lower pond
        lower = [j for j in range(args.ponds) if heights[j] < heights[i]]
        if lower:
            gap = [np.hypot(*((centers[j] - centers[i] + L / 2) % L - L / 2)) for j in lower]
            streams.append((i, lower[int(np.argmin(gap))]))
    return L, centers, scales, depths, shapes, kinds, heights, streams


def wobbles(rng, n):
    return [[(2, rng.uniform(0.05, 0.2), rng.uniform(0, 2 * np.pi)),
             (3, rng.uniform(0.0, 0.12), rng.uniform(0, 2 * np.pi)),
             (int(rng.integers(4, 7)), rng.uniform(0.0, 0.05), rng.uniform(0, 2 * np.pi))] for _ in range(n)]


def ground(rng, L, waves=6):
    """Smooth rolling ground: height(pos) from a few long waves that fit the (wrapping) box."""
    ks = [(1, 0), (0, 1), (1, 1), (1, -1), (2, 1), (1, 2), (2, -1), (-1, 2)]
    pick = rng.choice(len(ks), waves, replace=False)
    k = 2 * np.pi / L * np.array([ks[i] for i in pick], dtype=float)
    amp = rng.uniform(0.5, 1.0, waves) / (1 + np.hypot(k[:, 0], k[:, 1]) * L / (2 * np.pi))
    phase = rng.uniform(0, 2 * np.pi, waves)
    return lambda pos: (amp * np.cos(np.atleast_2d(pos) @ k.T + phase)).sum(-1)


def basin(args):
    """A river basin: ponds scattered over rolling ground; water runs downhill.

    The lowest pond is a big lake that hardly dries. Higher ponds are smaller
    and dry deeper. Each pond overflows (when nearly full) into the lowest
    nearest lower pond within reach (stream at most --reach across dry land,
    not through another pond), so streams merge on their way down like
    branches into a trunk. Water can't run over a ridge higher than the pond
    it leaves. A pond with no lower pond in reach has no outlet: it stays
    alone except when floods join it to a neighbour.
    """
    rng = np.random.default_rng(args.land_seed)
    n, R0 = args.ponds, args.r_wet
    L = args.spacing * np.sqrt(n) * 1.2
    wave = ground(rng, L)
    # shift the ground so its lowest spot is the middle of the box: the lake goes there, and
    # every pond stays clear of the edges (the world wraps around; the basin should not)
    g = np.linspace(0, L, 64, endpoint=False)
    grid = np.stack(np.meshgrid(g, g), -1).reshape(-1, 2)
    low = grid[np.argmin(wave(grid))]
    height = lambda pos: wave(np.atleast_2d(pos) + low - L / 2)
    edge = args.reach / 2   # dry land between any pond and the edge (no streams or floods across it)
    for attempt in range(200):
        scales = np.r_[args.lake, rng.uniform(0.55, 1.1, n - 1)]
        centers = [np.full(2, L / 2)]
        for i in range(1, n):
            ext = scales[i] * R0 * 1.4 + edge   # widest shoreline (wobbles add up to < 0.4)
            for _ in range(2000):
                c = rng.uniform(ext, L - ext, 2)
                if all(np.hypot(*((c - q + L / 2) % L - L / 2)) > (scales[i] + scales[j]) * R0 - 4
                       for j, q in enumerate(centers)):
                    centers.append(c)
                    break
            else:
                break
        if len(centers) == n:
            break
    else:
        raise ValueError("ponds don't fit: lower --ponds or raise --spacing")
    centers = np.array(centers)
    heights = height(centers)
    depths = np.r_[0.25, rng.uniform(0.6, 1.0, n - 1)]   # the lake hardly dries
    shapes = wobbles(rng, n)
    kinds = ["plain"] + [KINDS[k] for k in rng.integers(0, len(KINDS), n - 1)] if args.kinds is None \
        else args.kinds.split(",")
    # streams: overflow to the nearest lower pond within reach, not through another pond
    streams = []
    for i in range(1, n):
        best = None
        for j in range(n):
            if heights[j] >= heights[i]:
                continue
            ab = (centers[j] - centers[i] + L / 2) % L - L / 2
            gap = np.hypot(*ab) - (scales[i] + scales[j]) * R0   # dry land the stream crosses
            if gap > args.reach:
                continue
            e = ab / np.hypot(*ab)
            # water can't climb: the dry land between the two shores must be lower than the upper pond
            u = np.linspace(scales[i] * R0, np.hypot(*ab) - scales[j] * R0, 12)
            blocked = height(centers[i] + u[:, None] * e).max() > heights[i]
            for k in range(n):
                if k in (i, j):
                    continue
                d = (centers[k] - centers[i] + L / 2) % L - L / 2
                u = np.clip(d @ e, 0, np.hypot(*ab))
                if np.hypot(*(d - u * e)) < scales[k] * R0:
                    blocked = True
                    break
            if not blocked and (best is None or gap < best_gap):   # nearest way down
                best, best_gap = j, gap
        if best is not None:
            streams.append((i, best))
    return L, centers, scales, depths, shapes, kinds, heights, streams


def inside(pond, i, pts, t=0, margin=0.0):
    """Which points are in pond i's water, at least `margin` from its shore."""
    near, _, _, r, R, g, _ = pond.geom(pts, t)
    return (near == i) & (r < R * g - margin)


def area(pond, i, t=0):
    o = pond.outline(i, t, 360)
    x, y = o[:, 0], o[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


def free_points(world, pond, i, count, margin=2.0, gap=1.1, tries=2000):
    """Up to count random free points in pond i's water, clear of all beads and of each other."""
    pts, L = [], world.L
    tree = cKDTree(world._wrapped(), boxsize=L) if world.n else None
    for _ in range(count):
        for _ in range(tries):
            p, _ = pond.random_point(world.rng, world.t, margin, pond=i)
            if tree is not None and tree.query(p)[0] <= gap:
                continue
            if pts:
                d = np.array(pts) - p
                d -= L * np.round(d / L)
                if np.hypot(d[:, 0], d[:, 1]).min() <= gap:
                    continue
            pts.append(p)
            break
    return pts


def add_features(world, pond, kinds, rng):
    """Clay flakes, rock islands and hot springs. Returns hot spots [(x, y, radius, heat)]."""
    L, hot = world.L, []
    for i, kind in enumerate(kinds):
        c, R = pond.centers[i], pond.radii(0)[i]
        if kind == "clay":
            for _ in range(4):   # flakes of 7 charged mineral sites, frozen in the mud
                for _ in range(200):
                    p, _ = pond.random_point(rng, 0, 0.35 * R, pond=i)
                    a = rng.uniform(0, 2 * np.pi)
                    flake = (p + np.outer(np.arange(7) - 3, [np.cos(a), np.sin(a)])) % L
                    if world.n == 0 or cKDTree(world._wrapped(), boxsize=L).query(flake)[0].min() > 3:
                        break
                world.add_molecule(list(flake), ["mineral"] * 7)
                world.frozen[-7:] = True
        elif kind == "rocky":
            for _ in range(3):   # small round islands of rock
                p, _ = pond.random_point(rng, 0, 0.3 * R, pond=i)
                g = np.arange(-3, 3.01, 0.9)
                blob = np.array([[x, y] for x in g for y in g if x * x + y * y <= 2.8 ** 2]) + p
                tree = cKDTree(world._wrapped(), boxsize=L) if world.n else None
                if tree is not None:
                    blob = blob[tree.query(blob % L)[0] > 1.0]
                world.add_molecule(list(blob % L), ["rock"] * len(blob))
                world.frozen[-len(blob):] = True
        elif kind == "hot":
            p, _ = pond.random_point(rng, 0, 0.4 * R, pond=i)
            hot.append((p[0], p[1], 0.25 * R, 0.15))
    return hot


def build(args):
    rules = Rules()
    rules.bend_k = 50
    rules.bend_r0 = 2 * rules.bond_r0
    L, centers, scales, depths, shapes, kinds, heights, streams = landscape(args)
    world = World(box=L, rules=rules, kT=args.kT, seed=args.seed)
    world.repel_inside = (CARB,)
    timing = np.random.default_rng(args.land_seed + 2000)
    pond = Pond(L, args.r_wet, args.r_dry, args.period, rain=args.rain, seed=args.seed + 1000, centers=centers,
                phases=timing.uniform(-0.1, 0.1, args.ponds), periods=args.period * timing.uniform(0.85, 1.15, args.ponds),
                scales=scales, depths=depths, shapes=shapes, streams=None if args.no_streams else streams,
                stream_speed=args.stream_speed, shore_max=2.0,
                weather=None if args.no_weather else Weather(args.year, args.wind, args.storm_every, args.storm_length,
                                                             args.storm_wind, seed=args.seed + 3000,
                                                             current=args.current, surface=args.surface))
    pond.kinds, pond.heights = kinds, heights
    hot = add_features(world, pond, kinds, np.random.default_rng(args.land_seed + 4000))
    # lipids at random free spots inside each full pond (all 4 beads clear of everything)
    length = 3 * rules.bond_r0
    step = np.arange(4) * rules.bond_r0
    pond.lipids_wanted = []
    for i in range(args.ponds):
        want = int(args.density * area(pond, i))
        pond.lipids_wanted.append(want)
        placed, tree, fresh = 0, None, []
        for _ in range(200 * want):
            if placed == want:
                break
            if tree is None or len(fresh) >= 50:   # rebuild the search tree now and then
                tree, fresh = cKDTree(world._wrapped(), boxsize=L) if world.n else None, []
            p, _ = pond.random_point(world.rng, 0, length + 1.0, pond=i)
            a = world.rng.uniform(0, 2 * np.pi)
            beads = (p + np.outer(step, [np.cos(a), np.sin(a)])) % L
            if tree is not None and tree.query(beads)[0].min() < 1.1:
                continue
            if fresh:
                d = beads[:, None] - np.array(fresh)[None]
                d -= L * np.round(d / L)
                if np.hypot(d[..., 0], d[..., 1]).min() < 1.1:
                    continue
            add_lipid(world, *p, a, n_tails=3)
            fresh.extend(beads)
            placed += 1
        if placed < want:
            raise ValueError(f"pond {i}: room for only {placed} of {want} lipids")
    # ions: Na+ for the clay's charge, Ca2+ (+2 Cl-) for formose
    for i, kind in enumerate(kinds):
        if kind == "clay":
            for p in free_points(world, pond, i, 28):
                world.add_molecule([p], ["Na"])
        if args.sugar:
            for p in free_points(world, pond, i, args.sugar):
                world.add_molecule([p], ["sugar"])
        if args.chem:
            for name, count in (("Ca", args.ca), ("Cl", 2 * args.ca)):
                for p in free_points(world, pond, i, count):
                    world.add_molecule([p], [name])
    world.step(300, max_force=50.0)
    world.t = 0
    if pond.weather is not None:
        pond.weather.base_kT = args.kT
    if hot:
        spots = np.array(hot)

        def kT_field(pos):
            d = pos[:, None, :] - spots[None, :, :2]
            d -= L * np.round(d / L)
            r2 = np.einsum("ijk,ijk->ij", d, d)
            return world.kT + (spots[:, 3] * np.exp(-r2 / (2 * spots[:, 2] ** 2))).sum(1)

        world.kT_field = kT_field
    pond.hot = hot
    pond.home_beads = np.flatnonzero(world.type == INDEX["head"])
    pond.home = pond.where(world._wrapped()[pond.home_beads], 0)[0]
    return world, pond


class Chemistry(e05.Chemistry):
    def spots(self, world, n, pond, i, center, R):
        return free_points(world, pond, i, n, margin=1.0)


def measure(world, pond, args, chem, crust):
    s = e05.info(world, pond, 0.7, 0.7, chem, crust)
    heads = world._wrapped()[pond.home_beads]
    near, gap = pond.where(heads, world.t)[:2]
    wet = gap == 0
    s["pond_lipids"] = " ".join(str(int(((near == i) & wet).sum())) for i in range(len(pond.centers)))
    ves = np.zeros(len(pond.centers), dtype=int)
    if s["vesicle_xy"]:
        for k in pond.where(np.array(s["vesicle_xy"]).reshape(-1, 2), world.t)[0]:
            ves[k] += 1
    s["pond_vesicles"] = " ".join(map(str, ves))
    return s


def draw(view, world, pond, s):
    import matplotlib.patches as mpatches
    from alife.viz import WATER
    if getattr(view, "shores", None) is None:
        view.ax.set_facecolor(LAND)   # dry land; water is painted on per pond
        view.shores = []
        for i, kind in enumerate(pond.kinds):
            poly = mpatches.Polygon(pond.outline(i, 0), closed=True, fc=WATER, ec=KIND_COLORS[kind],
                                    lw=1.5, ls="--", zorder=0.5)
            view.ax.add_patch(poly)
            view.shores.append(poly)
            c = pond.centers[i]
            view.ax.text(c[0], c[1], kind, color=KIND_COLORS[kind], fontsize=14, alpha=0.35,
                         ha="center", va="center", zorder=1, weight="bold")
        for x, y, r, _ in pond.hot:
            view.ax.add_patch(mpatches.Circle((x, y), 1.5 * r, color="#F0997B", alpha=0.3, zorder=0.6, lw=0))
        view.streams = []
        for k in range(len(pond.streams)):
            corners = pond.stream_outline(k)
            bed = mpatches.Polygon(corners, closed=True, lw=0, zorder=0.45)
            view.ax.add_patch(bed)
            mid, e = corners.mean(0), (corners[1] - corners[0]) / 2
            sign = view.ax.annotate("", mid + 0.12 * e, mid - 0.12 * e, zorder=0.7,
                                    arrowprops=dict(arrowstyle="-|>", color="#185FA5", lw=1.5))
            view.streams.append((bed, sign))
        view.caption = view.ax.text(0.01, 0.99, "", transform=view.ax.transAxes, fontsize=9, va="top", zorder=5)
        view.wind = None
    for i, poly in enumerate(view.shores):
        poly.set_xy(pond.outline(i, world.t))
    for (bed, sign), f in zip(view.streams, pond.stream_flows(world.t)):
        bed.set_facecolor(WATER if f > 0 else "#D8C9A8")   # dry bed: darker sand
        sign.set_visible(f > 0)
    text = ""
    if s.get("season"):
        text = e05.e05_weather_text(s) + "\n"
        if view.wind is not None:
            view.wind.remove()
        a = np.radians(s["heading"])
        tip = np.array([0.93, 0.93]) + 0.05 * np.array([np.cos(a), np.sin(a)])
        view.wind = view.ax.annotate("wind", tip, 2 * np.array([0.93, 0.93]) - tip, xycoords="axes fraction",
                                     textcoords="axes fraction", ha="center", va="center", fontsize=8, zorder=6,
                                     arrowprops=dict(arrowstyle="-|>", color="k", lw=1 + s["wind"]))
    flows = pond.stream_flows(world.t)
    if flows:
        text += f"streams running {sum(f > 0 for f in flows)}/{len(flows)}  "
    text += f"ponds joined {s['links']}  lipids moved pond {s['travellers']}  vesicles per pond {s['pond_vesicles']}"
    if s.get("crust"):
        text += f"  crust {s['crust']}"
    if s.get("c1") or s.get("c2_4") or s.get("c5plus"):
        text += f"\nsugar carbons: C1 {s['c1']}  C2-4 {s['c2_4']}  C5+ {s['c5plus']}"
    view.caption.set_text(text)
    return view.update(s)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--layout", choices=("grid", "basin"), default="grid",
                   help="grid: ponds on a square grid; basin: scattered over rolling ground, streams run downhill")
    p.add_argument("--lake", type=float, default=1.5, help="basin: size of the lowest pond (x --r-wet)")
    p.add_argument("--reach", type=float, default=30.0, help="basin: longest stream across dry land")
    p.add_argument("--ponds", type=int, default=4)
    p.add_argument("--kinds", default=None, help="comma list per pond (default plain,clay,rocky,hot,...)")
    p.add_argument("--spacing", type=float, default=92.0, help="distance between pond centres")
    p.add_argument("--r-wet", type=float, default=42.0, help="typical full pond radius")
    p.add_argument("--r-dry", type=float, default=20.0, help="typical driest radius")
    p.add_argument("--period", type=int, default=20000, help="steps per dry + refill cycle")
    p.add_argument("--rain", type=float, default=1.0, help="rain splash speed while refilling")
    p.add_argument("--density", type=float, default=0.05, help="lipids per unit area of full pond")
    p.add_argument("--sugar", type=int, default=0, help="sugar beads per pond")
    p.add_argument("--chem", action="store_true")
    p.add_argument("--rain-carb", type=int, default=60, help="formaldehyde per refill per pond")
    p.add_argument("--ca", type=int, default=20, help="Ca2+ per pond (--chem)")
    p.add_argument("--meteor-every", type=int, default=50000)
    p.add_argument("--meteor-sugar", type=int, default=5)
    p.add_argument("--leak", action="store_true")
    p.add_argument("--crust", action="store_true")
    p.add_argument("--no-weather", action="store_true")
    p.add_argument("--no-streams", action="store_true", help="ponds never overflow into each other")
    p.add_argument("--stream-speed", type=float, default=1.0, help="top water speed in a full stream")
    p.add_argument("--year", type=int, default=80000)
    p.add_argument("--wind", type=float, default=1.5)
    p.add_argument("--storm-every", type=int, default=40000)
    p.add_argument("--storm-length", type=int, default=4000)
    p.add_argument("--storm-wind", type=float, default=2.0)
    p.add_argument("--current", type=float, default=0.3,
                   help="how hard wind and storms stir the whole pond (1 = as e05)")
    p.add_argument("--surface", type=float, default=0.3,
                   help="wind drag on floating lipids (piles them on the downwind shore; 0 = off)")
    p.add_argument("--max-force", type=float, default=300.0,
                   help="cap on any one push (safety net against beads flung into each other)")
    p.add_argument("--kT", type=float, default=0.6)
    p.add_argument("--land-seed", type=int, default=1, help="which landscape")
    p.add_argument("--seed", type=int, default=0, help="which history")
    p.add_argument("--steps", type=int, default=160000)
    p.add_argument("--log-every", type=int, default=400)
    p.add_argument("--video", default=None)
    p.add_argument("--fps", type=int, default=10)
    p.add_argument("--out", default="runs/e09_landscape")
    args = p.parse_args()

    world, pond = build(args)
    chem = Chemistry(args) if args.chem else None
    crust = Crust() if args.crust else None
    n_lip = len(pond.home_beads)
    print(f"landscape {args.land_seed}: {args.ponds} ponds {pond.kinds}, box {world.L:.0f}, "
          f"{n_lip} lipids {pond.lipids_wanted}, {world.n} beads", flush=True)
    out = os.path.join(args.out, f"{'land' if args.layout == 'grid' else 'basin'}{args.land_seed}_seed{args.seed}")
    os.makedirs(out, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    from alife.viz import Movie, View
    view = View(world, title=f"landscape {args.land_seed}", size=9)
    movie = None
    if args.video:
        view.legend()
        movie = Movie(view, args.video, fps=args.fps)
    log = open(os.path.join(out, "log.csv"), "w", newline="")
    writer = csv.DictWriter(log, FIELDS, extrasaction="ignore")
    writer.writeheader()
    start = time.time()
    while world.t < args.steps:
        e05.advance(world, pond, args.log_every, chem=chem, crust=crust, max_force=args.max_force)
        s = measure(world, pond, args, chem, crust)
        writer.writerow(s)
        log.flush()
        if movie:
            draw(view, world, pond, s)
            movie.grab()
        if world.t % (args.log_every * 25) == 0:
            draw(view, world, pond, s)
            view.save(os.path.join(out, f"step_{world.t:06d}.png"), s)
            print(f"step {world.t} {s.get('season', '')} lipids/pond {s['pond_lipids']} vesicles/pond "
                  f"{s['pond_vesicles']} links {s['links']} travellers {s['travellers']} C5+ {s['c5plus']} "
                  f"({time.time() - start:.0f}s)", flush=True)
    log.close()
    if movie:
        movie.close()
        print(f"video: {args.video} ({movie.frames} frames)")


if __name__ == "__main__":
    main()
