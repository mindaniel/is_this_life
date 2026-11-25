"""e26: water as beads. Does a wall hold in real (bead) water, and do peptides open a water channel in it?

Up to e25 water was invisible (implicit): oily beads simply attract each other. That made two things too cheap:
a water-loving bead buried in oil loses nothing, and a channel through the wall holds nothing. Here water is
beads ("wat"). Water sticks to water and to water-loving beads (heads, peptide polar stripe) with depth --ew;
it doesn't stick to oily tails. Tails keep their own attraction (van der Waals). Nothing says "make a pore".

A flat wall (half fatty acid, like the e24/e25 bag) spans the box in x-y; water fills above and below
(they join through the periodic z wall, so the only way across the wall itself is through it).
Logged every --every steps:
  cross   water beads that went fully through the wall (from one side's water to the other's)
  in_core water beads in the oily middle (|z - mid| < 1)
  thick   wall thickness (head layer to head layer)
  rough   how much the wall bends (spread of its local middle height)
Sides are measured from the wall's local middle (nearest tail beads in x-y), so riding over a bend is no crossing.
  D_wat   water diffusion (is it liquid, not frozen?)
  span / pep_cluster / clusters3 for peptides (as e25).

--bundle k (shortcut, flagged): peptides start as ready-made rings of k (water-loving stripes facing in, water in
the middle) = a pore put in by hand. Asks the other question: if a pore existed, would it stay and carry water?

  C:\\Apps\\Python314\\python.exe experiments/e26_water.py --ew 0.5 --peptides 0 --steps 10000
"""
import argparse, csv, os, sys, time
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from alife.chem.elements import INDEX   # noqa: E402
from alife.chem.molecules import add_lipid3, peptide3   # noqa: E402
from alife.chem.world import TAILS, World   # noqa: E402

WAT, PEP = INDEX["wat"], INDEX["pep"]
POLAR = ["head", "fa_head", "fa_acid", "fa_oh", "pep", "sugar", "carb", "Na", "K", "Cl", "Ca", "Mg"]
HEADS = [INDEX[h] for h in ("head", "fa_head", "fa_acid", "fa_oh")]
FIELDS = ["step", "cross", "in_core", "thick", "rough", "D_wat", "span", "pep_cluster", "clusters3"]


def wall(w, z0, nx, ny, mix, holes, clear=1.0):
    """Flat two-leaflet wall at height z0 on an nx x ny lattice (rows shifted: hexagon-like), skipping
    lipids whose head is within `clear` (in x-y) of a peptide spot."""
    Lx, Ly = w.box()[:2]
    for a in range(ny):
        for b in range(nx):
            xy = np.array([(b + 0.5 * (a % 2)) * Lx / nx, (a + 0.5) * Ly / ny])
            if len(holes):
                d = (holes - xy + np.array([Lx, Ly]) / 2) % np.array([Lx, Ly]) - np.array([Lx, Ly]) / 2
                if np.min(np.hypot(*d.T)) < clear:
                    continue
            for side in (1, -1):
                top = z0 + side * 2.5
                if w.rng.uniform() < mix:
                    add_lipid3(w, [xy[0], xy[1], top], [0, 0, -side], 2, head="fa_head", tail="stail")
                else:
                    add_lipid3(w, [xy[0], xy[1], top], [0, 0, -side], 2)


def fill_water(w, rho, z0, gap=0.85):
    """Water beads on a jittered cubic lattice at density rho, outside the wall and away from all beads."""
    box = w.box()
    a = rho ** (-1 / 3)
    n = np.maximum(1, np.round(box / a)).astype(int)
    g = np.stack(np.meshgrid(*[(np.arange(k) + 0.5) * box[i] / k for i, k in enumerate(n)], indexing="ij"), -1)
    g = g.reshape(-1, 3) + w.rng.uniform(-0.05, 0.05, (g.size // 3, 3))
    dz = (g[:, 2] - z0 + box[2] / 2) % box[2] - box[2] / 2
    g = g[np.abs(dz) > 2.9]
    g = g[cKDTree(w._wrapped(), boxsize=box).query(g % box)[0] > gap]
    for p in g:
        w.add_molecule([p], ["wat"])
    return len(g)


def fill_channel(w, centers, r, z0, half=2.9, gap=0.85):
    """Water beads inside each ready-made pore (within r of its axis, inside the wall)."""
    box = w.box()
    n = 0
    for c in centers:
        tree = cKDTree(w._wrapped(), boxsize=box)
        for x in np.arange(-r, r + 1e-9, 0.9):
            for y in np.arange(-r, r + 1e-9, 0.9):
                if x * x + y * y > r * r:
                    continue
                for z in np.arange(z0 - half, z0 + half + 1e-9, 0.9):
                    p = np.array([c[0] + x, c[1] + y, z]) % box
                    if tree.query(p)[0] > gap:
                        w.add_molecule([p], ["wat"])
                        tree = cKDTree(w._wrapped(), boxsize=box)
                        n += 1
    return n


def mid_z(w):
    z = w._wrapped()[np.isin(w.type, TAILS), 2]
    Lz = w.box()[2]
    ang = 2 * np.pi * z / Lz
    return (np.arctan2(np.sin(ang).mean(), np.cos(ang).mean()) * Lz / (2 * np.pi)) % Lz


def local_offset(w, idx, zm, k=16):
    """Height of beads idx above the wall's LOCAL middle (mean height of the k nearest tail beads in x-y), and the
    wall's roughness (spread of those local middles). Measuring from the average height only would count a water
    bead that rides over a bend of the wall as a crossing."""
    box = w.box()
    p = w._wrapped()
    tails = np.flatnonzero(np.isin(w.type, TAILS))
    tz = (p[tails, 2] - zm + box[2] / 2) % box[2] - box[2] / 2
    _, nb = cKDTree(p[tails, :2], boxsize=box[:2]).query(p[idx, :2], k=k)
    local = tz[nb].mean(1)
    dz = (p[idx, 2] - zm + box[2] / 2) % box[2] - box[2] / 2
    return dz - local, local


def peptide_stats(w, peps):
    span = 0
    for a, b in peps:
        ax = w.delta(np.array([a]), np.array([b]))[0]
        span += abs(ax[2]) / np.linalg.norm(ax) > 0.7
    owner = np.full(w.n, -1)
    for k, (a, b) in enumerate(peps):
        owner[np.flatnonzero(w.mol == w.mol[a])] = k
    polar = np.flatnonzero((owner >= 0) & (w.type == PEP))
    pr = cKDTree(w._wrapped()[polar], boxsize=w.box()).query_pairs(1.3, output_type="ndarray")
    i, j = owner[polar[pr[:, 0]]], owner[polar[pr[:, 1]]]
    keep = i != j
    g = connected_components(coo_matrix((np.ones(keep.sum()), (i[keep], j[keep])), shape=(len(peps),) * 2),
                             directed=False)[1]
    size = np.bincount(g)
    return {"span": int(span), "pep_cluster": int(size.max()), "clusters3": int((size >= 3).sum())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=float, default=20.0, help="box side in x and y")
    ap.add_argument("--Lz", type=float, default=20.0)
    ap.add_argument("--nx", type=int, default=18)
    ap.add_argument("--ny", type=int, default=20)
    ap.add_argument("--mix", type=float, default=0.5)
    ap.add_argument("--rho", type=float, default=0.75, help="water beads per unit volume")
    ap.add_argument("--ew", type=float, default=0.5, help="water-water and water-polar stickiness (tail-tail = 1)")
    ap.add_argument("--peptides", type=int, default=0)
    ap.add_argument("--length", type=int, default=7)
    ap.add_argument("--gap", type=float, default=3.0, help="min distance between scattered peptide spots")
    ap.add_argument("--bundle", type=int, default=0, help="start peptides as ready-made rings of this many (0 = alone)")
    ap.add_argument("--ring-r", type=float, default=2.5, help="ring radius (peptide centres) for --bundle")
    ap.add_argument("--kT", type=float, default=1.1)
    ap.add_argument("--steps", type=int, default=10000)
    ap.add_argument("--every", type=int, default=500)
    ap.add_argument("--sample", type=int, default=100, help="steps between crossing checks")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="runs/e26")
    ap.add_argument("--start", default="", help="carry on from a saved _state.pkl (--steps more)")
    ap.add_argument("--tag", default="", help="added to the output name")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    if args.start:
        w, ex = World.load(args.start)
        peps, water = [tuple(x) for x in ex["peps"]], ex["water"]
        lipids, n_wat, n_chan, k = len(np.unique(w.mol[np.isin(w.type, HEADS)])), len(water), 0, args.bundle or 1
        t_end = w.t + args.steps
    else:
        w = World(box=(args.L, args.L, args.Lz), kT=args.kT, seed=args.seed, dim=3)
        for p in POLAR:
            w.rules.set_pair_eps("wat", p, args.ew)
        w.rules.set_pair_eps("wat", "wat", args.ew)
        rng = np.random.default_rng(args.seed + 1000)
        z0 = args.Lz / 2
        dist2 = lambda a, b: np.hypot(*((a - b + args.L / 2) % args.L - args.L / 2))
        k = args.bundle or 1
        if args.peptides % k:
            raise SystemExit("--peptides must be a multiple of --bundle")
        min_gap = args.gap if k == 1 else 2 * args.ring_r + 4.0
        holes = []
        for _ in range(200000):   # spot per peptide (alone) or per ring (bundle), apart from each other
            if len(holes) == args.peptides // k:
                break
            h = rng.uniform(0, args.L, 2)
            if all(dist2(h, o) > min_gap for o in holes):
                holes.append(h)
        else:
            raise SystemExit("no room for the peptide spots: lower --gap")
        holes = np.array(holes).reshape(-1, 2)
        if k == 1:
            wall(w, z0, args.nx, args.ny, args.mix, holes)
        else:
            wall(w, z0, args.nx, args.ny, args.mix, holes, clear=args.ring_r + 1.0)
        lipids = w.n_mol
        peps = []
        for h in holes:
            if k == 1:
                s = rng.normal(size=2)
                peptide3(w, [h[0], h[1], z0], [0, 0, 1], [s[0], s[1], 0], args.length)
                peps.append((w.n - 2 * args.length, w.n - args.length - 1))
                continue
            a0 = rng.uniform(0, 2 * np.pi)
            for j in range(k):   # side vector points out: water-loving stripe A sits on the inner side
                u = np.array([np.cos(a0 + 2 * np.pi * j / k), np.sin(a0 + 2 * np.pi * j / k), 0.0])
                peptide3(w, [h[0] + args.ring_r * u[0], h[1] + args.ring_r * u[1], z0], [0, 0, 1], u, args.length)
                peps.append((w.n - 2 * args.length, w.n - args.length - 1))
        n_chan = fill_channel(w, holes, args.ring_r - 1.4, z0) if k > 1 else 0
        n_wat = fill_water(w, args.rho, z0) + n_chan
        water = np.flatnonzero(w.type == WAT)
        w.step(300, max_force=50)
        w.t = 0
        t_end = args.steps
    name = (f"wat_L{args.L:g}x{args.Lz:g}_mix{args.mix:g}_rho{args.rho:g}_ew{args.ew:g}_n{args.peptides}"
            f"_L{args.length}" + (f"_b{k}" if k > 1 else "") + f"_kT{args.kT:g}_s{args.seed}{args.tag}")
    if args.start:
        name = os.path.basename(args.start).replace("_state.pkl", "") + (args.tag or "_more")
    print(name, "lipids", lipids, "water", n_wat, "in channels", n_chan, "beads", w.n, flush=True)
    Lz = args.Lz
    zm = mid_z(w)
    # side counter: floor((z_unwrapped - mid)/Lz) changes only when a bead passes the wall's middle plane
    zu = w.pos[water, 2].copy()   # pos is unwrapped (World keeps it continuous)
    h, local = local_offset(w, water, zm)
    stable = np.where(np.abs(h) > 3.5, np.floor((zu - zm - local) / Lz), np.nan)   # nan = not yet settled on a side
    cross = 0
    t0 = time.time()
    with open(f"{args.out}/{name}.csv", "w", newline="") as fh:
        log = csv.DictWriter(fh, FIELDS)
        log.writeheader()
        p_ref, t_ref = w.pos[water].copy(), w.t
        while w.t < t_end:
            w.step(args.sample)
            zm = zm + ((mid_z(w) - zm + Lz / 2) % Lz - Lz / 2)   # follow the wall's drift without wrapping
            z = w.pos[water, 2]
            h, local = local_offset(w, water, zm)
            side = np.floor((z - zm - local) / Lz)
            off = np.abs(h)   # distance from the wall's local middle
            settled = off > 3.5   # back in open water
            moved = settled & ~np.isnan(stable) & (side != stable)
            cross += int(moved.sum())
            stable[settled] = side[settled]
            if w.t % args.every == 0:
                hz = (w._wrapped()[np.isin(w.type, HEADS), 2] - zm + Lz / 2) % Lz - Lz / 2
                thick = float(np.median(hz[hz > 0]) - np.median(hz[hz < 0]))
                dp = w.pos[water] - p_ref
                D = float(np.mean(np.sum(dp ** 2, 1)) / (6 * (w.t - t_ref) * w.dt))
                p_ref, t_ref = w.pos[water].copy(), w.t
                r = {"step": w.t, "cross": cross, "in_core": int((off < 1.0).sum()), "thick": round(thick, 2), "rough": round(float(np.std(local)), 2),
                     "D_wat": round(D, 4)}
                if peps:
                    r.update(peptide_stats(w, peps))
                log.writerow(r)
                fh.flush()
                print(" ".join(f"{k} {v}" for k, v in r.items()), f"({time.time() - t0:.0f}s)", flush=True)
    w.save(f"{args.out}/{name}_state.pkl", peps=peps, water=water)


if __name__ == "__main__":
    main()
