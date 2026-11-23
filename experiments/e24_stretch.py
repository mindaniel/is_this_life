"""e24: does a stretched wall open by itself? (first step toward gated channels)

A 3D fatty-acid vesicle (e15 build) with more sugar inside than outside (--osm-in vs --osm): the trapped
sugar pushes the wall out (osmotic pressure, like a cell in fresh water). Small neutral probe molecules
("carb" beads, formaldehyde-sized; no reactions run) start inside. Nothing tells them how to cross: no
Leak rule, they can only get out by physically squeezing through the wall.

Question: does more stretch -> thinner wall, more holes, more probes out? If yes, a stretch-gated
"channel" already exists in plain physics (real fatty-acid vesicles do this: tension opens brief pores).

Logged every --every steps: wall lipids, wall radius, area per lipid, hole share (directions from the
centre with no tail bead in the wall shell), probes and sugar still inside, pockets.

    C:\\Apps\\Python314\\python.exe experiments/e24_stretch.py --osm-in 0.2 --seed 0
"""
import argparse, csv, os, sys, time
import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from alife.chem.analysis import water_pockets3   # noqa: E402
from alife.chem.world import TAILS, World   # noqa: E402
from e15_split3d import biggest, mixed_vesicle, offset, scatter_fast   # noqa: E402

FIELDS = ["step", "wall", "R", "area", "holes", "probes_in", "sugar_in", "pockets", "biggest_pocket"]


def measure(w, sugar, probes, n_dir=400):
    mask, c, _ = biggest(w)
    tails = mask & np.isin(w.type, TAILS)
    d = offset(w, w._wrapped()[tails], c)
    r = np.linalg.norm(d, axis=1)
    R = float(np.median(r))
    wall = int(len(np.unique(w.mol[mask])))
    # holes: directions where no wall tail bead sits within the shell, cone half-width ~1 bead at R
    k = np.arange(n_dir) + 0.5
    z = 1 - 2 * k / n_dir
    phi = np.pi * (1 + 5 ** 0.5) * k
    u = np.stack([np.sqrt(1 - z ** 2) * np.cos(phi), np.sqrt(1 - z ** 2) * np.sin(phi), z], 1)
    near = cKDTree(d).query(u * R, distance_upper_bound=1.5)[0]
    holes = float(np.mean(near == np.inf))
    inside = lambda b: int((np.linalg.norm(offset(w, w._wrapped()[b], c), axis=1) < R).sum())
    vol, _ = water_pockets3(w)
    return {"step": w.t, "wall": wall, "R": round(R, 2), "area": round(4 * np.pi * R ** 2 / (wall / 2), 3),
            "holes": round(holes, 3), "probes_in": inside(probes), "sugar_in": inside(sugar),
            "pockets": len(vol), "biggest_pocket": int(max(vol)) if len(vol) else 0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--radius", type=float, default=10.0)
    ap.add_argument("--box", type=float, default=40.0)
    ap.add_argument("--mix", type=float, default=0.5, help="fatty-acid (short-tail) share of the wall")
    ap.add_argument("--osm", type=float, default=0.03, help="sugar per unit volume outside")
    ap.add_argument("--osm-in", type=float, default=0.03, help="sugar per unit volume inside (more = stretched)")
    ap.add_argument("--probes", type=int, default=40, help="small neutral molecules that start inside")
    ap.add_argument("--kT", type=float, default=1.1)
    ap.add_argument("--steps", type=int, default=60000)
    ap.add_argument("--every", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="runs/e24")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    w = World(box=args.box, kT=args.kT, seed=args.seed, dim=3)
    mid = np.full(3, args.box / 2)
    mixed_vesicle(w, mid, args.radius, args.mix)
    inner, outer = args.radius - 2.5, args.radius + 2.5
    r_in = inner - 1.0
    n0 = w.n
    n_in = int(round(args.osm_in * 4 / 3 * np.pi * r_in ** 3))
    n_out = int(round(args.osm * (args.box ** 3 - 4 / 3 * np.pi * (outer + 1.0) ** 3)))
    scatter_fast(w, n_in, lambda p: np.linalg.norm(offset(w, p, mid), axis=1) < r_in)
    sugar_in = np.arange(n0, w.n)
    scatter_fast(w, n_out, lambda p: np.linalg.norm(offset(w, p, mid), axis=1) > outer + 1.0)
    n1 = w.n
    scatter_fast(w, args.probes, lambda p: np.linalg.norm(offset(w, p, mid), axis=1) < r_in, name="carb")
    probes = np.arange(n1, w.n)
    w.step(300, max_force=50)
    w.t = 0
    name = f"stretch_R{args.radius:g}_mix{args.mix:g}_out{args.osm:g}_in{args.osm_in:g}_p{args.probes}_s{args.seed}"
    print(name, "sugar in", n_in, "out", n_out, "beads", w.n, flush=True)
    t0 = time.time()
    with open(f"{args.out}/{name}.csv", "w", newline="") as fh:
        log = csv.DictWriter(fh, FIELDS)
        log.writeheader()
        row = measure(w, sugar_in, probes)
        log.writerow(row)
        while w.t < args.steps:
            w.step(args.every)
            row = measure(w, sugar_in, probes)
            log.writerow(row)
            fh.flush()
            print(" ".join(f"{k} {v}" for k, v in row.items()), f"({time.time() - t0:.0f}s)", flush=True)
    w.save(f"{args.out}/{name}_state.pkl", sugar_in=sugar_in, probes=probes)


if __name__ == "__main__":
    main()
