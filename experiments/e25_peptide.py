"""e25: do peptides make pores in a protocell wall by themselves? (gated channels step 2)

Same 3D fatty-acid bag as e24, but no osmotic push (same sugar inside and out). Add short stiff
amphipathic peptides (alife.chem.molecules.peptide3: a water-loving stripe and an oily stripe, about as
long as the wall is thick; like alamethicin / magainin). Nothing says "make a pore": if the peptides
gather with their water-loving stripes facing each other, a water channel through the wall can appear.

--start wall: peptides start one by one, spanning the wall at random spots (shortcut: skips getting in).
--start water: they start in the water outside and must get in by themselves.
--peptides 0 = control bag. Probes (small "carb") and sugar start inside: do they get out faster?

Logs as e24 (stretch, holes, probes/sugar still inside, pockets) plus peptides spanning the wall /
lying in it / free, and peptide clusters (water-loving stripes touching).

    C:\\Apps\\Python314\\python.exe experiments/e25_peptide.py --peptides 40 --start wall --seed 0
"""
import argparse, csv, os, sys, time
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))
from alife.chem.elements import INDEX   # noqa: E402
from alife.chem.molecules import add_lipid3, peptide3, peptide_points   # noqa: E402
from alife.chem.world import World   # noqa: E402
from e15_split3d import biggest, fa_head, offset, scatter_fast   # noqa: E402
from e24_stretch import FIELDS as F24, measure   # noqa: E402

PEP = INDEX["pep"]
FIELDS = F24 + ["span", "flat", "free", "pep_cluster", "clusters3"]


def sphere_dirs(n):
    k = np.arange(n) + 0.5
    z = 1 - 2 * k / n
    phi = np.pi * (1 + 5 ** 0.5) * k
    return np.stack([np.sqrt(1 - z ** 2) * np.cos(phi), np.sqrt(1 - z ** 2) * np.sin(phi), z], 1)


def random_dirs(rng, n, min_angle):
    out = []
    while len(out) < n:
        d = rng.normal(size=3)
        d /= np.linalg.norm(d)
        if all(d @ o < np.cos(min_angle) for o in out):
            out.append(d)
    return np.array(out)


def perp(rng, u):
    s = np.cross(u, rng.normal(size=3))
    return s / np.linalg.norm(s)


def bag(w, mid, radius, mix, sites, clear=1.0):
    """e15 mixed_vesicle, but leaves room (no lipid head within `clear` of the line) at each site direction."""
    r0 = w.rules.bond_r0
    for side in (1, -1):
        r_head = radius + side * (0.5 + 2 * r0)
        for d in sphere_dirs(int(4 * np.pi * r_head ** 2 / 1.2)):
            p = r_head * d
            if len(sites):
                along = sites @ p
                if np.any((along > 0) & (np.linalg.norm(p - along[:, None] * sites, axis=1) < clear)):
                    continue
            if w.rng.uniform() < mix:
                add_lipid3(w, mid + p, -side * d, 2, head=fa_head(w, 0.0), tail="stail")
            else:
                add_lipid3(w, mid + p, -side * d, 2)


def peptide_stats(w, peps, c, R):
    """Spanning (in the wall, upright), flat (in the wall, lying), free; clusters by touching water-loving stripes."""
    x = offset(w, w._wrapped(), c)
    span = flat = free = 0
    for a, b in peps:   # a, b = first and last bead of stripe A
        m = (x[a] + x[b]) / 2
        ax = x[b] - x[a]
        ax /= np.linalg.norm(ax)
        r = np.linalg.norm(m)
        if abs(r - R) < 3.0:
            if abs(ax @ m / r) > 0.7:
                span += 1
            else:
                flat += 1
        else:
            free += 1
    owner = np.full(w.n, -1)
    for k, (a, b) in enumerate(peps):
        owner[np.flatnonzero(w.mol == w.mol[a])] = k
    polar = np.flatnonzero((owner >= 0) & (w.type == PEP))
    pr = cKDTree(w._wrapped()[polar], boxsize=w.L).query_pairs(1.3, output_type="ndarray")
    i, j = owner[polar[pr[:, 0]]], owner[polar[pr[:, 1]]]
    keep = i != j
    g = connected_components(coo_matrix((np.ones(keep.sum()), (i[keep], j[keep])), shape=(len(peps),) * 2),
                             directed=False)[1]
    size = np.bincount(g)
    return {"span": span, "flat": flat, "free": free, "pep_cluster": int(size.max()), "clusters3": int((size >= 3).sum())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--radius", type=float, default=10.0)
    ap.add_argument("--box", type=float, default=40.0)
    ap.add_argument("--mix", type=float, default=0.5)
    ap.add_argument("--osm", type=float, default=0.03, help="sugar per unit volume, inside and outside")
    ap.add_argument("--probes", type=int, default=40)
    ap.add_argument("--peptides", type=int, default=40)
    ap.add_argument("--start", choices=["wall", "water"], default="wall")
    ap.add_argument("--length", type=int, default=7, help="rungs per peptide (wall is ~6 beads thick)")
    ap.add_argument("--kT", type=float, default=1.1)
    ap.add_argument("--steps", type=int, default=60000)
    ap.add_argument("--every", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="runs/e25")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    w = World(box=args.box, kT=args.kT, seed=args.seed, dim=3)
    rng = np.random.default_rng(args.seed + 1000)
    mid = np.full(3, args.box / 2)
    R, L = args.radius, args.length
    wall_start = args.start == "wall" and args.peptides > 0
    sites = random_dirs(rng, args.peptides, 4.0 / R) if wall_start else np.zeros((0, 3))
    bag(w, mid, R, args.mix, sites)
    peps = []
    for u in sites:
        peptide3(w, mid + R * u, u, perp(rng, u), L)
        peps.append((w.n - 2 * L, w.n - L - 1))
    inner, outer = R - 2.5, R + 2.5
    r_in = inner - 1.0
    dist = lambda p: np.linalg.norm(offset(w, p, mid), axis=1)
    if args.start == "water":
        while len(peps) < args.peptides:
            c = rng.uniform(0, args.box, 3)
            u = rng.normal(size=3)
            u /= np.linalg.norm(u)
            s = perp(rng, u)
            pts = peptide_points(c, u, s, L, w.rules.bond_r0) % args.box
            if dist(pts).min() < outer + 2.0 or cKDTree(w._wrapped(), boxsize=w.L).query(pts)[0].min() < 1.1:
                continue
            peptide3(w, c, u, s, L)
            peps.append((w.n - 2 * L, w.n - L - 1))
    n0 = w.n
    scatter_fast(w, int(round(args.osm * 4 / 3 * np.pi * r_in ** 3)), lambda p: dist(p) < r_in)
    sugar_in = np.arange(n0, w.n)
    scatter_fast(w, int(round(args.osm * (args.box ** 3 - 4 / 3 * np.pi * (outer + 1.0) ** 3))),
                 lambda p: dist(p) > outer + 1.0)
    n1 = w.n
    scatter_fast(w, args.probes, lambda p: dist(p) < r_in, name="carb")
    probes = np.arange(n1, w.n)
    w.step(300, max_force=50)
    w.t = 0
    name = f"pep_R{R:g}_mix{args.mix:g}_osm{args.osm:g}_n{args.peptides}_{args.start}_L{L}_p{args.probes}_s{args.seed}"
    print(name, "beads", w.n, "sugar in", len(sugar_in), flush=True)
    t0 = time.time()

    def row():
        r = measure(w, sugar_in, probes)
        if peps:
            _, c, _ = biggest(w)
            r.update(peptide_stats(w, peps, c, r["R"]))
        return r

    with open(f"{args.out}/{name}.csv", "w", newline="") as fh:
        log = csv.DictWriter(fh, FIELDS)
        log.writeheader()
        log.writerow(row())
        while w.t < args.steps:
            w.step(args.every)
            r = row()
            log.writerow(r)
            fh.flush()
            print(" ".join(f"{k} {v}" for k, v in r.items()), f"({time.time() - t0:.0f}s)", flush=True)
    w.save(f"{args.out}/{name}_state.pkl", sugar_in=sugar_in, probes=probes, peps=peps)


if __name__ == "__main__":
    main()
