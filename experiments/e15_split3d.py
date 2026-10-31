"""e15: can a 3D vesicle grow and split? (fatty-acid feeding, Zhu & Szostak 2009)

In the lab, fatty-acid vesicles fed fresh fatty acid grow wall faster than
their volume: the solutes inside can't cross the membrane, so water can't
follow. With too much wall for a ball, the vesicle goes long (a thread), and
gentle shear or jiggling pinches the thread into small vesicles that keep
most of their contents. In 2D (e10) this can't happen: a 2D "tube" is two
parallel lines with no pinching drive.

Setup: one ready-made vesicle (--radius), a --mix fraction of its lipids
are short-tail ("stail", oiliness --short-oily): fatty-acid-like, they leave
and join membranes. Sugar beads (can't cross membranes) fill the water
inside and outside at the same concentration (--osm per unit volume): they
pin the inside volume, as in a real pond. Every --feed steps --batch loose
short-tail lipids are added at free spots (within --near of the vesicle's
surface, or anywhere if 0). Optional shear: water slides in bands
u_x = shear sin(2 pi z / L).

Logged every --every steps: lipid clusters (biggest, how many of 30+),
sealed pockets (vesicle insides) and their volumes, how long the biggest is
(shape: 1 = ball, larger = longer), loose lipids, lipids added.

    python experiments/e15_split3d.py --feed 200 --steps 200000 --video runs/video/e15.mp4 --fps 6
    python experiments/e15_split3d.py --feed 0     # control: no food
"""
import argparse
import csv
import os
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from alife.chem.analysis import cluster_labels, water_pockets3   # noqa: E402
from alife.chem.molecules import add_lipid3   # noqa: E402
from alife.chem.world import STAIL, TAILS, World   # noqa: E402

FIELDS = ["step", "biggest", "big_short", "clusters30", "loose", "pockets", "volumes", "length", "added"]


def mixed_vesicle(w, center, radius, mix, area=1.2):
    """vesicle3 with each lipid short-tail ("stail") with chance mix. Returns the count."""
    r0 = w.rules.bond_r0
    count = 0
    for side in (1, -1):   # outer leaflet: heads out; inner: heads in
        r_head = radius + side * (0.5 + 2 * r0)
        n = int(4 * np.pi * r_head ** 2 / area)
        k = np.arange(n) + 0.5
        z = 1 - 2 * k / n
        phi = np.pi * (1 + 5 ** 0.5) * k
        u = np.stack([np.sqrt(1 - z ** 2) * np.cos(phi), np.sqrt(1 - z ** 2) * np.sin(phi), z], 1)
        for d in u:
            add_lipid3(w, center + r_head * d, -side * d, 2, tail="stail" if w.rng.uniform() < mix else "tail")
        count += n
    return count


def scatter_fast(w, count, keep, gap=1.1, name="sugar"):
    """Same as scatter, batched: candidates far from all beads and from each other."""
    placed, empty = 0, 0
    while placed < count:
        tree = cKDTree(w._wrapped(), boxsize=w.L)
        p = w.rng.uniform(0, w.L, (200000, w.dim))
        p = p[keep(p)]
        p = p[tree.query(p)[0] > gap]
        if not len(p):
            empty += 1
            if empty > 20:
                raise RuntimeError("no room left")
            continue
        own = cKDTree(p, boxsize=w.L)
        take = []
        bad = np.zeros(len(p), bool)
        for i in range(len(p)):
            if bad[i]:
                continue
            take.append(i)
            bad[own.query_ball_point(p[i], gap)] = True
            if placed + len(take) == count:
                break
        for i in take:
            w.add_molecule([p[i]], [name])
        placed += len(take)
    return placed


def offset(w, pos, center):
    d = pos - center
    return d - w.L * np.round(d / w.L)


def centre(w, beads):
    """Mean position of beads in the periodic box (circular mean)."""
    ang = 2 * np.pi * w._wrapped()[beads] / w.L
    return np.angle(np.exp(1j * ang).mean(0)) % (2 * np.pi) / (2 * np.pi) * w.L


def biggest(w):
    """Bead mask, centre and 90th-percentile radius of the biggest lipid cluster."""
    lab = cluster_labels(w)
    tails = np.isin(w.type, TAILS)
    lipids = np.unique(w.mol[tails])
    big = np.bincount(lab[lipids]).argmax()
    mask = np.isin(w.mol, lipids[lab[lipids] == big])
    c = centre(w, mask)
    r = np.linalg.norm(offset(w, w._wrapped()[mask], c), axis=1)
    return mask, c, np.percentile(r, 90)


def feed(w, count, near, clear=1.5, tries=20000):
    """count loose short-tail lipids at free spots (within `near` outside the vesicle, or anywhere)."""
    r0 = w.rules.bond_r0
    _, c, R = biggest(w) if near else (None, None, None)
    added = 0
    for _ in range(count):
        tree = cKDTree(w._wrapped(), boxsize=w.L)
        for _ in range(tries):
            if near:
                d = w.rng.normal(size=3)
                d /= np.linalg.norm(d)
                p = (c + (R + 2 + near * w.rng.uniform()) * d) % w.L
            else:
                p = w.rng.uniform(0, w.L, 3)
            a = w.rng.normal(size=3)
            a /= np.linalg.norm(a)
            beads = (p + np.outer(np.arange(3) * r0, a)) % w.L
            if tree.query(beads)[0].min() > clear:
                add_lipid3(w, p, a, 2, tail="stail")
                added += 1
                break
    return added


def shear_flow(L, shear):
    k = 2 * np.pi / L

    def flow(pos):
        v = np.zeros_like(pos)
        v[:, 0] = shear * np.sin(k * pos[:, 2])
        return v
    return flow


def build(args):
    w = World(box=args.box, kT=args.kT, seed=args.seed, dim=3)
    w.rules.set_oily("stail", args.short_oily)
    mid = np.full(3, args.box / 2)
    n = mixed_vesicle(w, mid, args.radius, args.mix)
    inner, outer = args.radius - 2.5, args.radius + 2.5
    if args.osm:
        r_in = inner - 1.0
        n_in = int(round(args.osm * 4 / 3 * np.pi * r_in ** 3))
        n_out = int(round(args.osm * (args.box ** 3 - 4 / 3 * np.pi * (outer + 1.0) ** 3)))
        scatter_fast(w, n_in, lambda p: np.linalg.norm(offset(w, p, mid), axis=1) < r_in)
        scatter_fast(w, n_out, lambda p: np.linalg.norm(offset(w, p, mid), axis=1) > outer + 1.0)
        w.n_sugar_in = n_in
    w.step(300, max_force=50)
    w.t = 0
    w.n_start = w.n_mol
    if args.shear:
        w.flow = shear_flow(args.box, args.shear)
    return w, n


def measure(w):
    lab = cluster_labels(w)
    tails = np.isin(w.type, TAILS)
    lipids = np.unique(w.mol[tails])
    short = np.zeros(w.n_mol, bool)
    short[w.mol[w.type == STAIL]] = True
    size = np.bincount(lab[lipids], minlength=lab.max() + 1)
    big = size.argmax()
    vol, _ = water_pockets3(w)
    mask, c, _ = biggest(w)
    d = offset(w, w._wrapped()[mask & tails], c)
    ev = np.sort(np.linalg.eigvalsh(np.cov(d.T)))
    return {"step": w.t, "biggest": int(size[big]), "big_short": int(short[lipids[lab[lipids] == big]].sum()),
            "clusters30": int((size >= 30).sum()), "loose": int((size[lab[lipids]] == 1).sum()),
            "pockets": len(vol), "volumes": " ".join(str(int(v)) for v in sorted(vol, reverse=True)),
            "length": round(float(np.sqrt(ev[2] / ev[0])), 2), "added": w.n_mol - w.n_start}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--radius", type=float, default=10.0, help="bilayer middle radius of the start vesicle")
    p.add_argument("--box", type=float, default=50.0)
    p.add_argument("--mix", type=float, default=0.5, help="fraction of short-tail lipids in the start vesicle")
    p.add_argument("--short-oily", type=float, default=1.0, help="oiliness of short tails (lower = leave more)")
    p.add_argument("--osm", type=float, default=0.03, help="sugar beads per unit volume, inside and out")
    p.add_argument("--feed", type=int, default=200, help="add lipids every N steps (0 = never)")
    p.add_argument("--batch", type=int, default=1, help="lipids per feeding")
    p.add_argument("--near", type=float, default=4.0, help="feed within this distance of the vesicle (0 = anywhere)")
    p.add_argument("--shear", type=float, default=0.0)
    p.add_argument("--kT", type=float, default=1.1)
    p.add_argument("--steps", type=int, default=200000)
    p.add_argument("--every", type=int, default=1000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", default="runs/e15")
    p.add_argument("--video", default=None)
    p.add_argument("--fps", type=int, default=6)
    args = p.parse_args()

    name = (f"r{args.radius:g}_mix{args.mix:g}_o{args.short_oily:g}_osm{args.osm:g}_f{args.feed}x{args.batch}"
            f"_n{args.near:g}_sh{args.shear:g}_s{args.seed}")
    os.makedirs(args.out, exist_ok=True)
    w, n = build(args)
    print(f"{name}: {n} lipids, {w.n} beads, sugar inside {getattr(w, 'n_sugar_in', 0)}", flush=True)
    film = None
    if args.video:
        from e12_3d import Film
        film = Film(w, args.video, fps=args.fps)
    log = open(os.path.join(args.out, name + ".csv"), "w", newline="")
    writer = csv.DictWriter(log, FIELDS)
    writer.writeheader()
    t0 = time.time()
    while w.t < args.steps:
        for _ in range(args.every // (args.feed or args.every)):
            w.step(args.feed or args.every)
            if args.feed:
                feed(w, args.batch, args.near)
        row = measure(w)
        writer.writerow(row)
        log.flush()
        if film:
            film.grab(w, f"step {w.t}  wall {row['biggest']} lipids  vesicles {row['pockets']} [{row['volumes']}]  "
                         f"length {row['length']}")
        if w.t % (10 * args.every) == 0:
            print(f"step {w.t} {row} ({time.time() - t0:.0f}s)", flush=True)
        if row["pockets"] >= 2:
            np.savez_compressed(os.path.join(args.out, f"{name}_split_{w.t}.npz"),
                                pos=w.pos, type=w.type, mol=w.mol, L=w.L)
    log.close()
    np.savez_compressed(os.path.join(args.out, name + "_final.npz"), pos=w.pos, type=w.type, mol=w.mol, L=w.L)
    if film:
        film.close()


if __name__ == "__main__":
    main()
