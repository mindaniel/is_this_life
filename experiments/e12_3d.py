"""e12: lipids in 3D (World(dim=3), the Cooke-Kremer-Deserno model as published).

In 2D a vesicle is a ring and pops easily; real vesicles are hollow balls.
Starts:
    disc     a flat round bilayer patch (--radius). Its open edge costs
             energy; a big enough disc curls up and closes into a vesicle
    gas      lipids scattered at random (--lipids): self-assembly from scratch
    vesicle  a ready-made vesicle (--radius), to test how long it lasts

Logs clusters and sealed water pockets (vesicle insides, 3D detector).
--video films a 3D view and a slice through the middle (a vesicle shows as a ring).

    python experiments/e12_3d.py --start disc --radius 16 --steps 100000 --video runs/video/e12_disc.mp4
"""
import argparse
import csv
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from alife.chem.analysis import tail_clusters, water_pockets3   # noqa: E402
from alife.chem.molecules import bilayer_disc3, gas3, vesicle3   # noqa: E402
from alife.chem.world import HEAD, TAILS, World   # noqa: E402


def build(args):
    w = World(box=args.box, kT=args.kT, seed=args.seed, dim=3)
    mid = [args.box / 2] * 3
    if args.start == "disc":
        n = bilayer_disc3(w, mid, args.radius)
    elif args.start == "vesicle":
        n = vesicle3(w, mid, args.radius)
    else:
        gas3(w, args.lipids)
        n = args.lipids
    w.step(300, max_force=50)
    return w, n


def shape(w):
    """How round the biggest cluster is: 0 = flat sheet, ~1 = ball (from the spread along 3 axes)."""
    tails = np.isin(w.type, TAILS)
    from alife.chem.analysis import cluster_labels
    lab = cluster_labels(w)[w.mol]
    big = np.bincount(lab[tails]).argmax()
    p = w.pos[(lab == big) & tails]
    ang = 2 * np.pi * p / w.L
    c = np.angle(np.exp(1j * ang).mean(0)) % (2 * np.pi) / (2 * np.pi) * w.L
    d = p - c
    d -= w.L * np.round(d / w.L)
    ev = np.sort(np.linalg.eigvalsh(np.cov(d.T)))
    return ev[0] / ev[2]


class Film:
    def __init__(self, w, path, fps=15):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from alife.viz import Movie

        class V:
            pass
        self.fig = plt.figure(figsize=(12, 6))
        self.a3 = self.fig.add_subplot(1, 2, 1, projection="3d")
        self.a2 = self.fig.add_subplot(1, 2, 2)
        v = V()
        v.fig = self.fig
        self.movie = Movie(v, path, fps=fps)

    def grab(self, w, title):
        """Close-up of the biggest cluster: left = 3D cut in half (shows the hollow inside), right = thin slice."""
        p = w._wrapped()
        head = w.type == HEAD
        from alife.chem.analysis import cluster_labels
        lab = cluster_labels(w)[w.mol]
        big = np.bincount(lab[~head]).argmax()
        mine = lab == big
        c = np.angle(np.exp(2j * np.pi * p[mine] / w.L).mean(0)) % (2 * np.pi) / (2 * np.pi) * w.L
        d = p - c
        d -= w.L * np.round(d / w.L)
        axes = np.linalg.eigh(np.cov(d[mine].T))[1]   # columns: thinnest, middle, widest direction
        u, depth, v = d @ axes[:, 2], d @ axes[:, 1], d @ axes[:, 0]
        r = min(np.abs(np.c_[u, depth, v][mine]).max() + 2, w.L / 2)   # zoom to the cluster
        half = mine & (depth < 0)   # front half cut away so you see inside
        self.a3.clear()
        self.a3.scatter(u[half & ~head], depth[half & ~head], v[half & ~head], s=4, c="#EF9F27", alpha=0.5, depthshade=True)
        self.a3.scatter(u[half & head], depth[half & head], v[half & head], s=7, c="#D85A30", depthshade=True)
        for f in (self.a3.set_xlim, self.a3.set_ylim, self.a3.set_zlim):
            f(-r, r)
        self.a3.view_init(elev=15, azim=90)
        self.a3.set_title("biggest cluster, cut in half (heads red, tails orange)")
        cut = np.abs(depth) < 1.5
        self.a2.clear()
        self.a2.scatter(u[cut & ~head], v[cut & ~head], s=10, c="#EF9F27")
        self.a2.scatter(u[cut & head], v[cut & head], s=14, c="#D85A30")
        self.a2.set_xlim(-r, r)
        self.a2.set_ylim(-r, r)
        self.a2.set_aspect("equal")
        self.a2.set_title("thin slice through it: sheet = line, vesicle = ring")
        self.fig.suptitle(title, y=0.02, va="bottom", fontsize=10)
        self.movie.grab()

    def close(self):
        self.movie.close()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--start", choices=("disc", "gas", "vesicle"), default="disc")
    p.add_argument("--radius", type=float, default=16.0)
    p.add_argument("--lipids", type=int, default=2000)
    p.add_argument("--box", type=float, default=50.0)
    p.add_argument("--kT", type=float, default=1.1, help="Cooke model: fluid bilayer ~0.9-1.3 with wc 1.6")
    p.add_argument("--steps", type=int, default=100000)
    p.add_argument("--every", type=int, default=1000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", default="runs/e12")
    p.add_argument("--video", default=None)
    args = p.parse_args()
    name = f"{args.start}_r{args.radius:g}_kT{args.kT:g}_s{args.seed}" if args.start != "gas" else \
        f"gas_n{args.lipids}_box{args.box:g}_kT{args.kT:g}_s{args.seed}"
    os.makedirs(args.out, exist_ok=True)
    w, n = build(args)
    print(f"{name}: {n} lipids, {w.n} beads, box {args.box}", flush=True)
    film = Film(w, args.video) if args.video else None
    rows = []
    t0 = time.time()
    while w.t < args.steps:
        w.step(args.every)
        cl = tail_clusters(w)
        vol, _ = water_pockets3(w)
        row = dict(step=w.t, clusters=len(cl), largest=int(cl[0]), roundness=round(shape(w), 3),
                   pockets=len(vol), volume=int(vol[0]) if len(vol) else 0, energy=round(w.energy / n, 3))
        rows.append(row)
        if film:
            film.grab(w, f"{name}  step {w.t}  largest {row['largest']}  sealed pockets {len(vol)}"
                         f"{'  volume ' + str(row['volume']) if len(vol) else ''}")
        if len(rows) % 10 == 0:
            print(f"step {w.t} {row} ({time.time() - t0:.0f}s)", flush=True)
    with open(os.path.join(args.out, name + ".csv"), "w", newline="") as f:
        out = csv.DictWriter(f, fieldnames=list(rows[0]))
        out.writeheader()
        out.writerows(rows)
    np.savez_compressed(os.path.join(args.out, name + "_final.npz"), pos=w.pos, type=w.type, mol=w.mol, L=w.L)
    if film:
        film.close()


if __name__ == "__main__":
    main()
