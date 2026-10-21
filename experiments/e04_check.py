"""Check e04 vesicle detections: are they real, and how long do they live?

Reads runs/<name>/seed*/log.csv (pockets column) and the saved frames
(frame_*.npz, written whenever a vesicle was detected).

  lifetimes: a pocket at the next log within `--match` units of one now is
             the same vesicle. Prints how many distinct vesicles and how long.
  gallery:   zoomed pictures of every detected pocket, to judge by eye.

    python experiments/e04_check.py runs/e04_world
    python experiments/e04_check.py runs/e04_fatty --max-crops 48
"""
import argparse
import csv
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from alife.chem.analysis import lining, water_pockets
from alife.chem.elements import COLOR
from alife.chem.world import TAIL, World


def read_pockets(log_path):
    """[(step, [(x, y, area), ...]), ...] from one log."""
    rows = []
    with open(log_path) as f:
        for r in csv.DictReader(f):
            ps = [tuple(float(v) for v in p.split(",")) for p in r["pockets"].split(";") if p]
            rows.append((int(r["step"]), ps))
    return rows


def track(rows, L, match):
    """Link pockets between consecutive logs. Returns list of tracks [(step, x, y, area), ...]."""
    open_tracks, done = [], []
    for step, ps in rows:
        new = []
        taken = set()
        for tr in open_tracks:
            _, x0, y0, _ = tr[-1]
            best, bd = None, match
            for k, (x, y, a) in enumerate(ps):
                if k in taken:
                    continue
                d = np.array([x - x0, y - y0])
                d -= L * np.round(d / L)
                if np.hypot(*d) < bd:
                    best, bd = k, np.hypot(*d)
            if best is None:
                done.append(tr)
            else:
                taken.add(best)
                tr.append((step, *ps[best]))
                new.append(tr)
        for k, (x, y, a) in enumerate(ps):
            if k not in taken:
                new.append([(step, x, y, a)])
        open_tracks = new
    return done + open_tracks


def load_frame(path):
    """A World holding a saved frame (positions only), for the detectors."""
    f = np.load(path)
    w = World(box=float(f["L"]))
    w.pos, w.type, w.mol, w.bonds = f["pos"].astype(float), f["type"], f["mol"], f["bonds"]
    w.n_mol = int(w.mol.max()) + 1
    w.vel = np.zeros_like(w.pos)
    w.frozen = np.zeros(len(w.pos), dtype=bool)
    return w


def score_at(world, x, y):
    """(lining score, coverage) of the detected pocket nearest (x, y)."""
    _, where, masks = water_pockets(world, centers=True, masks=True)
    if not where:
        return 0.0, 0.0
    d = np.array(where) - np.array([x, y])
    d -= world.L * np.round(d / world.L)
    return lining(world, masks[int(np.argmin(np.hypot(*d.T)))])


def draw_crop(ax, frame, x, y, half):
    """Beads and bonds within `half` of (x, y), periodic, centred."""
    from matplotlib.collections import LineCollection
    from matplotlib.patches import Circle

    L = float(frame["L"])
    p = frame["pos"] - np.array([x, y])
    p -= L * np.round(p / L)
    near = np.all(np.abs(p) < half + 2, axis=1)
    b = frame["bonds"]
    keep = near[b[:, 0]] & near[b[:, 1]]
    seg = np.stack([p[b[keep, 0]], p[b[keep, 1]]], axis=1)
    seg = seg[np.all(np.abs(seg[:, 0] - seg[:, 1]) < L / 2, axis=1)]
    ax.add_collection(LineCollection(seg, colors="#888780", linewidths=0.6))
    t = frame["type"][near]
    ax.scatter(*p[near].T, c=np.array(COLOR)[t], s=np.where(t == TAIL, 8, 14), zorder=3)
    ax.add_patch(Circle((0, 0), 3, fill=False, ec="#E24B4A", lw=1, ls="--", zorder=4))
    ax.set_xlim(-half, half)
    ax.set_ylim(-half, half)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor("#E6F1FB")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("run", help="e.g. runs/e04_world")
    p.add_argument("--box", type=float, default=192.0)
    p.add_argument("--match", type=float, default=8.0, help="max drift between logs for the same vesicle")
    p.add_argument("--half", type=float, default=14.0, help="crop half-width")
    p.add_argument("--max-crops", type=int, default=36)
    p.add_argument("--min-lining", type=float, default=0.8, help="heads facing in (1 = all)")
    p.add_argument("--min-cover", type=float, default=0.9, help="lining lipids per unit of pocket edge")
    args = p.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    all_tracks = []
    for seed_dir in sorted(glob.glob(os.path.join(args.run, "seed*"))):
        rows = read_pockets(os.path.join(seed_dir, "log.csv"))
        if not rows:
            continue
        every = rows[1][0] - rows[0][0] if len(rows) > 1 else 1
        tracks = track(rows, args.box, args.match)
        life = [(tr[-1][0] - tr[0][0] + every) for tr in tracks]
        n_det = sum(len(ps) for _, ps in rows)
        print(f"{os.path.basename(seed_dir)}: {rows[-1][0]} steps, {n_det} detections, "
              f"{len(tracks)} vesicles, lifetimes {sorted(life, reverse=True)[:8]}")
        all_tracks += [(seed_dir, tr, t) for tr, t in zip(tracks, life)]

    if not all_tracks:
        print("no vesicles")
        return
    life = np.array([t for _, _, t in all_tracks])
    print(f"total {len(all_tracks)} vesicles; median life {np.median(life):.0f} steps, "
          f"longest {life.max()}, living >= 5000 steps: {np.sum(life >= 5000)}")

    # score each vesicle at its middle frame: does the wall look like a membrane?
    scores = np.zeros((len(all_tracks), 2))
    for k, (seed_dir, tr, _) in enumerate(all_tracks):
        step, x, y, _ = tr[len(tr) // 2]
        path = os.path.join(seed_dir, f"frame_{step:06d}.npz")
        if os.path.exists(path):
            scores[k] = score_at(load_frame(path), x, y)
    real = (scores[:, 0] >= args.min_lining) & (scores[:, 1] >= args.min_cover)
    print(f"membrane-lined (lining >= {args.min_lining}, cover >= {args.min_cover}): "
          f"{np.sum(real)} of {len(all_tracks)}; their lifetimes {sorted(life[real], reverse=True)}")

    # gallery: best-covered first, each shown at its middle frame
    order = np.argsort(-scores[:, 1])[:args.max_crops]
    n = len(order)
    cols = min(6, n)
    rows_n = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows_n, cols, figsize=(2.6 * cols, 2.8 * rows_n), squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")
    for ax, k in zip(axes.ravel(), order):
        seed_dir, tr, _ = all_tracks[k]
        step, x, y, a = tr[len(tr) // 2]
        path = os.path.join(seed_dir, f"frame_{step:06d}.npz")
        if not os.path.exists(path):
            continue
        ax.axis("on")
        draw_crop(ax, np.load(path), x, y, args.half)
        ax.set_title(f"{os.path.basename(seed_dir)} t={step} area {a:.0f} life {life[k]}\n"
                     f"lining {scores[k, 0]:.2f}  cover {scores[k, 1]:.2f}",
                     fontsize=8, color="#1D9E75" if real[k] else "black")
    fig.tight_layout()
    out = os.path.join(args.run, "vesicle_gallery.png")
    fig.savefig(out, dpi=110)
    print("gallery:", out)


if __name__ == "__main__":
    main()
