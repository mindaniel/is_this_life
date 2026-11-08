"""e15 snapshot check: are 2+ pockets really separate vesicles?

Per pocket: volume, which lipid cluster walls it (touching tails = one cluster),
sugar from the start vesicle inside it (sugar molecule ids follow the lipids,
inside sugar first). Per pair of pockets: head (H) / tail (T) / sugar (s) / empty
(.) along the line between them, one letter per 0.5: one bilayer reads H TTTT H,
two vesicles pressed together read H TTTT H H TTTT H (s between if water).
Start sugar near the pockets (within 12 of a centre) vs far out (beyond 20).
With --png: a cut through the two biggest pockets (blue = start sugar).

    python experiments/e15_check.py runs/e15c/box80/*_final.npz --png check.png
"""
import argparse
import os
import sys

import numpy as np
from scipy import ndimage

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from alife.chem.analysis import cluster_labels   # noqa: E402
from alife.chem.elements import INDEX   # noqa: E402
from alife.chem.world import TAILS, World   # noqa: E402

SUGAR, HEAD = INDEX["sugar"], INDEX["head"]


def load(f, n_lipids, radius):
    z = np.load(f)
    L = float(z["L"])
    w = World(box=L, dim=3)
    w.pos, w.type, w.mol = z["pos"] % L, z["type"], z["mol"]
    w.vel, w.frozen, w.n_mol = np.zeros_like(w.pos), np.zeros(len(w.pos), bool), int(w.mol.max()) + 1
    osm = float(f.split("osm")[1].split("_")[0])
    n_in = int(round(osm * 4 / 3 * np.pi * (radius - 3.5) ** 3))   # as e15 build(): r_in = radius - 2.5 - 1
    w.start_sugar = (w.type == SUGAR) & (w.mol < n_lipids + n_in)
    return w


def pockets(w, min_volume=30, one_cluster=0.8):
    """water_pockets3 with each pocket's wall cluster, start-sugar count and centre."""
    m, r = int(np.ceil(w.L)), 1
    has_tail = np.zeros(w.n_mol, bool)
    has_tail[w.mol[np.isin(w.type, TAILS)]] = True
    lipid = has_tail[w.mol]
    lab = cluster_labels(w)
    wall, owner = np.zeros((m,) * 3, bool), np.full((m,) * 3, -1)
    idx = w._wrapped()[lipid].astype(int) % m
    wall[tuple(idx.T)] = True
    owner[tuple(idx.T)] = lab[w.mol[lipid]]
    ball = np.sqrt((np.mgrid[-r:r + 1, -r:r + 1, -r:r + 1] ** 2).sum(0)) <= r
    wall = ndimage.binary_dilation(np.pad(wall, r, mode="wrap"), ball)[r:-r, r:-r, r:-r]
    labels, n = ndimage.label(~wall)
    parent = np.arange(n + 1)

    def root(a):
        while parent[a] != a:
            a = parent[a]
        return a
    for ax in range(3):
        a_, b_ = np.take(labels, 0, axis=ax).ravel(), np.take(labels, -1, axis=ax).ravel()
        ok = (a_ > 0) & (b_ > 0)
        for a, b in set(zip(a_[ok], b_[ok])):
            parent[root(a)] = root(b)
    region = np.array([root(a) for a in range(n + 1)])[labels]
    vols = np.bincount(region.ravel(), minlength=n + 1)
    vols[0] = 0
    sea = np.argmax(vols)
    sidx = w._wrapped()[w.start_sugar].astype(int) % m
    sreg = region[tuple(sidx.T)]
    size = np.bincount(lab[np.unique(w.mol[np.isin(w.type, TAILS)])], minlength=lab.max() + 1)
    out = []
    for k in np.flatnonzero(vols >= min_volume):
        if k == sea:
            continue
        mask = region == k
        border = ndimage.binary_dilation(np.pad(mask, 2, mode="wrap"), iterations=2)[2:-2, 2:-2, 2:-2]
        walls = owner[border & (owner >= 0)]
        if not len(walls) or np.bincount(walls).max() < one_cluster * len(walls):
            continue
        top = np.bincount(walls).argmax()
        ang = 2 * np.pi * np.argwhere(mask) / m
        cen = np.angle(np.exp(1j * ang).mean(0)) % (2 * np.pi) / (2 * np.pi) * m + 0.5
        out.append({"vol": int(vols[k]), "cluster": int(top), "lipids": int(size[top]),
                    "sugar": int((sreg == k).sum()), "centre": cen})
    return sorted(out, key=lambda p: -p["vol"])


def line(w, a, b):
    """Unit vector a->b (periodic), distance, and every bead's (along, sideways) offset from a."""
    e = ((b - a + w.L / 2) % w.L) - w.L / 2
    D = np.linalg.norm(e)
    e /= D
    d = w._wrapped() - a
    d -= w.L * np.round(d / w.L)
    t = d @ e
    return e, D, t, d - np.outer(t, e)


def layers(w, a, b):
    _, D, t, perp = line(w, a, b)
    r = np.linalg.norm(perp, axis=1)
    tail, head, sug = np.isin(w.type, TAILS), w.type == HEAD, w.type == SUGAR
    s = ""
    for x in np.arange(0, D, 0.5):
        sl = (np.abs(t - x) < 0.25) & (r < 3)
        h, tt, sg = (sl & head).sum(), (sl & tail).sum(), (sl & sug).sum()
        s += "H" if h > tt and h >= sg and h else "T" if tt >= h and tt >= sg and tt else "s" if sg else "."
    return D, s


def picture(ax, w, ps, title):
    e, D, t, perp = line(w, ps[0]["centre"], ps[1]["centre"])
    tail = np.isin(w.type, TAILS)
    near = tail & (np.linalg.norm(perp, axis=1) < 30)
    _, vec = np.linalg.eigh(np.cov(perp[near].T))
    u = vec[:, -1]
    v = np.cross(e, u)
    cut = np.abs(perp @ v) < 1.2
    y = perp @ u
    for m, col, sz in ((tail, "#EF9F27", 6), (w.type == HEAD, "#D85A30", 8),
                       ((w.type == SUGAR) & ~w.start_sugar, "#BBBBBB", 6), (w.start_sugar, "#2A6FB0", 12)):
        ax.scatter(t[cut & m] - D / 2, y[cut & m], s=sz, c=col)
    ax.set_xlim(-30, 30), ax.set_ylim(-30, 30), ax.set_aspect("equal")
    ax.set_title(title, fontsize=8)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("snaps", nargs="+")
    p.add_argument("--lipids", type=int, default=2225, help="lipids in the start vesicle")
    p.add_argument("--radius", type=float, default=10.0)
    p.add_argument("--png", default=None)
    args = p.parse_args()
    shown = []
    for f in args.snaps:
        w = load(f, args.lipids, args.radius)
        ps = pockets(w)
        name = os.path.basename(f)[:-4]
        cen = np.array([q["centre"] for q in ps]) if ps else np.zeros((0, 3))
        d = w._wrapped()[w.start_sugar][:, None, :] - cen[None]
        d -= w.L * np.round(d / w.L)
        dist = np.linalg.norm(d, axis=2).min(1) if len(ps) else np.full(w.start_sugar.sum(), np.inf)
        print(f"{name}: start sugar {w.start_sugar.sum()}, in pockets {sum(q['sugar'] for q in ps)}, "
              f"within 12 of a pocket {int((dist < 12).sum())}, beyond 20 {int((dist > 20).sum())}")
        for q in ps:
            print(f"   pocket vol {q['vol']:4d}  wall cluster {q['cluster']} ({q['lipids']} lipids)  "
                  f"start sugar {q['sugar']}  centre {np.round(q['centre'])}")
        for i in range(len(ps)):
            for j in range(i + 1, len(ps)):
                D, s = layers(w, ps[i]["centre"], ps[j]["centre"])
                print(f"   {i}-{j} dist {D:4.1f}: {s}")
        if len(ps) >= 2:
            shown.append((w, ps, name))
    if args.png and shown:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, len(shown), figsize=(5 * len(shown), 5), squeeze=False)
        for ax, (w, ps, name) in zip(axes[0], shown):
            picture(ax, w, ps, f"{name}\ncut through the two biggest pockets {[q['vol'] for q in ps]}")
        fig.savefig(args.png, dpi=70, bbox_inches="tight")


if __name__ == "__main__":
    main()
