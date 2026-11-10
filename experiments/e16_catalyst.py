"""e16: a catalyst inside: does the vesicle that makes two-chain lipid win?

Two fatty-acid vesicles side by side (A at x = Lx/4, B at 3 Lx/4), same size,
same sugar inside and out. A also holds --cats catalyst beads ("cat": a short
peptide / ribozyme, too big to leave); B holds the same number of extra sugar
beads instead (same osmotic push, no chemistry). The catalyst joins two fatty
acids of A's inner wall into one two-chain lipid (LipidJoin in reactions.py).
Lab result to compare (Budin & Szostak 2011; Adamala & Szostak 2013): two-chain
lipids stay put, so the vesicle that makes them holds its fatty acids and takes
up loose ones; its neighbour shrinks. Nothing is coded to make A grow.

Logged every --every steps, for A and B (the lipid cluster lining each one's
marker beads: A's catalysts, B's extra sugar): fatty-acid chains in the wall,
two-chain lipids in it, sealed pockets; loose lipids; joins so far.

    python experiments/e16_catalyst.py --steps 200000 --video runs/video/e16_s0.mp4
    python experiments/e16_catalyst.py --cats 0      # control: no catalyst anywhere
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
from alife.chem.elements import INDEX   # noqa: E402
from alife.chem.reactions import LipidJoin   # noqa: E402
from alife.chem.world import TAILS, World   # noqa: E402
from e15_split3d import feed, mixed_vesicle, offset, scatter_fast, shear_flow   # noqa: E402

FIELDS = ["step", "A_chains", "A_double", "B_chains", "B_double", "same_cluster", "pockets", "volumes",
          "loose", "joins", "added"]
HEADS = (INDEX["head"], INDEX["fa_head"])


def build(args):
    box = (args.box_x, args.box, args.box)
    w = World(box=box, kT=args.kT, seed=args.seed, dim=3)
    w.rules.set_oily("stail", args.short_oily)
    L = np.array(box)
    centres = [np.array([L[0] / 4, L[1] / 2, L[2] / 2]), np.array([3 * L[0] / 4, L[1] / 2, L[2] / 2])]
    n = sum(mixed_vesicle(w, c, args.radius, args.mix, args.charged) for c in centres)
    r_in, outer = args.radius - 3.5, args.radius + 2.5
    inside = [lambda p, c=c: np.linalg.norm(offset(w, p, c), axis=1) < r_in for c in centres]
    n_in = int(round(args.osm * 4 / 3 * np.pi * r_in ** 3))
    for keep in inside:
        scatter_fast(w, n_in, keep)
    scatter_fast(w, args.cats, inside[0], name="cat")            # A: catalyst
    scatter_fast(w, args.cats, inside[1])                        # B: same number of plain sugar beads
    w.marker_B = np.zeros(w.n, bool)
    w.marker_B[-args.cats:] = args.cats > 0
    n_out = int(round(args.osm * (np.prod(L) - 2 * 4 / 3 * np.pi * (outer + 1.0) ** 3)))
    scatter_fast(w, n_out, lambda p: np.all([np.linalg.norm(offset(w, p, c), axis=1) > outer + 1.0
                                             for c in centres], axis=0))
    w.marker_B = np.concatenate([w.marker_B, np.zeros(w.n - len(w.marker_B), bool)])
    w.step(300, max_force=50)
    w.t = 0
    w.n_start = w.n_mol
    if args.shear:
        w.flow = shear_flow(args.box, args.shear)
    return w, n


def wall_of(w, markers, lab, lipid_beads):
    """Lipid cluster lining the marker beads: majority vote of each marker's nearest lipid bead."""
    if not markers.any():
        return -1
    tree = cKDTree(w._wrapped()[lipid_beads], boxsize=w.L)
    _, k = tree.query(w._wrapped()[markers])
    return int(np.bincount(lab[w.mol[lipid_beads][k]]).argmax())


def measure(w, rx):
    lab = cluster_labels(w)
    tails = np.isin(w.type, TAILS)
    lipids = np.unique(w.mol[tails])
    lipid_beads = np.flatnonzero(np.isin(w.mol, lipids))
    size = np.bincount(lab[lipids], minlength=lab.max() + 1)
    heads = np.isin(w.type, HEADS)
    double = np.zeros(w.n_mol, bool)
    double[w.mol[w.bonds[rx.head_bonds(w)][:, 0]]] = True
    row = {"step": w.t}
    cat = w.type == INDEX["cat"]
    walls = {"A": wall_of(w, cat, lab, lipid_beads), "B": wall_of(w, w.marker_B, lab, lipid_beads)}
    for k, c in walls.items():
        mine = heads & (lab[w.mol] == c) if c >= 0 else np.zeros(w.n, bool)
        row[f"{k}_chains"] = int(mine.sum())                         # one head per fatty-acid chain
        row[f"{k}_double"] = int(np.unique(w.mol[mine & double[w.mol]]).size)
    row["same_cluster"] = int(walls["A"] == walls["B"] and walls["A"] >= 0)
    vol, _ = water_pockets3(w)
    row["pockets"] = len(vol)
    row["volumes"] = " ".join(str(int(v)) for v in sorted(vol, reverse=True))
    row["loose"] = int((size[lab[lipids]] == 1).sum())
    row["joins"] = rx.joins - rx.breaks
    row["added"] = w.n_mol - w.n_start
    return row


class Film:
    """Left: a thin slice through both vesicles (z middle), whole box in x-y. Right: wall size over time."""

    def __init__(self, w, path, fps=6):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from alife.viz import Movie

        class V:
            pass
        self.fig, (self.a, self.b) = plt.subplots(1, 2, figsize=(13, 5), gridspec_kw={"width_ratios": [2, 1]})
        v = V()
        v.fig = self.fig
        self.movie = Movie(v, path, fps=fps)
        self.hist = []

    def grab(self, w, row, rx):
        p = w._wrapped()
        L = np.broadcast_to(w.L, (3,))
        cut = np.abs(p[:, 2] - L[2] / 2) < 1.5
        tails, heads = np.isin(w.type, TAILS), np.isin(w.type, HEADS)
        double = np.zeros(w.n_mol, bool)
        double[w.mol[w.bonds[rx.head_bonds(w)][:, 0]]] = True
        two = double[w.mol]
        self.a.clear()
        for m, col, sz in ((w.type == INDEX["sugar"], "#DDDDDD", 4), (tails & ~two, "#EF9F27", 6),
                           (heads & ~two, "#D85A30", 8), (two & (tails | heads), "#7B2FBE", 10),
                           (w.type == INDEX["cat"], "#2E8B57", 22)):
            self.a.scatter(p[cut & m, 0], p[cut & m, 1], s=sz, c=col)
        cat = w.type == INDEX["cat"]   # catalysts out of the slice too (faint), so you see where A is
        self.a.scatter(p[cat & ~cut, 0], p[cat & ~cut, 1], s=10, c="#2E8B57", alpha=0.3)
        self.a.set_xlim(0, L[0]), self.a.set_ylim(0, L[1]), self.a.set_aspect("equal")
        self.a.set_title("slice through both vesicles: A (left) has catalyst (green); "
                         "purple = two-chain lipid", fontsize=9)
        self.hist.append((row["step"] / 1000, row["A_chains"], row["B_chains"], row["A_double"]))
        h = np.array(self.hist)
        self.b.clear()
        self.b.plot(h[:, 0], h[:, 1], c="#2E8B57", label="A wall (catalyst)")
        self.b.plot(h[:, 0], h[:, 2], c="#888888", label="B wall (no catalyst)")
        self.b.plot(h[:, 0], 2 * h[:, 3], c="#7B2FBE", ls="--", label="chains in two-chain lipids (A)")
        self.b.set_xlabel("step (thousands)"), self.b.set_ylabel("fatty-acid chains")
        self.b.legend(fontsize=8)
        self.fig.suptitle(f"step {w.t}", y=0.02, va="bottom", fontsize=10)
        self.movie.grab()

    def close(self):
        self.movie.close()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--radius", type=float, default=10.0)
    p.add_argument("--box", type=float, default=48.0, help="y and z size")
    p.add_argument("--box-x", type=float, default=80.0, help="x size (the two vesicles sit along x)")
    p.add_argument("--mix", type=float, default=0.5, help="fraction of fatty acids (short tail) in the walls")
    p.add_argument("--short-oily", type=float, default=0.85)
    p.add_argument("--charged", type=float, default=0.5, help="fraction of fatty-acid heads charged COO-")
    p.add_argument("--osm", type=float, default=0.15)
    p.add_argument("--cats", type=int, default=20, help="catalyst beads in A (B gets as many sugar beads)")
    p.add_argument("--p-join", type=float, default=0.2, help="join chance per check with a catalyst near")
    p.add_argument("--no-cat", type=float, default=0.0, help="join chance factor without a catalyst")
    p.add_argument("--p-break", type=float, default=0.0, help="chance per check that a two-chain lipid splits")
    p.add_argument("--feed", type=int, default=0, help="add --batch loose fatty acids every N steps (0 = never)")
    p.add_argument("--batch", type=int, default=1)
    p.add_argument("--shear", type=float, default=0.0)
    p.add_argument("--kT", type=float, default=1.1)
    p.add_argument("--steps", type=int, default=200000)
    p.add_argument("--every", type=int, default=1000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", default="runs/e16")
    p.add_argument("--video", default=None)
    p.add_argument("--fps", type=int, default=6)
    args = p.parse_args()

    name = (f"cats{args.cats}_p{args.p_join:g}_nc{args.no_cat:g}_br{args.p_break:g}_mix{args.mix:g}"
            f"_q{args.charged:g}_f{args.feed}x{args.batch}_sh{args.shear:g}_s{args.seed}")
    os.makedirs(args.out, exist_ok=True)
    w, n = build(args)
    rx = LipidJoin(p_join=args.p_join, no_cat=args.no_cat, p_break=args.p_break)
    print(f"{name}: {n} lipids, {w.n} beads", flush=True)
    film = None
    if args.video:
        film = Film(w, args.video, fps=args.fps)
    log = open(os.path.join(args.out, name + ".csv"), "w", newline="")
    writer = csv.DictWriter(log, FIELDS)
    writer.writeheader()
    t0 = time.time()
    while w.t < args.steps:
        for _ in range(args.every // rx.every):
            w.step(rx.every)
            rx.react(w)
            if args.feed and w.t % args.feed == 0:
                feed(w, args.batch, 0, clear=1.0, charged=args.charged)
        row = measure(w, rx)
        writer.writerow(row)
        log.flush()
        if film:
            film.grab(w, row, rx)
        if w.t % (10 * args.every) == 0:
            print(f"step {w.t} {row} ({time.time() - t0:.0f}s)", flush=True)
    log.close()
    np.savez_compressed(os.path.join(args.out, name + "_final.npz"), pos=w.pos, type=w.type, mol=w.mol, L=w.L,
                        bonds=w.bonds, marker_B=w.marker_B)
    if film:
        film.close()


if __name__ == "__main__":
    main()
