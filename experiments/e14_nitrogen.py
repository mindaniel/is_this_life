"""Experiment 14: nitrogen chemistry alone in a small box of water.

HCN and ammonia rain into ponds (lightning, impacts), next to formaldehyde.
See alife/chem/nitrogen.py for the rules: Strecker (-> glycine), cyanohydrin
side product, HCN polymer (dark, absorbs light), and hydrolysis.

Checks that the rules behave like the real chemistry:
  - ammonia picks glycine over glycolonitrile      --nh3 0
  - ammonia (base) speeds the HCN polymer          --nh3 0 --carb 0
  - dilute HCN is lost to water, concentrated
    HCN polymerises                                --box 80 (4x more water)

    python experiments/e14_nitrogen.py --name base --seed 0
    python experiments/e14_nitrogen.py --plot      # compare all runs
    python experiments/e14_nitrogen.py --video runs/video/e14_base.mp4 --log-every 250 --out runs/e14_film
"""
import argparse
import csv
import glob
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from alife.chem import World, add_ions   # noqa: E402
from alife.chem.nitrogen import CARB, HCN, NH3, Nitrogen, census   # noqa: E402

MAX_LEN = 8   # HCN chains 2..7 logged separately, 8 = 8 or longer
DARK = 5      # HCN chains this long or longer count as dark polymer (pigment-like)
FIELDS = (["step", "hcn", "nh3", "carb", "glycine", "cyano"] + [f"h{n}" for n in range(2, MAX_LEN + 1)]
          + ["dark_hcn", "hydrolysed"])


def build(args):
    world = World(box=args.box, kT=args.kT, seed=args.seed)
    world.repel_inside = (CARB, HCN, NH3)   # floppy products: their beads must not overlap
    for name, count in (("hcn", args.hcn), ("nh3", args.nh3), ("carb", args.carb)):
        if count:
            add_ions(world, name, count)
    world.step(200, max_force=50.0)
    world.t = 0
    return world, Nitrogen(every=args.every)


def info(world, chem):
    c = census(world)
    h = c.pop("hcn_chains")
    row = dict(step=world.t, **c)
    for n in range(2, MAX_LEN):
        row[f"h{n}"] = int(h[n]) if n < len(h) else 0
    row[f"h{MAX_LEN}"] = int(h[MAX_LEN:].sum())
    row["dark_hcn"] = int((h * np.arange(len(h)))[DARK:].sum())
    row["hydrolysed"] = chem.hydrolysed
    return row


RULES = "\n".join([
    "Rules (beads only react when they touch):",
    "1  formaldehyde + HCN + ammonia -> glycine (amino acid)",
    "2  formaldehyde + HCN, no ammonia -> glycolonitrile (dead end)",
    "3  HCN chain end + HCN chain end -> longer chain",
    "    (faster with ammonia nearby; first HCN+HCN is slow)",
    "4  lone HCN + water -> formic acid + ammonia (lost)"])


class Film:
    """Left: the box, each bead coloured by what it has become, bonds as lines.
    Right: counts over time and the rules."""

    def __init__(self, world, path):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from alife.viz import Movie

        self.fig = plt.figure(figsize=(14, 7))
        self.ax = self.fig.add_axes([0.02, 0.05, 0.47, 0.88])
        self.plot = self.fig.add_axes([0.56, 0.46, 0.41, 0.46])
        self.fig.text(0.66, 0.06, RULES, fontsize=9, family="monospace", va="bottom")
        self.size = (7 * 0.47 * 72 / world.L) ** 2   # bead diameter 1 -> marker area in points^2
        self.hist, self.movie = [], Movie(self, path, fps=15)

    def grab(self, world, row):
        t, ax = world.type, self.ax
        size = np.bincount(world.mol, minlength=world.n_mol)[world.mol]
        n_h = np.bincount(world.mol[t == HCN], minlength=world.n_mol)[world.mol]
        role = np.full(world.n, "other", dtype=object)
        role[(size == 1) & (t == HCN)] = "free HCN"
        role[(size == 1) & (t == NH3)] = "ammonia"
        role[(size == 1) & (t == CARB)] = "formaldehyde"
        role[size == 3] = "glycine"
        role[(size == 2) & (n_h == 1)] = "glycolonitrile"
        role[(n_h == size) & (size > 1)] = "HCN chain"
        look = {"free HCN": ("#2A6FB0", 0.8), "ammonia": ("#9DCBEA", 0.75), "formaldehyde": ("#8E5BB5", 0.8),
                "glycine": ("#1F9E4A", 1.0), "glycolonitrile": ("#AAAAAA", 0.8), "HCN chain": ("#3B2314", 1.0)}
        ax.clear()
        p = world._wrapped()
        for a, b in world.bonds:
            d = p[b] - p[a]
            if np.abs(d).max() < world.L / 2:   # skip bonds drawn across the wrap
                ax.plot(*np.c_[p[a], p[b]], color=look.get(role[a], ("#888888",))[0], lw=2.5, zorder=1)
        for name, (color, d) in look.items():
            m = role == name
            ax.scatter(p[m, 0], p[m, 1], s=self.size * d * d, c=color, edgecolors="none", zorder=2,
                       label=f"{name} ({m.sum() if name != 'glycine' else m.sum() // 3})")
        ax.set_xlim(0, world.L)
        ax.set_ylim(0, world.L)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"step {world.t}   (box ~{world.L * 0.9:.0f} nm across)", loc="left")
        if getattr(self, "key", None):
            self.key.remove()
        self.key = self.fig.legend(*ax.get_legend_handles_labels(), loc="lower left", bbox_to_anchor=(0.5, 0.04),
                                   fontsize=10, frameon=False, title="what each bead is now")
        self.hist.append(row)
        self.plot.clear()
        steps = [r["step"] for r in self.hist]
        for key, label, color in (("hcn", "free HCN", "#2A6FB0"), ("dark_hcn", "HCN in long chains (5+): dark polymer", "#3B2314"),
                                  ("glycine", "glycine", "#1F9E4A"), ("cyano", "glycolonitrile", "#AAAAAA"),
                                  ("hydrolysed", "HCN lost to water", "#D85A30")):
            self.plot.plot(steps, [r[key] for r in self.hist], color=color, lw=2, label=label)
        self.plot.set_xlim(0, max(steps[-1], 1))
        self.plot.set_xlabel("step")
        self.plot.set_ylabel("molecules (HCN beads for chains)")
        self.plot.legend(fontsize=9)
        self.movie.grab()


def run(args):
    world, chem = build(args)
    film = Film(world, args.video) if args.video else None
    out = os.path.join(args.out, f"{args.name}_seed{args.seed}")
    os.makedirs(out, exist_ok=True)
    print(f"{args.name}: {args.hcn} HCN, {args.nh3} NH3, {args.carb} formaldehyde, box {args.box:g}, kT {args.kT}",
          flush=True)
    with open(os.path.join(out, "log.csv"), "w", newline="") as fh:
        writer = csv.DictWriter(fh, FIELDS)
        writer.writeheader()
        start = time.time()
        while world.t < args.steps:
            for _ in range(args.log_every // chem.every):
                world.step(chem.every)
                chem.react(world)
            row = info(world, chem)
            if film:
                film.grab(world, row)
            writer.writerow(row)
            fh.flush()
            if world.t % (args.steps // 10) == 0:
                print(f"step {world.t} {row} ({time.time() - start:.0f}s)", flush=True)
    if film:
        film.movie.close()
    np.savez_compressed(os.path.join(out, "final.npz"), pos=world.pos, type=world.type, mol=world.mol,
                        bonds=world.bonds, L=world.L)


def plot(args):
    """Where the HCN went over time: still free, polymer, glycine, lost to water."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    runs = {}
    for path in sorted(glob.glob(os.path.join(args.out, "*_seed*", "log.csv"))):
        rows = list(csv.DictReader(open(path)))
        if rows:
            runs.setdefault(os.path.basename(os.path.dirname(path)).split("_seed")[0], []).append(rows)
    keys = [("hcn", "free HCN"), ("dark_hcn", f"HCN in dark polymer (chains {DARK}+)"),
            ("glycine", "glycine"), ("hydrolysed", "HCN lost to water")]
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.5))
    colors = plt.cm.tab10.colors
    for c, (name, logs) in enumerate(runs.items()):
        for k, rows in enumerate(logs):
            t = [int(r["step"]) for r in rows]
            for ax, (key, _) in zip(axes, keys):
                ax.plot(t, [int(r[key]) for r in rows], color=colors[c % 10], lw=0.8, alpha=0.7,
                        label=name if k == 0 else None)
        end = {key: np.mean([int(rows[-1][key]) for rows in logs]) for key, _ in keys + [("cyano", "")]}
        print(f"{name:12s} {len(logs)} runs, at end: free HCN {end['hcn']:.0f}  dark {end['dark_hcn']:.0f}  "
              f"glycine {end['glycine']:.0f}  glycolonitrile {end['cyano']:.0f}  lost {end['hydrolysed']:.0f}")
    for ax, (_, title) in zip(axes, keys):
        ax.set_title(title)
        ax.set_xlabel("step")
    axes[0].legend(fontsize=8)
    path = os.path.join(args.out, "compare.png")
    fig.savefig(path, dpi=100, bbox_inches="tight")
    print("saved", path)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--hcn", type=int, default=300)
    p.add_argument("--nh3", type=int, default=100)
    p.add_argument("--carb", type=int, default=100, help="formaldehyde")
    p.add_argument("--box", type=float, default=40.0)
    p.add_argument("--kT", type=float, default=0.6)
    p.add_argument("--every", type=int, default=50, help="steps between reaction checks")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--steps", type=int, default=100000)
    p.add_argument("--log-every", type=int, default=500)
    p.add_argument("--name", default="base")
    p.add_argument("--plot", action="store_true")
    p.add_argument("--out", default="runs/e14_nitrogen")
    p.add_argument("--video", default=None, help="film the run (.mp4 or .gif)")
    args = p.parse_args()
    plot(args) if args.plot else run(args)


if __name__ == "__main__":
    main()
