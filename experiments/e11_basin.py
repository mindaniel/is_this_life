"""e11: the river basin in real units (landscape level, alife/basin.py).

Ponds of 4-25 m on hills, hourly rain/evaporation/streams/floods for years.
Asks: which ponds ever reach the lipid concentration needed for vesicles
(CVC), how often they dry and rewet, and whether sugar is there at the same
time. The answers say what conditions the bead world should simulate.

    python experiments/e11_basin.py --years 10 --seeds 0 1 2
    python experiments/e11_basin.py --years 2 --video runs/video/e11_basin.mp4
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from alife.basin import Basin   # noqa: E402

COLOR = {"plain": "#4a90c8", "hot": "#d9603b", "rocky": "#7a7a7a"}


def draw_map(ax, b, lip=None, flow=None, grow=1.0, crop=False):
    """Ground contours, ponds (size = water now, colour = lipid / CVC), streams.

    grow > 1 draws ponds bigger than life (they are small on the map);
    crop=True zooms to where the ponds are and writes each pond's numbers.
    """
    import matplotlib.colors as mc
    from matplotlib import cm
    from matplotlib.patches import Circle
    g = np.linspace(0, b.size, 120)
    X, Y = np.meshgrid(g, g)
    Z = b.height(np.stack([X, Y], -1).reshape(-1, 2)).reshape(X.shape)
    ax.contourf(X, Y, Z, 14, cmap="YlOrBr_r", alpha=0.45)
    ax.contour(X, Y, Z, 14, colors="k", linewidths=0.3, alpha=0.3)
    for i, j in enumerate(b.streams):
        if j >= 0:
            a, c = b.centers[i], b.centers[j]
            e = (c - a) / np.hypot(*(c - a))
            ax.annotate("", c - e * grow * b.radius[j], a + e * grow * b.radius[i],
                        arrowprops=dict(arrowstyle="->", color="#1f5fa0", lw=1.2))
    if flow is not None:   # floods this moment: thin dashed lines
        for i, j in zip(*np.nonzero(flow > 1.0)):
            if b.streams[i] != j:
                ax.plot(*b.centers[[i, j]].T, "--", color="#1f5fa0", lw=0.8)
    area = b.area() if lip is not None else b.area_max
    norm = mc.LogNorm(1e-3, 10)
    for i in range(len(b.radius)):
        ax.add_patch(Circle(
            b.centers[i], grow * b.radius[i], fill=False, ls=":", ec=COLOR[b.kind[i]], lw=1))
        r = grow * np.sqrt(area[i] / np.pi)
        fc = "#cfe3f5" if lip is None else cm.viridis(norm(max(lip[i] / b.cvc, 1e-3)))
        ax.add_patch(Circle(
            b.centers[i], r, fc=fc, ec=COLOR[b.kind[i]], lw=1.5))
        if crop and lip is not None:
            ax.text(b.centers[i, 0], b.centers[i, 1] - grow * b.radius[i] - 2,
                    f"{i} {b.kind[i]}\n{100 * b.v[i] / b.v_max[i]:.0f}% full\nlipid {lip[i]:.2g} mM",
                    ha="center", va="top", fontsize=8, bbox=dict(fc="white", ec="none", alpha=0.7, pad=1))
        else:
            ax.text(*b.centers[i], f"{i}\n{b.kind[i]}", ha="center", va="center", fontsize=7)
    if crop:
        pad = grow * b.radius.max() + 12
        lo = np.maximum(b.centers.min(0) - pad, 0)
        hi = np.minimum(b.centers.max(0) + pad, b.size)
        ax.set_xlim(lo[0], hi[0])
        ax.set_ylim(lo[1], hi[1])
    else:
        ax.set_xlim(0, b.size)
        ax.set_ylim(0, b.size)
    ax.set_aspect("equal")
    ax.set_xlabel("metres")
    return norm


def summary(b, log, years):
    n = len(b.radius)
    print(f"  pond kind   size(m) depth  wet%  dry->wet/yr  vesicle h/yr  rewet-with-vesicles  "
          f"vesicles+sugar h/yr  max lipid mM  max sugar mM")
    for i in range(n):
        hours = years * 365 * 24
        print(f"  {i:4d} {b.kind[i]:5s} {2 * b.radius[i]:7.0f} {b.depth[i]:5.2f} "
              f"{100 * b.wet_hours[i] / hours:5.0f} {b.rewets[i] / years:12.1f} {b.ves_hours[i] / years:13.0f} "
              f"{b.events[i]:20d} {b.sweet_hours[i] / years:20.0f} {log['lipid_mM'][:, i].max():13.1f} "
              f"{log['sugar_mM'][:, i].max():13.2f}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--land-seed", type=int, default=3)
    p.add_argument("--seeds", type=int, nargs="+", default=[0], help="weather seeds")
    p.add_argument("--years", type=float, default=10)
    p.add_argument("--ponds", type=int, default=9)
    p.add_argument("--cvc", type=float, default=15.0, help="mM lipid needed for vesicles")
    p.add_argument("--vent-lipid", type=float, default=0.05, help="mM fatty acid in hot-spring water")
    p.add_argument("--dust", type=float, default=1e-9, help="mol lipid per m^2 per day from dust/meteorites")
    p.add_argument("--rain-hcho", type=float, default=0.01, help="mM formaldehyde in rain")
    p.add_argument("--formose", type=float, default=0.05, help="formose speed (per mM per day)")
    p.add_argument("--out", default="runs/e11")
    p.add_argument("--video", default=None, help="film the first seed (mp4); one frame per 2 days")
    args = p.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(args.out, exist_ok=True)
    for seed in args.seeds:
        b = Basin(seed=args.land_seed, ponds=args.ponds, cvc=args.cvc, vent_lipid=args.vent_lipid,
                  dust=args.dust, rain_hcho=args.rain_hcho, formose=args.formose)
        b.rng = np.random.default_rng(1000 + seed)   # the land stays, the weather changes
        movie = None
        if args.video and seed == args.seeds[0]:
            movie = film(b, args.video, plt)
        log = b.run(years=args.years, callback=movie)
        if movie:
            movie.close()
        out = os.path.join(args.out, f"basin{args.land_seed}_seed{seed}")
        np.savez_compressed(out + ".npz", **log, kind=b.kind, streams=b.streams, centers=b.centers,
                            radius=b.radius, rewets=b.rewets, events=b.events, ves_hours=b.ves_hours,
                            sweet_hours=b.sweet_hours, wet_hours=b.wet_hours)
        print(f"== land {args.land_seed} weather {seed}: streams {b.streams.tolist()}")
        summary(b, log, args.years)
        plot(b, log, out + ".png", plt)


def plot(b, log, path, plt):
    fig, axs = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw=dict(width_ratios=[1, 1.4]))
    draw_map(axs[0], b)
    axs[0].set_title("basin: dark ground = low; dotted = full pond; arrows = streams")
    days = log["hour"] / 24
    for i in range(len(b.radius)):
        axs[1].semilogy(days, np.maximum(log["lipid_mM"][:, i], 1e-4), lw=0.8,
                        color=plt.cm.tab10(i), label=f"{i} {b.kind[i]}")
    axs[1].axhline(b.cvc, color="k", ls="--", lw=1)
    axs[1].text(days[0], b.cvc * 1.2, "CVC: vesicles above this line", fontsize=8)
    axs[1].set_xlabel("day")
    axs[1].set_ylabel("lipid in water (mM)")
    axs[1].legend(fontsize=7, ncol=3)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


class film:
    """Callback for Basin.run: one frame every 2 days (map + lipid lines so far)."""

    def __init__(self, b, path, plt):
        from alife.viz import Movie

        class V:
            pass
        self.plt = plt
        self.fig, self.axs = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw=dict(width_ratios=[1.2, 1]))
        from matplotlib import cm
        import matplotlib.colors as mc
        self.fig.colorbar(cm.ScalarMappable(mc.LogNorm(1e-3, 10), "viridis"), ax=self.axs[0], shrink=0.6,
                          label="lipid / CVC  (1 and above = vesicles can form)")
        v = V()
        v.fig = self.fig
        self.movie = Movie(v, path, fps=20)
        self.b = b

    def __call__(self, b, log):
        if (b.hour - 1) % 48 >= 6:
            return
        for ax in self.axs:
            ax.clear()
        lip = log["lipid_mM"][-1]
        draw_map(self.axs[0], b, lip, log["flow"][-1], grow=2.0, crop=True)
        day = b.hour / 24
        season = "wet season" if b.season() > 0.5 else "dry season"
        self.axs[0].set_title(f"year {day // 365 + 1:.0f} day {day % 365:.0f}  {season}\n"
                              "ponds drawn 2x size; dotted = full, filled = water now; arrows = streams",
                              fontsize=10)
        d = np.array(log["hour"]) / 24
        L = np.array(log["lipid_mM"])
        recent = d > d[-1] - 365   # last year only, so it stays readable
        for i in range(len(b.radius)):
            self.axs[1].semilogy(d[recent], np.maximum(L[recent, i], 1e-4), lw=1, color=self.plt.cm.tab10(i),
                                 label=f"{i} {b.kind[i]}")
        self.axs[1].axhline(b.cvc, color="k", ls="--", lw=1)
        self.axs[1].set_ylim(1e-4, 300)
        self.axs[1].set_xlabel("day")
        self.axs[1].set_ylabel("lipid in water (mM); dashed = CVC")
        self.axs[1].set_title("last 365 days")
        self.axs[1].legend(fontsize=8, ncol=3, loc="upper left")
        self.movie.grab()

    def close(self):
        self.movie.close()


if __name__ == "__main__":
    main()
