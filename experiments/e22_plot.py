"""e22 plot: per planet distance, tar made, tar sunk into walls, loose fatty acids, and the tar pump
(tar-rich bag minus the other, inner charged share, day minus night). Reads runs/e22/e22_summary.csv."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import sys
VAR = "e22w" if "--warm" in sys.argv else "e22"
t = pd.read_csv(f"runs/e22/{VAR}_summary.csv")
old = pd.read_csv("runs/e22/e22_summary.csv") if VAR == "e22w" else None   # kicked heat, for comparison
if old is not None:
    old["pump"] = old.day - old.night
t["pump"] = t.day - t.night
ok = t.end >= 80000
INK, MUTED, C, BAD = "#1f1f1f", "#6b6b6b", "#2a6fdb", "#c2410c"
fig, ax = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
panels = [("dark", "Tar made (dark beads at end)"), ("tail", "Tar sunk into walls (beads touching tails)"),
          ("loose", "Loose fatty acids (mean)"), ("pump", "Tar pump: day minus night\n(tar-rich bag's inside, charged share)")]
for a, (col, title) in zip(ax.flat, panels):
    if old is not None:
        o = old[old.end >= 80000].groupby("d")[col].mean()
        a.plot(o.index, o.values, color="#b0b0b0", lw=1.5, ls="--", zorder=1, label="heat as kick (first run)")
    m = t[ok].groupby("d")[col].mean()
    a.plot(m.index, m.values, color=C, lw=2, zorder=2)
    a.scatter(t.d[ok], t[col][ok], s=36, color=C, edgecolor="white", lw=1.5, zorder=3)
    a.set_xlim(0.62, 1.58)
    a.set_title(title, fontsize=10, color=INK, loc="left")
    a.grid(axis="y", color="#e6e6e6", lw=0.8)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    a.tick_params(colors=MUTED, labelsize=9)
ax[1, 1].axhline(0, color=MUTED, lw=0.8)
if old is not None:
    ax[0, 0].legend(frameon=False, fontsize=8)
for a in ax[1]:
    a.set_xlabel("distance from sun (Earth = 1)", color=MUTED)
    a.set_xticks([0.7, 0.85, 1.0, 1.2, 1.5], ["0.7\nhot, bright", "0.85", "1\nEarth", "1.2", "1.5\ncold, dim"])
crash = ", ".join(f"{n} at {d:g}" for d, n in (~ok).groupby(t.d).sum().items() if n)
fig.suptitle(("e22, heat as warm spot" if VAR == "e22w" else "e22") + ": the same pond on planets closer to and further from the sun (dots = finished runs, line = mean)"
             + (f"\ncrashed before the end (too hot at noon): {crash}" if crash else ""), fontsize=11, color=INK)
fig.tight_layout()
fig.savefig(f"runs/e22/{VAR}_zone.png", dpi=130)
print("saved")
