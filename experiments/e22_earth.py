"""e22 Earth seeds (d 1, warm photon heat): tar pump per seed (day vs night) and pump vs how unevenly
the tar sits in the two bags. Reads runs/e22/e22w_summary.csv (run e22_sum.py --warm first)."""
import ast
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

t = pd.read_csv("runs/e22/e22w_summary.csv")
t = t[(t.d == 1) & (t.end >= 80000)].sort_values("s").reset_index(drop=True)
b = t.bags.map(ast.literal_eval)
t["gap"] = [abs(x[0] - x[1]) for x in b]          # tar beads in walls: rich bag minus other
r = np.corrcoef(t.gap, t.day)[0, 1]
print(t[["s", "bags", "gap", "day", "night"]].to_string(index=False))
print(f"day mean {t.day.mean():.3f} sd {t.day.std():.3f}  night mean {t.night.mean():.3f} sd {t.night.std():.3f}"
      f"  seeds day < -0.02: {(t.day < -0.02).sum()}/{len(t)}  r(gap, day) {r:.2f}")
INK, MUTED, DAY, NIGHT = "#1f1f1f", "#6b6b6b", "#2a6fdb", "#9a9a9a"
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
x = np.arange(len(t))
a = ax[0]
a.bar(x - 0.2, t.day, 0.38, color=DAY, label="day")
a.bar(x + 0.2, t.night, 0.38, color=NIGHT, label="night")
a.set_xticks(x, [f"seed {s}" for s in t.s], fontsize=8, rotation=45)
a.set_title("Tar-rich bag's inside minus the other bag's\n(charged share; below 0 = more acid)", loc="left", fontsize=10, color=INK)
a.legend(frameon=False, fontsize=9)
a = ax[1]
a.scatter(t.gap, t.day, s=44, color=DAY, edgecolor="white", lw=1.5, zorder=3, label="day")
a.scatter(t.gap, t.night, s=44, color=NIGHT, edgecolor="white", lw=1.5, zorder=3, label="night")
k = np.polyfit(t.gap, t.day, 1); g = np.linspace(0, t.gap.max(), 2)
a.plot(g, np.polyval(k, g), color=DAY, lw=1.5, ls="--", zorder=2)
a.set_xlabel("how much more tar sits in one bag's wall than the other's (beads)", color=MUTED)
a.set_title(f"More uneven tar = stronger day pump (r = {r:.2f})", loc="left", fontsize=10, color=INK)
a.legend(frameon=False, fontsize=9)
for a in ax:
    a.axhline(0, color=MUTED, lw=0.8)
    a.grid(axis="y", color="#e6e6e6", lw=0.8)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    a.tick_params(colors=MUTED, labelsize=9)
fig.suptitle(f"e22, Earth (d = 1), 9 seeds: the tar pump. Day mean {t.day.mean():.3f} ± {t.day.std():.3f}, "
             f"night {t.night.mean():.3f} ± {t.night.std():.3f}", fontsize=11, color=INK)
fig.tight_layout()
fig.savefig("runs/e22/e22w_earth.png", dpi=130)
print("saved")
