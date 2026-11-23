"""e24 plot: wall stretch (area per lipid) and probes still inside, per inside-sugar level (3 seeds each)."""
import glob, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

levels = ["0.03", "0.1", "0.2", "0.4"]
cols = ["#9AA5B1", "#5B8DB8", "#2E6A9E", "#173F66"]   # one hue, light -> dark = more sugar inside
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
for o, c in zip(levels, cols):
    for k, f in enumerate(sorted(glob.glob(f"runs/e24/stretch_*_in{o}_p40_s*.csv"))):
        d = pd.read_csv(f)
        lab = f"inside {o} ({float(o) / 0.03:.0f}x outside)" if k == 0 else None
        ax[0].plot(d.step / 1000, d.area, color=c, lw=1.5, label=lab)
        ax[1].plot(d.step / 1000, d.holes * 100, color=c, lw=1.5)
        ax[2].plot(d.step / 1000, d.probes_in, color=c, lw=1.5)
ax[0].set(title="Wall stretch: area per lipid", xlabel="steps (thousands)", ylabel="area per lipid")
ax[1].set(title="Holes in the wall (% of directions)", xlabel="steps (thousands)", ylim=(-0.5, 5))
ax[2].set(title="Probes still inside (of 40)", xlabel="steps (thousands)", ylim=(30, 41))
ax[0].legend(frameon=False, fontsize=8, loc="upper center", ncol=2)
ax[0].set_ylim(1.34, 1.56)
for a in ax:
    a.spines[["top", "right"]].set_visible(False); a.grid(alpha=0.25)
fig.suptitle("e24: more sugar inside stretches the wall up to ~7%, but it never opens (1 probe of 480 got out)")
fig.tight_layout()
fig.savefig("runs/e24/e24_stretch.png", dpi=120)
