"""e24 plot: wall stretch (area per lipid), holes, sugar and probes still inside, per inside-sugar level (3 seeds each)."""
import glob
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

levels = ["0.03", "0.1", "0.2", "0.4", "0.8", "1.6"]
# blues = wall holds (0.03-0.4), oranges = wall tears (0.8, 1.6); darker = more sugar inside
cols = ["#B8C2CC", "#8FB0CF", "#5B8DB8", "#2E6A9E", "#C0632B", "#7A2E12"]
fig, ax = plt.subplots(1, 4, figsize=(17, 4))
for o, c in zip(levels, cols):
    for k, f in enumerate(sorted(glob.glob(f"runs/e24/stretch_*_in{o}_p40_s*.csv"))):
        d = pd.read_csv(f)
        x = d.step / 1000
        ax[0].plot(x, d.area, color=c, lw=1.5, label=f"sugar inside {o} (outside 0.03)" if k == 0 else None)
        ax[1].plot(x, d.holes * 100, color=c, lw=1.5)
        ax[2].plot(x, 100 * d.sugar_in / d.sugar_in[0], color=c, lw=1.5)
        ax[3].plot(x, d.probes_in, color=c, lw=1.5)
ax[0].set(title="Wall stretch: area per lipid", xlabel="steps (thousands)", ylabel="area per lipid")
ax[1].set(title="Holes in the wall (% of directions)", xlabel="steps (thousands)")
ax[2].set(title="Sugar still inside (% of start)", xlabel="steps (thousands)", ylim=(0, 105))
ax[3].set(title="Probes still inside (of 40)", xlabel="steps (thousands)", ylim=(0, 41))
ax[0].legend(frameon=False, fontsize=8, loc="upper right")
for a in ax:
    a.spines[["top", "right"]].set_visible(False); a.grid(alpha=0.25)
fig.suptitle("e24: with up to 0.4 sugar inside the wall only stretches (~7%); at 0.8 and 1.6 it tears, sugar pours out, "
             "the bag shrinks and the wall closes again")
fig.tight_layout()
fig.savefig("runs/e24/e24_stretch.png", dpi=120)
