"""e21 plot: tar sunk into walls over time; inside charge of the tar-rich bag minus the other (light vs dark)."""
import glob, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
C = ["#2a78d6", "#eb6834", "#1baf7a"]
ink = "#333333"
rich = {0: "A", 1: "A", 2: "B"}   # bag with more tar in its wall at the end (e21_bags.py), same in all 3 versions
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
for k, (lab, pat) in enumerate({"tar 0.5, light": "tar5o0.5_glow", "tar 0.8, light": "tar5o0.8_glow",
                                 "tar 0.8, no light": "tar5o0.8_s"}.items()):
    for s in range(3):
        d = pd.read_csv([f for f in glob.glob(f"runs/e21/dry{s}_then_*.csv") if pat in f][0])
        sign = 1 if rich[s] == "A" else -1
        ax[0].plot(d.step, d.dark_tail, color=C[k], lw=1.3, alpha=0.8, label=lab if s == 0 else None)
        ax[1].plot(d.step, sign * (d.A_qin - d.B_qin), color=C[k], lw=1.5, alpha=0.85, label=lab if s == 0 else None)
for a, (t, y) in zip(ax, [("Tar beads sunk into lipid walls (touching tails)", "beads"),
                          ("Inside charge: tar-rich bag minus the other (lower = more acid)", "difference")]):
    a.set_title(t, loc="left", fontsize=11, color=ink); a.set_ylabel(y); a.set_xlabel("step")
    for sp in ("top", "right"): a.spines[sp].set_visible(False)
    a.grid(alpha=0.25, lw=0.6)
for t0 in (30000, 50000, 70000):
    ax[1].axvspan(t0, t0 + 10000, color="#999999", alpha=0.12, lw=0)
ax[1].axhline(0, color="#888888", lw=0.8)
ax[0].legend(frameon=False, fontsize=9)
fig.suptitle("e21: HCN tar made by the pond sinks into bag walls; in light, the bag with more tar goes acid (grey = night)",
             x=0.01, ha="left", color=ink, fontsize=11)
fig.tight_layout(); fig.savefig("runs/e21/e21_tar.png", dpi=130)
