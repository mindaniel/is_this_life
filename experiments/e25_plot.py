r"""Plot e25 peptide runs: peptides in the wall (upright), biggest peptide group, probes and sugar still inside.

    C:\Apps\Python314\python.exe experiments/e25_plot.py
"""
import glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

SETUPS = [("n40_wall", "peptides start in the wall", "#1f77b4"),
          ("n40_water", "peptides start in the water", "#ff7f0e"),
          ("n0_wall", "control: no peptides", "#7f7f7f")]
fig, ax = plt.subplots(1, 4, figsize=(17, 4.2))
for key, label, col in SETUPS:
    for k, f in enumerate(sorted(glob.glob(f"runs/e25/pep_*_{key}_*_s?.csv"))):
        d = pd.read_csv(f)
        lab = label if k == 0 else None
        t = d.step / 1000
        if "n40" in key:
            ax[0].plot(t, d.span, color=col, lw=1.5, label=lab)
            ax[1].plot(t, d.pep_cluster, color=col, lw=1.5, label=lab)
        ax[2].plot(t, d.probes_in, color=col, lw=1.5, label=lab)
        ax[3].plot(t, 100 * d.sugar_in / d.sugar_in[0], color=col, lw=1.5, label=lab)
titles = ["peptides standing across the wall (of 40)", "biggest peptide group (touching)",
          "probes still inside (of 40)", "sugar still inside (% of start)"]
for a, ti in zip(ax, titles):
    a.set_title(ti, fontsize=10)
    a.set_xlabel("thousand steps")
    a.grid(alpha=0.3)
ax[0].set_ylim(-1, 41); ax[1].set_ylim(0, 10); ax[2].set_ylim(0, 41); ax[3].set_ylim(0, 105)
ax[2].legend(fontsize=8, loc="lower left")
fig.suptitle("e25: do peptides gather into pores? (3 seeds each)  No: they stay one by one, nothing leaks", fontsize=11)
fig.tight_layout()
fig.savefig("runs/e25/e25_peptide.png", dpi=110)
print("saved runs/e25/e25_peptide.png")
