r"""Plot e26 explicit-water runs: water through the wall, water in the oily middle, biggest peptide group, wall thickness.

    C:\Apps\Python314\python.exe experiments/e26_plot.py
"""
import glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

SETUPS = [("_n0_L7_kT", "control: no peptides", "#7f7f7f"),
          ("_n16_L7_kT", "16 peptides, start alone", "#1f77b4"),
          ("_n32_L7_kT", "32 peptides, start alone", "#2ca02c"),
          ("_n12_L7_b6_", "2 ready-made rings of 6 (shortcut)", "#d62728")]
fig, ax = plt.subplots(1, 4, figsize=(17, 4.2))
for key, label, col in SETUPS:
    for k, f in enumerate(sorted(glob.glob(f"runs/e26full/wat_*{key}*_s?.csv"))):
        d = pd.read_csv(f)
        try:   # 30k carry-on, sides measured from the wall's local middle (e26_water --start ... --tag _local)
            e = pd.read_csv(f.replace(".csv", "_local.csv"))
            e["cross"] += d.cross.iloc[-1]
            d = pd.concat([d, e], ignore_index=True)
        except FileNotFoundError:
            pass
        lab = label if k == 0 else None
        t = d.step / 1000
        ax[0].plot(t, d.cross, color=col, lw=1.5, label=lab)
        ax[1].plot(t, d.in_core, color=col, lw=1.2, label=lab)
        if "n0_" not in key:
            ax[2].plot(t, d.pep_cluster, color=col, lw=1.2, label=lab)
        ax[3].plot(t, d.thick, color=col, lw=1.2, label=lab)
titles = ["water beads gone through the wall (total)", "water beads in the oily middle",
          "biggest peptide group (touching)", "wall thickness (head to head)"]
for a, ti in zip(ax, titles):
    a.set_title(ti, fontsize=10)
    a.set_xlabel("thousand steps")
    a.grid(alpha=0.3)
ax[0].legend(fontsize=8, loc="upper left")
for a in ax:
    a.axvline(100, color="k", lw=0.8, ls=":")
ax[0].text(101, 5, "recount:\nlocal wall\nmiddle", fontsize=7)
fig.suptitle("e26: water as beads, flat wall (3 seeds each)", fontsize=11)
fig.tight_layout()
fig.savefig("runs/e26full/e26_water.png", dpi=110)
print("saved runs/e26full/e26_water.png")
