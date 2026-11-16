"""e19 + e20 plots. e19: inside charged share A - B (light pump vs no pump). e20: dark polymer beads, light absorbed."""
import glob, re, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
ink, mute = "#333333", "#888888"
def style(ax, title, ylab):
    ax.set_title(title, loc="left", fontsize=11, color=ink)
    ax.set_ylabel(ylab, color=ink); ax.set_xlabel("step", color=ink)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.grid(alpha=0.25, lw=0.6)

# e19
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
groups = {"no pump": r"sealH_s\d", "pump, push 5 kT": "pump30_g5", "pump, push 10 kT": "pump30_g10"}
for k, (lab, pat) in enumerate(groups.items()):
    for j, f in enumerate(sorted(f for f in glob.glob("runs/e19/*.csv") if re.search(pat, f))):
        d = pd.read_csv(f); d = d[d.step <= 60000]
        ax[0].plot(d.step, d.A_qin - d.B_qin, color=C[k], lw=1.5, alpha=0.85, label=lab if j == 0 else None)
        ax[1].plot(d.step, d.A_chains - d.B_chains, color=C[k], lw=1.5, alpha=0.85, label=lab if j == 0 else None)
for a in ax:
    for t0 in (10000, 30000, 50000):
        a.axvspan(t0, t0 + 10000, color="#999999", alpha=0.12, lw=0)
style(ax[0], "Inside of bag A vs B: charged heads (lower = more acid in A)", "A_qin - B_qin")
style(ax[1], "Fatty acids in wall: A - B", "A - B")
ax[0].legend(frameon=False, fontsize=9); ax[0].text(15000, ax[0].get_ylim()[1] * 0.9, "night", color=mute, ha="center", fontsize=8)
fig.suptitle("e19: hand-placed pigment pump in bag A (3 seeds each; grey = dark)", x=0.01, ha="left", color=ink)
fig.tight_layout(); fig.savefig("runs/e19/e19_pump.png", dpi=130)

# e20
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
kinds = {"dark (no light)": ("dry", lambda f: "glow" not in f), "light + heat": ("dry", lambda f: "_g0_" in f),
         "light + heat + acid": ("dry", lambda f: "_g5_" in f), "wet (dilute HCN), light + heat": ("wet", lambda f: True)}
for k, (lab, (pre, ok)) in enumerate(kinds.items()):
    for s in range(3):
        f = [f for f in glob.glob(f"runs/e20/{pre}{s}_then_*.csv") if ok(f)][0]
        d = pd.read_csv(f)
        ax[0].plot(d.step, d.dark, color=C[k], lw=1.5, alpha=0.85, label=lab if s == 0 else None)
        ax[1].plot(d.step, d.dark_wall, color=C[k], lw=1.2, alpha=0.7, label=lab if s == 0 else None)
style(ax[0], "Beads in dark HCN polymer (chains of 5+)", "beads")
style(ax[1], "Dark polymer beads touching a bag wall", "beads")
ax[0].legend(frameon=False, fontsize=9)
fig.suptitle("e20: pond-made pigment next to bags (3 seeds each; joined at step 20k)", x=0.01, ha="left", color=ink)
fig.tight_layout(); fig.savefig("runs/e20/e20_glow.png", dpi=130)
