"""e24 pore picture: wall seen from outside + slice through the biggest hole, at 3 times of one run.
Needs states at step 0 and 5000 (rerun e24_stretch.py --osm-in 1.6 --seed 0 --steps 0 / 5000 --out DIR/t0, DIR/t5k).
    python experiments/e24_pore.py DIR"""
import sys, numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "experiments")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from alife.chem.world import World, TAILS
from e15_split3d import biggest, offset
S = sys.argv[1]
name = "stretch_R10_mix0.5_out0.03_in1.6_p40_s0_state.pkl"
paths = [(f"{S}/t0/{name}", "step 0: just filled (already torn while filling)"), (f"{S}/t5k/{name}", "step 5,000: torn"),
         (f"runs/e24/{name}", "step 60,000: sugar gone, healed")]
fig, ax = plt.subplots(2, 3, figsize=(15, 10))
for col, (p, title) in enumerate(paths):
    w, ex = World.load(p)
    mask, c, _ = biggest(w)
    x = offset(w, w._wrapped(), c)
    tail = mask & np.isin(w.type, TAILS); head = mask & ~np.isin(w.type, TAILS)
    sugar = np.zeros(w.n, bool); sugar[ex["sugar_in"]] = True
    probe = np.zeros(w.n, bool); probe[ex["probes"]] = True
    R = np.median(np.linalg.norm(x[tail], axis=1))
    # biggest hole: direction whose point on the wall is farthest from any tail bead
    k = np.arange(4000) + 0.5; z = 1 - 2 * k / 4000; ph = np.pi * (1 + 5 ** 0.5) * k
    u = np.stack([np.sqrt(1 - z**2) * np.cos(ph), np.sqrt(1 - z**2) * np.sin(ph), z], 1)
    gap = cKDTree(x[tail]).query(u * R)[0]
    h = u[gap.argmax()]
    n = np.cross(h, [0, 0, 1.0] if abs(h[2]) < 0.9 else [1.0, 0, 0]); n /= np.linalg.norm(n)
    e2 = np.cross(n, h)
    # top: the wall seen from outside, looking straight at the biggest hole (front half only)
    a = ax[0, col]
    front = (x @ h) > 0
    for m, cl, s in ((tail & front, "#E0A030", 6), (head & front, "#3C78B4", 6)):
        a.scatter(x[m] @ n, x[m] @ e2, s=s, c=cl, lw=0)
    a.set_title(f"{title}\nwall seen from outside (front half)  R {R:.1f}")
    # bottom: thin slice through the centre and the hole
    a = ax[1, col]
    sl = np.abs(x @ n) < 1.2
    for m, cl, s, lab in ((tail & sl, "#E0A030", 14, "lipid tails (oily)"), (head & sl, "#3C78B4", 14, "lipid heads"),
                          (sugar & sl, "#9AA0A6", 8, "sugar that started inside"), (probe & sl, "#D62728", 30, "probes")):
        a.scatter(x[m] @ h, x[m] @ e2, s=s, c=cl, lw=0, label=lab)
    a.annotate("biggest hole", (R * 1.25, 0), color="k", fontsize=9, va="center")
    a.set_title("slice through the middle (2 beads thick)")
    if col == 0: a.legend(loc="lower left", fontsize=8, frameon=False, markerscale=1.5)
    print(title, "R", round(R, 2), "biggest gap", round(gap.max(), 2), "sugar inside", int((np.linalg.norm(x[sugar], axis=1) < R).sum()))
for a in ax.ravel():
    a.set_aspect("equal"); a.set_xlim(-19, 19); a.set_ylim(-19, 19); a.set_xticks([]); a.set_yticks([])
    for s in a.spines.values(): s.set_visible(False)
fig.suptitle("e24, sugar inside 1.6 (seed 0): the stretched wall tears, sugar streams out, the wall heals", fontsize=13)
fig.tight_layout(); fig.savefig("runs/e24/e24_pore.png", dpi=110)
