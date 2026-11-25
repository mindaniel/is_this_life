r"""Side view of an e26 wall: thin slab (in y) through the spot with the most water in the oily middle.

    C:\Apps\Python314\python.exe experiments/e26_look.py runs/e26full/<run>_state.pkl out.png
"""
import sys
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))
from alife.chem.elements import ELEMENTS, INDEX   # noqa: E402
from alife.chem.world import World   # noqa: E402
from e26_water import mid_z   # noqa: E402

w, ex = World.load(sys.argv[1])
box = w.box()
zm = mid_z(w)
p = w._wrapped()
dz = (p[:, 2] - zm + box[2] / 2) % box[2] - box[2] / 2
core = (w.type == INDEX["wat"]) & (np.abs(dz) < 1.5)
# densest y band of core water
ys = np.linspace(0, box[1], 80, endpoint=False)
cnt = [np.sum(np.abs((p[core, 1] - y + box[1] / 2) % box[1] - box[1] / 2) < 1.5) for y in ys]
y0 = ys[int(np.argmax(cnt))]
band = np.abs((p[:, 1] - y0 + box[1] / 2) % box[1] - box[1] / 2) < 1.5
color = np.array([e[4] for e in ELEMENTS], dtype=object)
color[INDEX["pep"]] = "#2CA02C"   # peptide water-loving stripe: green (pep blue looks like water)
color[INDEX["pho"]] = "#5B2C06"
fig, ax = plt.subplots(figsize=(8, 4.6))
for show, size, alpha in ((w.type == INDEX["wat"], 18, 0.35), (w.type != INDEX["wat"], 30, 1.0)):
    m = band & show & (np.abs(dz) < 7)
    ax.scatter(p[m, 0], dz[m], s=size, c=color[w.type[m]], alpha=alpha, edgecolors="none")
ax.set_aspect("equal")
ax.set_xlabel("x (bead sizes)"); ax.set_ylabel("height from wall middle")
ax.set_title("side slice (3 beads thick): water blue, heads/tails as usual, peptide green (water-loving) / dark brown (oily)", fontsize=9)
fig.tight_layout()
fig.savefig(sys.argv[2], dpi=120)
print("saved", sys.argv[2], "y0", y0)
