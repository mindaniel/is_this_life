"""Checks for World.save / load / merge.

1. Save, load, carry on: same positions as the unbroken run (bit for bit).
2. Merge two worlds, each with a vesicle sitting across its periodic wall: both come out
   whole (cut at the emptiest slab), and the joined world runs.

    python experiments/check_save.py
"""
import os, sys, tempfile
import numpy as np
sys.path.insert(0, "experiments"); sys.path.insert(0, ".")
from alife.chem.world import World
from alife.chem.analysis import cluster_labels, water_pockets3
from e15_split3d import mixed_vesicle, scatter_fast

path = os.path.join(tempfile.gettempdir(), "check_save.pkl")


def bag(seed, centre):
    w = World(box=40.0, kT=1.1, seed=seed, dim=3)
    w.rules.set_oily("stail", 0.85)
    mixed_vesicle(w, np.array(centre), 8.0, 0.5, 0.5)
    near = lambda p: np.linalg.norm((p - centre) - 40.0 * np.round((p - centre) / 40.0), axis=1) < 4.0
    scatter_fast(w, 100, near)   # sugar inside, so the bag holds water
    w.step(300, max_force=50); w.step(500)
    return w


def biggest(w):
    lab = cluster_labels(w)
    return np.sort(np.bincount(lab[w.mol[np.isin(w.type, (1, 13))]]))[::-1][:3]


# 1
w = bag(1, [20.0, 20, 20])
w.save(path, note="hi")
w.step(500)
w2, extra = World.load(path)
w2.step(500)
print("1. carry on after load: max |pos difference| =", np.abs(w.pos - w2.pos).max(), " extra:", extra)

# 2
a = bag(2, [1.0, 20, 20])     # across the x wall
b = bag(3, [39.0, 2, 20])     # across the x and y walls
print("2. before: tail beads per clump  a", biggest(a), " b", biggest(b))
m = World.merge(a, b, axis=0, seed=5)
print("   merged box", m.box(), "beads", m.n, "=", a.n + b.n, " bonds", len(m.bonds), "=", len(a.bonds) + len(b.bonds))
m.step(200, max_force=50); m.step(2000)
v = water_pockets3(m)[0]
print("   after 2200 steps: tail beads per clump", biggest(m), " sealed pockets (volumes)", v)
