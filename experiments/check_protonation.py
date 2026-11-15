"""Checks for reactions.Protonation: kernel dE == engine energy change (fast and numpy), and lone
fatty acids titrate like the textbook (charged = 1 / (1 + 10^(pKa - pH))).

    python experiments/check_protonation.py
"""
import numpy as np, sys
sys.path.insert(0, "experiments"); sys.path.insert(0, ".")
from alife.chem.world import World
from alife.chem.reactions import Protonation, FA_HEAD, FA_ACID
from alife.chem.elements import INDEX
from e15_split3d import mixed_vesicle

# 1. dE from the kernel == change in the engine's total energy (fast and numpy)
w = World(box=30.0, kT=1.1, seed=1, dim=3)
w.rules.set_oily("stail", 0.85)
mixed_vesicle(w, np.array([15.0, 15, 15]), 8.0, 0.5, 0.5, acid="fa_acid")
w.step(300, max_force=50); w.step(200)
pr = Protonation(ph=7.0, pka=7.0)   # mu = 0: dE is pure energy
dE = pr.run(w, apply=False)
heads = np.flatnonzero(np.isin(w.type, (FA_HEAD, FA_ACID)))
print("heads", len(heads), "charged", pr.charged(w))
worst = 0
for fastmode in (True, False):
    w.fast = fastmode
    for h in w.rng.choice(len(heads), 15, replace=False):
        i = heads[h]
        w._pairs = None; e0 = w.forces()[1]
        old = w.type[i]; w.type[i] = FA_ACID if old == FA_HEAD else FA_HEAD
        w._pairs = None; e1 = w.forces()[1]
        w.type[i] = old
        worst = max(worst, abs((e1 - e0) - dE[h]))
    print("fast" if fastmode else "numpy", "max |dE_kernel - dE_engine| =", worst)
w.fast = True

# 2. lone fatty acids far apart: fraction charged = 1 / (1 + 10^(pKa - pH))
from alife.chem.molecules import add_lipid3
for ph in (3.8, 4.8, 5.8):
    w = World(box=60.0, kT=1.1, seed=2, dim=3)
    g = np.arange(5) * 12 + 6
    for x in g:
        for y in g:
            for z in g:
                add_lipid3(w, np.array([x, y, z], float), np.array([1.0, 0, 0]), 2, head="fa_acid", tail="stail")
    pr = Protonation(ph=ph, pka=4.8)
    f = []
    for k in range(400):
        pr.run(w)
        if k >= 50: f.append(pr.charged(w))
    print(f"pH {ph}: charged {np.mean(f):.3f}  ideal {1 / (1 + 10 ** (4.8 - ph)):.3f}")
