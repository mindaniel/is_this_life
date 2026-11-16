"""Checks for reactions.Protonation: kernel dE == engine energy change (fast and numpy), and lone
fatty acids titrate like the textbook (charged = 1 / (1 + 10^(pKa - pH))). Also: an H+ hand-over
between two heads (dE == engine change), sealed-inside heads (`facing`), and Burn keeping the H+ count.

    python experiments/check_protonation.py
"""
import numpy as np, sys
sys.path.insert(0, "experiments"); sys.path.insert(0, ".")
from alife.chem import fast
from alife.chem.world import World
from alife.chem.reactions import Burn, Protonation, FA_HEAD, FA_ACID, FUEL, WASTE
from alife.chem.elements import INDEX
from e15_split3d import mixed_vesicle, scatter_fast

# 1. dE from the kernel == change in the engine's total energy (fast and numpy)
w = World(box=30.0, kT=1.1, seed=1, dim=3)
w.rules.set_oily("stail", 0.85)
mixed_vesicle(w, np.array([15.0, 15, 15]), 8.0, 0.5, 0.5, acid="fa_acid", alcohol=0.2)
w.step(300, max_force=50); w.step(200)
pr = Protonation(ph=7.0, pka=7.0)   # mu = 0: dE is pure energy
dE = pr.run(w, apply=False)
heads = np.flatnonzero(np.isin(w.type, (FA_HEAD, FA_ACID)))
print("heads", len(heads), "charged", pr.charged(w), "alcohols", int((w.type == INDEX["fa_oh"]).sum()))
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

# 3. H+ hand-over: flip a COOH and a nearby COO- together; kernel energy == engine change
heads, (sp, sn), (cpt, cn) = pr.lists(w)
r = w.rules
where = np.full(w.n, -1); where[heads] = np.arange(len(heads))
worst, tried = 0, 0
for a in range(len(heads)):
    if w.type[heads[a]] != FA_ACID:
        continue
    for j in sn[sp[a]:sp[a + 1]]:
        b = where[j]
        if b < 0 or w.type[j] != FA_HEAD:
            continue
        typ = w.type.copy()
        args = (sp, sn, cpt, cn)
        d = fast._flip_energy(w.pos, w.box(), a, heads[a], FA_HEAD, *args, typ, r.sigma, r.rc, r.eps, r.wc,
                              r.charge, r.coulomb, r.debye, r.charge_cutoff, 0.0)
        typ[heads[a]] = FA_HEAD
        d = fast._flip_energy(w.pos, w.box(), b, j, FA_ACID, *args, typ, r.sigma, r.rc, r.eps, r.wc,
                              r.charge, r.coulomb, r.debye, r.charge_cutoff, d)
        w._pairs = None; e0 = w.forces()[1]
        w.type[heads[a]], w.type[j] = FA_HEAD, FA_ACID
        w._pairs = None; e1 = w.forces()[1]
        w.type[heads[a]], w.type[j] = FA_ACID, FA_HEAD
        worst = max(worst, abs((e1 - e0) - d))
        tried += 1
        break
    if tried >= 15:
        break
print(f"hand-over: {tried} pairs, max |dE_kernel - dE_engine| = {worst:.2e}")
w._pairs = None
n_h = (w.type == FA_ACID).sum()
pr2 = Protonation(ph=5.0, hand_over=True)
pr2.closed = heads                     # all sealed: only hand-overs, so the H+ count can't change
for _ in range(20):
    pr2.run(w)
print("all sealed, 20 passes: H+ on heads", n_h, "->", (w.type == FA_ACID).sum(), "hops", pr2.hops)

# 4. heads facing the sealed inside: about the inner leaflet's share (inner area / total)
inner = Protonation.facing(w)
rin, rout = 8.0 - 1.5, 8.0 + 1.5
print(f"facing inside: {len(inner)}/{len(heads)} = {len(inner) / len(heads):.2f} "
      f"(inner leaflet share ~{rin ** 2 / (rin ** 2 + rout ** 2):.2f})")

# 5. Burn: each burnt fuel puts exactly one H+ on a head
w.type[heads] = FA_HEAD
scatter_fast(w, 30, lambda p: np.linalg.norm(p - 15.0, axis=1) < 4.0, name="fuel")
scatter_fast(w, 10, lambda p: np.linalg.norm(p - 15.0, axis=1) < 4.0, name="cat")
w._pairs = None
bu = Burn(p_burn=0.5)
for _ in range(40):
    w.step(50); bu.react(w)
print("burns", bu.burns, "waste", (w.type == WASTE).sum(), "COOH", (w.type == FA_ACID).sum(),
      "fuel left", (w.type == FUEL).sum())

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
