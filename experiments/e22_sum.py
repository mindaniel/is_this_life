"""e22 summary per planet distance: bags sealed, tar made and in walls, and the tar pump (tar-rich bag minus
the other's inner charged share, day vs night). Per-bag tar from the final state (as e21_bags); crashed
runs have no final state and are summarised up to the crash."""
import glob, os, sys, numpy as np, pandas as pd
sys.path.insert(0, "experiments"); sys.path.insert(0, ".")
from scipy.spatial import cKDTree
from alife.chem.world import World, TAILS
from alife.chem.analysis import cluster_labels
from alife.chem.nitrogen import dark_beads
from e16_catalyst import wall_of


def tar_per_bag(f):
    w, _ = World.load(f)
    lab = cluster_labels(w)
    lip = np.flatnonzero(np.isin(w.type, TAILS))
    p, dark = w._wrapped(), dark_beads(w, 5)
    out = []
    for m in w.markers:
        t = lip[lab[w.mol[lip]] == wall_of(w, m, lab, lip)]
        out.append(int((cKDTree(p[t], boxsize=w.L).query(p[dark], distance_upper_bound=1.5)[0] < np.inf).sum()))
    return out


rows = []
for f in sorted(glob.glob("runs/e22/sun*_s?.csv")):
    d = pd.read_csv(f)
    dist, seed = float(f.split("_sun")[1].split("_")[0]), int(f[-5])
    st = f.replace(".csv", "_state.pkl")
    bags = tar_per_bag(st) if os.path.exists(st) else None
    rich = "A" if bags is None or bags[0] >= bags[1] else "B"
    sw = (d.A_qin - d.B_qin) * (1 if rich == "A" else -1)
    late = d.step > 40000
    e = d.iloc[-1]
    rows.append(dict(d=dist, s=seed, end=int(e.step), kT=round(d.kT.mean(), 3), sealed=int(d.pockets.min()) >= 2,
                     A_in=round(d.A_in.min(), 2), B_in=round(d.B_in.min(), 2), loose=int(d.loose.mean()),
                     fa=int(e.A_fa + e.B_fa), hcn0=int(d.hcn.iloc[0]), hcn=int(e.hcn), dark=int(e.dark), longest=int(e.longest),
                     tail=int(d.dark_tail[late].mean()) if late.any() else int(d.dark_tail.mean()), bags=bags,
                     day=round(sw[late & (d.light == 1)].mean(), 3) if late.any() else np.nan,
                     night=round(sw[late & (d.light != 1)].mean(), 3) if late.any() else np.nan,
                     dip=round(sw[d.light == 1].min(), 3), acid=int(e.glow_acid)))
t = pd.DataFrame(rows)
pd.set_option("display.width", 250)
print(t.to_string(index=False))
t.to_csv("runs/e22/e22_summary.csv", index=False)
