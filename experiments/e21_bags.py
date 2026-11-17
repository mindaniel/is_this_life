"""e21: tar beads in each bag's wall at the end (touching that wall's tails), and the light/dark swing."""
import glob, sys, numpy as np, pandas as pd
sys.path.insert(0, "experiments"); sys.path.insert(0, ".")
from scipy.spatial import cKDTree
from alife.chem.world import World, TAILS
from alife.chem.analysis import cluster_labels
from alife.chem.nitrogen import dark_beads
from e16_catalyst import wall_of
for f in sorted(glob.glob("runs/e21/*_state.pkl")):
    w, _ = World.load(f)
    lab = cluster_labels(w)
    lip = np.flatnonzero(np.isin(w.type, TAILS))
    p = w._wrapped(); dark = dark_beads(w, 5)
    out = []
    for m in w.markers:
        c = wall_of(w, m, lab, lip)
        t = lip[lab[w.mol[lip]] == c]
        out.append(int((cKDTree(p[t], boxsize=w.L).query(p[dark], distance_upper_bound=1.5)[0] < np.inf).sum()))
    d = pd.read_csv(f.replace("_state.pkl", ".csv"))
    late = d[d.step > 40000]
    sw = (late.A_qin - late.B_qin)
    lit = sw[late.light == 1].mean() if "glow" in f else np.nan
    print(f.split("sealH_")[1][:-10], "tar in wall A/B", out, f"A-B qin lit {lit:+.3f} night {sw[late.light != 1].mean() if 'glow' in f else sw.mean():+.3f}")
