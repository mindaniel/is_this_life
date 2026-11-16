import glob, re, numpy as np, pandas as pd
tags = {"calm": r"sealH_s\d", "heat": "heat20", "pump10": "pump30_g10", "pump5": "pump30_g5", "holdacid": "shock10000x-3",
        "holdalk": "shock10000x3x", "random": "srate", "rain": "_rain"}
for t, pat in tags.items():
    for f in sorted(glob.glob("runs/e19/*.csv")):
        if not re.search(pat, f): continue
        d = pd.read_csv(f)
        e = d.iloc[-1]
        lit = d[d.light == 1] if d.light.notna().any() else d
        print(f"{t:8s} s{f[-5]} end {int(e.step):6d} A_in {e.A_in:.2f} B_in {e.B_in:.2f} pk_min {d.pockets.min()} "
              f"Ach {int(e.A_chains)} Bch {int(e.B_chains)} qin A-B mean {np.mean(d.A_qin - d.B_qin):+.3f} lit {np.mean(lit.A_qin - lit.B_qin):+.3f} "
              f"left A/B {d.A_left.mean():.0f}/{d.B_left.mean():.0f} water {e.water} pumps {e.pumps} burns {e.burns} vol {e.volumes}")
