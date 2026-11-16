import glob, numpy as np, pandas as pd
for t in ("dark", "heat", "acid", "wet"):
    for s in range(3):
        pre = "wet" if t == "wet" else "dry"
        fs = [f for f in glob.glob(f"runs/e20/{pre}{s}_then_*.csv")
              if (t == "dark") == ("glow" not in f) and (t != "acid" or "_g5_" in f) and (t != "heat" or "_g0_" in f)]
        d = pd.read_csv(fs[0]); e = d.iloc[-1]; b = d.iloc[0]
        print(f"{t:5s} s{s} {int(b.step)}->{int(e.step)} hcn {int(b.hcn)}->{int(e.hcn)} dark {int(b.dark)}->{int(e.dark)} longest {int(e.longest)} "
              f"wall max {int(d.dark_wall.max())} mean {d.dark_wall.mean():.1f} abs {e.absorbed} acid {e.glow_acid} "
              f"A_in {e.A_in:.2f} B_in {e.B_in:.2f} pk_min {d.pockets.min()} qin A/B {d.A_qin.mean():.3f}/{d.B_qin.mean():.3f} "
              f"q A/B {d.A_q.mean():.3f}/{d.B_q.mean():.3f} left {d.A_left.mean():.0f}/{d.B_left.mean():.0f} ch {int(e.A_chains)}/{int(e.B_chains)}")
