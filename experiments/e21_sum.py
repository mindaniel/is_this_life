"""e21 summary: tar in walls, and inner charge vs the no-light tar run (same seed, same start)."""
import glob, numpy as np, pandas as pd
def get(s, tag):
    pat = {"tar0.5": "tar5o0.5_glow", "tar0.8": "tar5o0.8_glow", "tar0.8dark": "tar5o0.8_s"}[tag]
    return pd.read_csv([f for f in glob.glob(f"runs/e21/dry{s}_then_*.csv") if pat in f][0])
for tag in ("tar0.5", "tar0.8", "tar0.8dark"):
    for s in range(3):
        d, c = get(s, tag), get(s, "tar0.8dark")
        e = d.iloc[-1]
        late = d.step > 50000
        print(f"{tag:10s} s{s} dark {int(e.dark)} longest {int(e.longest)} wall {d.dark_wall[late].mean():.0f} tail {d.dark_tail[late].mean():.0f} "
              f"(max {int(d.dark_tail.max())}) abs {e.absorbed} acid {e.glow_acid} A_in {e.A_in:.2f} B_in {e.B_in:.2f} pk_min {d.pockets.min()} "
              f"qin A/B {d.A_qin[late].mean():.3f}/{d.B_qin[late].mean():.3f} vs dark {c.A_qin[late].mean():.3f}/{c.B_qin[late].mean():.3f} "
              f"q {d.A_q[late].mean():.3f} ch {int(e.A_chains)}/{int(e.B_chains)} loose {int(d.loose[late].mean())} vol {e.volumes}")
