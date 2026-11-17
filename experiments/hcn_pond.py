"""A small pond where HCN is concentrated (a drying puddle, or ice where HCN gathers between
crystals): HCN + ammonia make the dark HCN polymer here. Saved with World.save, to be joined to
a bag world with merge_worlds.py (two waters meeting).

--like state.pkl copies the other world's box cross-section, temperature and rules, so the two
fit together at the seam. Sugar osmolyte at --osm (the bag world's outside level), so joining
doesn't dilute the bags' surroundings.

Shortcut (flagged): HCN starts concentrated (2000 per 30^3 ~ 0.07 per bead volume; real ponds
reach that only by drying or freezing); there is no drying step here.

    python experiments/hcn_pond.py --like runs/e20/bags_state.pkl --x 24 --hcn 4000 --nh3 1300 \
        --steps 20000 --out runs/e20/pond.pkl
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from alife.chem.nitrogen import HCN, NH3, Nitrogen, census, dark_beads   # noqa: E402
from alife.chem.world import World   # noqa: E402
from e15_split3d import scatter_fast   # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--like", default=None, help="saved world to match (cross-section, kT, rules)")
    p.add_argument("--x", type=float, default=24.0, help="pond length along x (the seam axis)")
    p.add_argument("--box", type=float, default=48.0, help="y and z size (ignored with --like)")
    p.add_argument("--kT", type=float, default=1.1, help="temperature (ignored with --like)")
    p.add_argument("--kT0", type=float, default=None, help="temperature the reaction rates are calibrated at "
                                                          "(default: this pond's; set it to let heat speed them)")
    p.add_argument("--hcn", type=int, default=4000)
    p.add_argument("--nh3", type=int, default=1300)
    p.add_argument("--osm", type=float, default=0.15, help="sugar osmolyte per unit volume")
    p.add_argument("--dark-len", type=int, default=5)
    p.add_argument("--steps", type=int, default=20000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    if args.like:
        other, _ = World.load(args.like)
        box, kT, rules = (args.x, *other.box()[1:]), other.kT, other.rules
    else:
        box, kT, rules = (args.x, args.box, args.box), args.kT, None
    w = World(box=box, rules=rules, kT=kT, seed=args.seed, dim=3)
    w.repel_inside = (HCN, NH3)   # floppy chains: their beads must not overlap
    anywhere = lambda q: np.ones(len(q), bool)
    scatter_fast(w, int(round(args.osm * np.prod(box))), anywhere)
    scatter_fast(w, args.hcn, anywhere, name="hcn")
    scatter_fast(w, args.nh3, anywhere, name="nh3")
    w.step(200, max_force=50)
    w.t = 0
    nit = Nitrogen(every=50, kT0=args.kT0 or w.kT)   # rates as calibrated, at kT0
    t0 = time.time()
    while w.t < args.steps:
        w.step(50)
        nit.react(w)
        if w.t % 2000 == 0:
            c = census(w)
            h = c["hcn_chains"]
            print(f"step {w.t} lone HCN {c['hcn']} NH3 {c['nh3']} dark beads {len(dark_beads(w, args.dark_len))} "
                  f"longest {int(np.flatnonzero(h)[-1]) if h[2:].any() else 1} lost {nit.hydrolysed} "
                  f"({time.time() - t0:.0f}s)", flush=True)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    w.save(args.out, args=vars(args))


if __name__ == "__main__":
    main()
