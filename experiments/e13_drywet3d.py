"""e13: wet-dry cycles in 3D, the way a real pond dries (implicit water).

Drying = water leaves = the box shrinks slowly, squeezing the lipids into
stacked sheets (as on a dried pond rim). Rewetting = the box grows back and
the sheets can peel off and close into vesicles (Deamer's mechanism).
Each cycle: wet -> dry (shrink) -> stay dry -> rewet (grow) -> wet.

The lipid amount is given as a real concentration (mM) with bead size
--sigma nm, so the e11 landscape model can say which numbers to use: the
critical vesicle concentration of short fatty acids is ~10-40 mM.

    python experiments/e13_drywet3d.py --mM 15 --cycles 3
    python experiments/e13_drywet3d.py --mM 15 --cycles 0      # control: always wet, same length
    python experiments/e13_drywet3d.py --mM 15 --cycles 1 --rewet-box 120 --last-wet 500000   # flood
    python experiments/e13_drywet3d.py --mM 15 --cycles 1 --rewet-box 120 --no-dry   # control: flood, never dried
"""
import argparse
import csv
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from alife.chem.analysis import cluster_labels, tail_clusters, water_pockets3   # noqa: E402
from alife.chem.molecules import gas3   # noqa: E402
from alife.chem.world import World   # noqa: E402

PER_NM3_PER_MM = 6.022e-4   # molecules per nm^3 in 1 mM


def schedule(args):
    """(phase, steps) list for the whole run."""
    out = [("wet", args.wet)]
    for _ in range(args.cycles):
        out += [("drying", args.dry), ("dry", args.hold), ("rewetting", args.dry), ("wet", args.wet)]
    if args.last_wet:
        out[-1] = ("wet", args.last_wet)
    if args.cycles == 0:   # control: same total length, always wet
        out = [("wet", args.wet + 3 * (2 * args.dry + args.hold + args.wet))]
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mM", type=float, default=15.0, help="lipid concentration when wet")
    p.add_argument("--sigma", type=float, default=0.9, help="bead size in nm")
    p.add_argument("--box", type=float, default=60.0, help="wet box side (beads)")
    p.add_argument("--packing", type=float, default=0.5, help="dry: beads fill this fraction of the box")
    p.add_argument("--cycles", type=int, default=3)
    p.add_argument("--wet", type=int, default=20000)
    p.add_argument("--dry", type=int, default=30000, help="steps to dry (and to rewet)")
    p.add_argument("--hold", type=int, default=10000, help="steps staying dry")
    p.add_argument("--rewet-box", type=float, default=None,
                   help="box side after rewetting (more water than at the start = a flood)")
    p.add_argument("--last-wet", type=int, default=None, help="steps of the final wet phase")
    p.add_argument("--no-dry", action="store_true",
                   help="control: same schedule, but the box doesn't shrink (no drying), only the flood")
    p.add_argument("--kT", type=float, default=1.1)
    p.add_argument("--every", type=int, default=1000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", default="runs/e13")
    p.add_argument("--video", default=None)
    p.add_argument("--fps", type=int, default=15, help="video frames per second (one frame per --every steps)")
    args = p.parse_args()

    n = int(round(args.mM * PER_NM3_PER_MM * (args.box * args.sigma) ** 3))
    w = World(box=args.box, kT=args.kT, seed=args.seed, dim=3)
    gas3(w, n)
    w.step(300, max_force=50)
    L_wet = args.rewet_box or w.L
    L_dry = w.L if args.no_dry else (w.n / args.packing) ** (1 / 3)
    name = f"mM{args.mM:g}_c{args.cycles}{'_flood' + format(L_wet, 'g') if args.rewet_box else ''}{'_nodry' if args.no_dry else ''}_s{args.seed}"
    os.makedirs(args.out, exist_ok=True)
    print(f"{name}: {n} lipids ({w.n} beads), box {w.L:.0f} wet -> {L_dry:.1f} dry "
          f"(water volume down {(w.L / L_dry) ** 3:.0f}x), rewet to {L_wet:.0f}", flush=True)
    film = None
    if args.video:
        from e12_3d import Film
        film = Film(w, args.video, fps=args.fps)
    fh, out, t0 = open(os.path.join(args.out, name + ".csv"), "w", newline=""), None, time.time()
    snap = lambda tag: np.savez_compressed(os.path.join(args.out, f"{name}_{tag}.npz"),
                                           pos=w.pos, type=w.type, mol=w.mol, L=w.L)
    had_pocket, chunk = False, 100
    for phase, steps in schedule(args):
        target = {"drying": L_dry, "rewetting": L_wet}.get(phase)
        f = (target / w.L) ** (chunk / steps) if target else 1.0
        for k in range(0, steps, chunk):
            flood = phase == "rewetting" and args.rewet_box
            if f != 1.0:   # a flood moves whole clumps (water flows around membranes)
                w.rescale(f, cluster_labels(w) if flood else None)
            if flood:   # clumps moving apart can push a stray lipid into a membrane: let it ease out
                w.step(10, max_force=50)
                w.step(chunk - 10, max_force=None)
            else:
                w.step(chunk, max_force=None)
            if w.t % args.every < chunk:
                cl = tail_clusters(w)
                vol, _ = water_pockets3(w)
                row = dict(step=w.t, phase=phase, box=round(w.L, 2), clusters=len(cl), largest=int(cl[0]),
                           pockets=len(vol), volume=int(vol.sum()), biggest=int(vol[0]) if len(vol) else 0)
                if out is None:   # written as we go: a killed run keeps its log
                    out = csv.DictWriter(fh, fieldnames=list(row))
                    out.writeheader()
                out.writerow(row)
                fh.flush()
                if len(vol) and (not had_pocket or w.t % (10 * args.every) < chunk):
                    snap(f"pocket_{w.t}")   # so a detection can be checked by eye later
                had_pocket = len(vol) > 0
                if film:
                    film.grab(w, f"{name}  step {w.t}  {phase}  box {w.L:.0f}  sealed pockets {len(vol)}")
                if w.t % (10 * args.every) < chunk:
                    print(f"step {w.t} {row} ({time.time() - t0:.0f}s)", flush=True)
    fh.close()
    snap("final")
    if film:
        film.close()


if __name__ == "__main__":
    main()
