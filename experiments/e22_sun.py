"""e22: a "habitable zone" for our pond. The same world (two bags + a concentrated HCN puddle that
makes oily tar, as e21) on planets at different distances from the sun: closer = brighter and
hotter, further = dimmer and colder (alife/chem/sun.py). Smooth day, temperature follows.

Per distance d and seed: (1) HCN puddle made at that planet's mean temperature for 20k steps
(rates calibrated at Earth's, so heat speeds them); (2) joined to the e20 bag world (made at
Earth's temperature, then moved: the planet changes at the join, flagged); (3) 60k steps with the
sun, tar 0.8, light heat + photo-acid.
Distances 0.7 (~Venus), 0.85, 1 (Earth), 1.2, 1.5 (~Mars) -> mean kT 1.31 .. 0.90.
Logs in runs/stdout/e22_*.txt. --warm: photon heat as a warm spot (LightPump warm), not a kick;
the kicked runs crashed at strong light (0.7, one 0.85). Logs then in runs/stdout/e22w_*.txt.

    python experiments/e22_sun.py --seeds 0 1 2
"""
import argparse
import os

from e20_glow import BASE, BAGS, run_all

KT_EARTH = 1.1


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    p.add_argument("--d", type=float, nargs="+", default=[0.7, 0.85, 1.0, 1.2, 1.5])
    p.add_argument("--pre", type=int, default=20000)
    p.add_argument("--steps", type=int, default=60000)
    p.add_argument("--out", default="runs/e22")
    p.add_argument("--warm", action="store_true", help="photon heat as a warm spot (e16 --glow-warm)")
    a = p.parse_args()
    pre = os.path.join(a.out, "pre")
    os.makedirs(pre, exist_ok=True)
    kT = {d: round(KT_EARTH * d ** -0.5, 4) for d in a.d}
    ponds, joins, runs = [], [], []
    for d in a.d:
        for s in a.seeds:
            pond = os.path.join(pre, f"hcnpond_d{d:g}_s{s}.pkl")
            joined = os.path.join(pre, f"sun{d:g}_{s}_state.pkl")
            if not os.path.exists(pond):
                ponds.append((f"pond_d{d:g}_s{s}", ["experiments/hcn_pond.py", "--x", "24", "--hcn", "4000", "--nh3", "1300",
                                                    "--kT", str(kT[d]), "--kT0", str(KT_EARTH), "--steps", str(a.pre),
                                                    "--seed", str(s), "--out", pond]))
            if not os.path.exists(joined):
                joins.append((f"join_d{d:g}_s{s}", ["experiments/merge_worlds.py",
                                                    os.path.join("runs", "e20", "pre", BAGS.format(s=s)), pond, joined,
                                                    "--seed", str(s), "--kT", str(kT[d])]))
            runs.append((f"sun{d:g}_s{s}", ["experiments/e16_catalyst.py", "--start", joined, *BASE, "--steps", str(a.steps),
                                            "--seed", str(s), "--out", a.out, "--tar-len", "5", "--tar-oily", "0.8",
                                            "--glow", "--glow-heat", "10", "--glow-gain", "5",
                                            "--sun-d", str(d), "--kT-earth", str(KT_EARTH)]
                         + (["--glow-warm"] if a.warm else [])))
    run_all(ponds, a.out, "e22")
    run_all(joins, a.out, "e22")
    run_all(runs, a.out, "e22w" if a.warm else "e22")
    print("e22 done")


if __name__ == "__main__":
    main()
