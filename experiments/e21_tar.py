"""e21: the pond's HCN polymer as oily tar (long HCN polymer doesn't dissolve). Does it sink into
the bag walls by itself, and does light on it then change a bag?

Starts from e20's joined worlds (runs/e20/pre/dry{s}_state.pkl: bags + concentrated HCN puddle).
HCN chains of 5+ turn into tar beads (nitrogen.Nitrogen tar_len); tar oiliness is a guess, so two
values. 60k steps:
  tar0.5     - tar oiliness 0.5, light: heat + photo-acid
  tar0.8     - tar oiliness 0.8, light: heat + photo-acid
  tar0.8dark - tar oiliness 0.8, no light (control)
Logs in runs/stdout/e21_*.txt.

    python experiments/e21_tar.py --seeds 0 1 2
"""
import argparse
import os
import sys

from e20_glow import BASE, run_all


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    p.add_argument("--steps", type=int, default=60000)
    p.add_argument("--out", default="runs/e21")
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    jobs = []
    for s in a.seeds:
        start = ["experiments/e16_catalyst.py", "--start", os.path.join("runs", "e20", "pre", f"dry{s}_state.pkl"),
                 *BASE, "--steps", str(a.steps), "--seed", str(s), "--out", a.out, "--tar-len", "5"]
        light = ["--glow", "--glow-heat", "10", "--glow-gain", "5"]
        jobs += [(f"tar0.5_s{s}", start + ["--tar-oily", "0.5"] + light),
                 (f"tar0.8_s{s}", start + ["--tar-oily", "0.8"] + light),
                 (f"tar0.8dark_s{s}", start + ["--tar-oily", "0.8"])]
    run_all(jobs, a.out, "e21")
    print("e21 done")


if __name__ == "__main__":
    main()
