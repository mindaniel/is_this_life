"""e20: pigment made by the pond itself (HCN polymer) next to two bags, in light.

Per seed: (1) the e16 bag world runs 20k steps and is saved; a concentrated HCN puddle
(hcn_pond.py, 4000 HCN + 1300 NH3 in 24 x 48 x 48) runs 20k steps at the same time, and a plain
puddle (sugar only) too. (2) Joined (merge_worlds.py): bags + HCN puddle ("dry": HCN was
concentrated by drying), bags + plain puddle ("wet": rain kept it dilute; the same 4000 HCN +
1300 NH3 are spread through all the water outside the bags instead). (3) Carry on 60k steps:
  dark  - HCN puddle, no light
  heat  - HCN puddle, dark chains absorb light and warm their surroundings
  acid  - as heat, plus light puts H+ on a nearby COO- (photo-acid, speculative)
  wet   - dilute HCN everywhere, light + heat
Each stage runs one process per job; logs in runs/stdout/e20_*.txt.

    python experiments/e20_glow.py --seeds 0 1 2
"""
import argparse
import os
import subprocess
import sys

PY = sys.executable
BASE = ["--p-join", "0", "--fa-only", "--length", "1.2", "--fa-lengths", "1:0.25,2:0.5,3:0.25", "--ph", "5",
        "--hbond", "1", "--sealed-h"]
BAGS = "cats20_p0_nc0_br0_mix0.5_q0.5_f0x1_sh0_fa_len1.2_mix1-2-3_ph5_pka4.8_hb1_sealH_s{s}_state.pkl"


def run_all(jobs, out):
    procs = []
    for name, cmd in jobs:
        f = open(os.path.join("runs", "stdout", f"e20_{name}.txt"), "w")
        procs.append((name, subprocess.Popen([PY] + cmd, stdout=f, stderr=subprocess.STDOUT), f))
    bad = []
    for name, pr, f in procs:
        if pr.wait():
            bad.append(name)
        f.close()
    if bad:
        sys.exit(f"failed: {bad}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    p.add_argument("--pre", type=int, default=20000)
    p.add_argument("--steps", type=int, default=60000)
    p.add_argument("--out", default="runs/e20")
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    os.makedirs(os.path.join("runs", "stdout"), exist_ok=True)
    o = a.out
    stage1 = []
    for s in a.seeds:
        if not os.path.exists(os.path.join(o, "pre", BAGS.format(s=s))):
            stage1.append((f"bags_s{s}", ["experiments/e16_catalyst.py", *BASE, "--steps", str(a.pre), "--seed", str(s),
                                          "--out", os.path.join(o, "pre")]))
        for kind, n_h, n_n in (("hcnpond", 4000, 1300), ("plainpond", 0, 0)):
            if not os.path.exists(os.path.join(o, "pre", f"{kind}_s{s}.pkl")):
                stage1.append((f"{kind}_s{s}", ["experiments/hcn_pond.py", "--x", "24", "--hcn", str(n_h),
                                                "--nh3", str(n_n), "--steps", str(a.pre), "--seed", str(s),
                                                "--out", os.path.join(o, "pre", f"{kind}_s{s}.pkl")]))
    run_all(stage1, o)
    merges = []
    for s in a.seeds:
        bags = os.path.join(o, "pre", BAGS.format(s=s))
        for kind, tag in (("hcnpond", "dry"), ("plainpond", "wet")):
            merges.append((f"join_{tag}_s{s}", ["experiments/merge_worlds.py", bags,
                                                os.path.join(o, "pre", f"{kind}_s{s}.pkl"),
                                                os.path.join(o, "pre", f"{tag}{s}_state.pkl"), "--seed", str(s)]))
    run_all(merges, o)
    branches = []
    for s in a.seeds:
        dry, wet = (os.path.join(o, "pre", f"{t}{s}_state.pkl") for t in ("dry", "wet"))
        common = BASE + ["--steps", str(a.steps), "--seed", str(s), "--out", o]
        branches += [(f"dark_s{s}", ["experiments/e16_catalyst.py", "--start", dry, *common]),
                     (f"heat_s{s}", ["experiments/e16_catalyst.py", "--start", dry, "--glow", "--glow-heat", "10", *common]),
                     (f"acid_s{s}", ["experiments/e16_catalyst.py", "--start", dry, "--glow", "--glow-heat", "10",
                                     "--glow-gain", "5", *common]),
                     (f"wet_s{s}", ["experiments/e16_catalyst.py", "--start", wet, "--hcn", "4000", "--nh3", "1300",
                                    "--glow", "--glow-heat", "10", *common])]
    run_all(branches, o)
    print("e20 done")


if __name__ == "__main__":
    main()
