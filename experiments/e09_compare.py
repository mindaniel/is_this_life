"""Compare e09 runs: vesicles, clumping and sugar per setup, averaged over seeds.

python experiments/e09_compare.py runs/e09_pile runs/e09_nopile
"""
import csv
import glob
import os
import sys

import numpy as np


def load(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def stats(rows):
    ves = np.array([float(r["vesicles"]) for r in rows])
    real = np.array([float(r["real"]) for r in rows])
    largest = np.array([float(r["largest"]) for r in rows])
    c5 = np.array([float(r["c5plus"]) for r in rows])
    ponds = np.array([[float(x) for x in r["pond_vesicles"].split()] for r in rows])
    return dict(logs=len(rows), seen=(ves > 0).mean(), ves=ves.mean(), real=real.mean(), peak=ves.max(),
                largest=largest.mean(), c5=c5.mean(), per_pond=ponds.sum(axis=0))


for folder in sys.argv[1:]:
    print(f"== {folder}")
    all_ = []
    for f in sorted(glob.glob(os.path.join(folder, "*", "log.csv"))):
        s = stats(load(f))
        all_.append(s)
        print(f"  {os.path.basename(os.path.dirname(f)):14s} logs {s['logs']:3d}  with vesicle {s['seen']:.0%}  "
              f"mean vesicles {s['ves']:.2f} (real {s['real']:.2f}) peak {s['peak']:.0f}  "
              f"largest cluster {s['largest']:.0f}  C5+ {s['c5']:.1f}  per pond {s['per_pond'].astype(int)}")
    if all_:
        m = {k: np.mean([s[k] for s in all_]) for k in ("seen", "ves", "real", "largest", "c5")}
        print(f"  MEAN           with vesicle {m['seen']:.0%}  mean vesicles {m['ves']:.2f} (real {m['real']:.2f})  "
              f"largest cluster {m['largest']:.0f}  C5+ {m['c5']:.1f}")
