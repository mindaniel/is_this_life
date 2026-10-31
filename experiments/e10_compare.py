"""Summarise e10 runs: how long the bag stayed closed, wall growth, splits.

python experiments/e10_compare.py runs/e10_scan
"""
import csv
import glob
import os
import re
import sys
from collections import defaultdict

import numpy as np

groups = defaultdict(list)
for f in sorted(glob.glob(os.path.join(sys.argv[1], "*", "log.csv"))):
    name = os.path.basename(os.path.dirname(f))
    rows = list(csv.DictReader(open(f)))
    if not rows:
        continue
    closed = [int(r["step"]) for r in rows if int(r["pockets"]) > 0]
    groups[re.sub(r"_seed\d+$", "", name)].append(dict(
        last=max(closed or [0]), split=max(int(r["pockets"]) for r in rows),
        wall0=int(rows[0]["biggest"]), wall_max=max(int(r["biggest"]) for r in rows),
        wall_end=int(rows[-1]["biggest"]), two_end=int(rows[-1]["biggest_2"]),
        loose=int(rows[-1]["loose"]), steps=int(rows[-1]["step"]),
        area_max=max(max([int(a) for a in r["areas"].split()] or [0]) for r in rows)))

for g, runs in groups.items():
    f = lambda k: " ".join(f"{r[k]:>6}" for r in runs)
    print(f"== {g}  ({len(runs)} seeds, {runs[0]['steps']} steps)")
    print(f"  bag closed until  {f('last')}")
    print(f"  most bags at once {f('split')}")
    print(f"  wall start / max / end  " + "  ".join(f"{r['wall0']}/{r['wall_max']}/{r['wall_end']}" for r in runs))
    print(f"  2-tail in wall at end {f('two_end')}   loose at end {f('loose')}   biggest bag area {f('area_max')}")
