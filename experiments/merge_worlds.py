"""Join two saved worlds side by side (two waters meeting), then let the seam settle.

The box sides across the seam must match. Each world is cut open at its emptiest slab, so
bags and clumps stay whole; b goes after a along --axis. a's rules, time and bag markers come
first (e16 reads a's two bags as A and B). Then 200 capped steps (beads at the seam may touch).

    python experiments/merge_worlds.py runs/e20/bags_state.pkl runs/e20/pond.pkl runs/e20/joined.pkl
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from alife.chem.world import World   # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument("a")
    p.add_argument("b")
    p.add_argument("out")
    p.add_argument("--axis", type=int, default=0)
    p.add_argument("--seed", type=int, default=0, help="random stream of the joined world")
    p.add_argument("--kT", type=float, default=None, help="set both worlds to this temperature first")
    args = p.parse_args()
    a, _ = World.load(args.a)
    b, _ = World.load(args.b)
    if args.kT:
        a.kT = b.kT = args.kT
    w = World.merge(a, b, axis=args.axis, seed=args.seed)
    t = w.t
    w.step(200, max_force=50)
    w.t = t
    w.save(args.out, joined=(args.a, args.b))
    print(f"joined: box {w.box()}, {w.n} beads ({a.n} + {b.n}), t {w.t}")


if __name__ == "__main__":
    main()
