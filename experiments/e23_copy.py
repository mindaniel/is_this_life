"""e23: do RNA-like chains copy themselves? (alife/chem/copying.py)

A 3D box with a few template chains of X/Y (random order) and a soup of loose X and Y monomers.
Only pairing, joining and breaking rules (no copy rule). Questions:
  - do joins happen more often between monomers held side by side on a chain (templated joins)?
  - do new chains carry the complement of a template's order (information copied)? Compared with the
    same new chains shuffled (chance level).
  - do warm/cool cycles (strands melt apart, then re-pair) let copies be made again and again?

  C:\\Apps\\Python314\\python.exe experiments/e23_copy.py --seed 0 --steps 50000
"""
import argparse, os, sys, time
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, ".")
from alife.chem.world import World
from alife.chem.copying import Copying, BB, add_monomer, add_strand, base_of, sequence


def place(w, rng, n, gap=1.1, tries=200000):
    """n free spots at least `gap` from every bead (and each other)."""
    out = []
    p = w._wrapped()
    tree = cKDTree(p, boxsize=w.L) if len(p) else None
    for _ in range(tries):
        x = rng.uniform(0, w.L, 3)
        if tree is not None and tree.query(x)[0] < gap:
            continue
        if out and np.min(np.linalg.norm((np.array(out) - x + w.L / 2) % w.L - w.L / 2, axis=1)) < gap:
            continue
        out.append(x)
        if len(out) == n:
            return np.array(out)
    raise RuntimeError("box too full")


def add_template(w, rng, seq):
    for _ in range(1000):
        start, d = rng.uniform(0, w.L, 3), rng.normal(size=3)
        d /= np.linalg.norm(d)
        side = np.cross(d, rng.normal(size=3)); side /= np.linalg.norm(side)
        pts = (start + np.outer(np.arange(len(seq)), d)) % w.L
        if w.n and cKDTree(w._wrapped(), boxsize=w.L).query(np.vstack([pts, (pts + side) % w.L]))[0].min() < 1.5:
            continue
        add_strand(w, start, d, side, seq)
        return
    raise RuntimeError("no room for template")


def comp(s):
    return s.translate(str.maketrans("XY", "YX"))


def copy_score(s, temps):
    """Best fraction of s matching the complement of any template stretch, read either way."""
    best = 0.0
    for t in temps:
        c = comp(t)
        for r in (c, c[::-1]):
            for o in range(len(r) - len(s) + 1):
                best = max(best, sum(a == b for a, b in zip(s, r[o:o + len(s)])) / len(s))
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--box", type=float, default=20.0)
    ap.add_argument("--templates", type=int, default=6)
    ap.add_argument("--tlen", type=int, default=10)
    ap.add_argument("--monomers", type=int, default=300)
    ap.add_argument("--steps", type=int, default=50000)
    ap.add_argument("--kT", type=float, default=1.1)
    ap.add_argument("--eps", type=float, default=4.0, help="X-Y pair energy")
    ap.add_argument("--eps-mis", type=float, default=1.5, help="X-X / Y-Y pair energy")
    ap.add_argument("--k-pair", type=float, default=30.0, help="stiffness of a base pair's springs")
    ap.add_argument("--p-join", type=float, default=0.02)
    ap.add_argument("--stack", type=float, default=0.0, help="oiliness of bases (stacking: flat bases attract)")
    ap.add_argument("--p-break", type=float, default=0.0)
    ap.add_argument("--cycle", type=int, default=0, help="warm/cool period in steps (0 = steady)")
    ap.add_argument("--hot", type=float, default=1.6, help="kT in the warm half of a cycle")
    ap.add_argument("--log", type=int, default=1000)
    ap.add_argument("--out", default="runs/e23")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    w = World(box=args.box, kT=args.kT, dim=3, seed=args.seed)
    for b in ("nx", "ny"):
        w.rules.set_oily(b, args.stack)
    temps = ["".join(rng.choice(["X", "Y"], args.tlen)) for _ in range(args.templates)]
    for s in temps:
        add_template(w, rng, s)
    tbeads = np.arange(w.n)   # template beads (first ones)
    pts = place(w, rng, args.monomers, gap=2.3)
    for k, x in enumerate(pts):
        d = rng.normal(size=3)
        add_monomer(w, x, d / np.linalg.norm(d), "nx" if k % 2 == 0 else "ny")
    chem = Copying(w, eps_match=args.eps, eps_mis=args.eps_mis, p_join=args.p_join, p_break=args.p_break, k_pair=args.k_pair)
    w.step(200, max_force=50)
    name = (f"copy_t{args.templates}x{args.tlen}_m{args.monomers}_b{args.box:g}_kT{args.kT:g}_e{args.eps:g}"
            f"_mis{args.eps_mis:g}_k{args.k_pair:g}_st{args.stack:g}_j{args.p_join:g}_br{args.p_break:g}"
            + (f"_cyc{args.cycle}h{args.hot:g}" if args.cycle else "") + f"_s{args.seed}")
    print(name, "templates:", temps, flush=True)
    rows, t0 = [], time.time()
    while w.t < args.steps:
        if args.cycle:
            w.kT = args.hot if (w.t // (args.cycle // 2)) % 2 else args.kT
        w.step(chem.every)
        v = np.abs(w.vel).max(axis=1)
        if v.max() > 30:   # blowing up: say which beads
            k = np.argsort(v)[-4:]
            raise RuntimeError(f"step {w.t} fast beads {k} types {w.type[k]} v {v[k].round(1)} partner {chem.partner[k]} "
                               f"bonds {[w.bonds[(w.bonds == x).any(axis=1)].tolist() for x in k]}")
        chem.react()
        if w.t % args.log < chem.every:
            size = np.bincount(w.mol[w.type == BB], minlength=w.n_mol)
            is_t = np.zeros(w.n_mol, dtype=bool); is_t[w.mol[tbeads]] = True
            new = size[(size >= 2) & ~is_t]
            paired = (chem.partner[:w.n] >= 0)
            on_t = paired & np.isin(chem.partner[:w.n], tbeads) & ~np.isin(np.arange(w.n), tbeads)
            rows.append(dict(step=w.t, kT=w.kT, joins=chem.joins, templated=chem.templated, ok=chem.templated_ok, breaks=chem.breaks,
                             paired=round(paired[w.type != BB].mean(), 3), on_templates=int(on_t.sum()),
                             monomers=int((size == 1).sum()), new_chains=int(len(new)),
                             longest_new=int(new.max()) if len(new) else 0, mean_new=round(new.mean(), 2) if len(new) else 0))
            r = rows[-1]
            print(f"{w.t} kT {w.kT:.2f} joins {r['joins']} templated {r['templated']} ok {r['ok']} paired {r['paired']} "
                  f"on_t {r['on_templates']} mono {r['monomers']} new {r['new_chains']} longest {r['longest_new']} "
                  f"({time.time() - t0:.0f}s)", flush=True)
    pd.DataFrame(rows).to_csv(f"{args.out}/{name}.csv", index=False)
    # which new chains copy a template? compare with the same chains shuffled
    size = np.bincount(w.mol[w.type == BB], minlength=w.n_mol)
    is_t = np.zeros(w.n_mol, dtype=bool); is_t[w.mol[tbeads]] = True
    out, bo = [], base_of(w)[0]
    for m in np.flatnonzero((size >= 3) & ~is_t):
        s = sequence(w, m, bo)
        sh = np.mean([copy_score("".join(rng.permutation(list(s))), temps) for _ in range(20)])
        out.append(dict(mol=m, len=len(s), seq=s, score=copy_score(s, temps), shuffled=round(sh, 3)))
    o = pd.DataFrame(out)
    o.to_csv(f"{args.out}/{name}_chains.csv", index=False)
    if len(o):
        print(o.sort_values("len", ascending=False).head(15).to_string(index=False))
        print(f"new chains >=3: {len(o)}  match {o.score.mean():.3f}  shuffled {o.shuffled.mean():.3f}")
    w.save(f"{args.out}/{name}_state.pkl", partner=chem.partner, temps=temps)


if __name__ == "__main__":
    main()
