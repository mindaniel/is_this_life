"""Molecule builders."""
import numpy as np
from scipy.spatial import cKDTree

from .world import HEAD, TAIL, TYPE_NAMES


def add_lipid(world, x, y, angle, n_tails=2, head="head"):
    """Head bead, then a straight chain of tail beads pointing along angle.

    head="head" is a neutral lipid; head="fa_head" is a fatty acid (charge -1),
    the simple soap-like lipid early Earth could make.
    """
    r0 = world.rules.bond_r0
    direction = np.array([np.cos(angle), np.sin(angle)])
    n = 1 + n_tails
    positions = [np.array([x, y]) + k * r0 * direction for k in range(n)]
    types = [TYPE_NAMES.index(head) if isinstance(head, str) else head] + [TAIL] * n_tails
    bonds = [(k, k + 1) for k in range(n - 1)]
    bends = [(k, k + 2) for k in range(n - 2)]
    world.add_molecule(positions, types, bonds, bends)


def fill_with_lipids(world, n_lipids, n_tails=2, spacing=None, head="head"):
    """Scatter lipids on a grid with random orientation. Grid avoids overlaps."""
    length = n_tails * world.rules.bond_r0
    spacing = spacing or length + 1.2
    per_row = int(world.L / spacing + 1e-9)
    free = np.arange(per_row ** 2)
    if world.n:
        # skip grid cells that already hold something (rocks, ions)
        centers = (np.stack([free % per_row, free // per_row], axis=1) + 0.5) * spacing
        dist, _ = cKDTree(world._wrapped(), boxsize=world.L).query(centers)
        free = free[dist > 0.5 * spacing + 0.5]
    if len(free) < n_lipids:
        raise ValueError(f"box {world.L} fits only {len(free)} lipids at spacing {spacing}")
    cells = world.rng.permutation(free)[:n_lipids]
    for c in cells:
        cx = (c % per_row + 0.5) * spacing
        cy = (c // per_row + 0.5) * spacing
        angle = world.rng.uniform(0, 2 * np.pi)
        # centre the lipid in its cell
        x = cx - 0.5 * length * np.cos(angle)
        y = cy - 0.5 * length * np.sin(angle)
        add_lipid(world, x, y, angle, n_tails, head)


def add_ions(world, name, count, min_gap=1.0, tries=200):
    """Drop single-bead ions ("Na", "Ca", ...) at random free spots."""
    placed = 0
    for _ in range(count):
        for _ in range(tries):
            p = world.rng.uniform(0, world.L, 2)
            if not world.n:
                break
            dist, _ = cKDTree(world._wrapped(), boxsize=world.L).query(p)
            if dist >= min_gap:
                break
        else:
            raise ValueError(f"no free room for {name} ion {placed + 1}")
        world.add_molecule([p], [name])
        placed += 1


def add_mineral(world, y, spacing=1.0):
    """Flat frozen mineral surface across the box at height y.

    Clay surfaces carry negative charge sites. Each site is one frozen bead.
    """
    xs = np.arange(0, world.L, spacing)
    world.add_molecule([[x, y] for x in xs], ["mineral"] * len(xs), frozen=True)


def add_rock_wall(world, start, end, gaps=(), spacing=0.9):
    """Frozen rock beads in a line from start to end.

    Beads sit closer than their size, so nothing slips between them.
    gaps: list of (from, to) fractions along the wall left open, e.g.
    [(0.45, 0.55)] leaves a hole in the middle. Small gap = semi-closed pool.
    """
    start, end = np.asarray(start, float), np.asarray(end, float)
    length = np.linalg.norm(end - start)
    fracs = np.arange(0, length, spacing) / length
    keep = [f for f in fracs if not any(a <= f <= b for a, b in gaps)]
    if keep:
        world.add_molecule([start + f * (end - start) for f in keep], ["rock"] * len(keep), frozen=True)


def net_charge(world):
    return float(np.sum(world.rules.charge[world.type]))
