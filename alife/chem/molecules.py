"""Molecule builders."""
import numpy as np
from scipy.spatial import cKDTree

from .world import HEAD, TAIL, TYPE_NAMES


def add_lipid(world, x, y, angle, n_tails=2, head="head", tail=TAIL):
    """Head bead, then a straight chain of tail beads pointing along angle.

    head="head" is a neutral lipid; head="fa_head" is a fatty acid (charge -1),
    the simple soap-like lipid early Earth could make.
    """
    r0 = world.rules.bond_r0
    direction = np.array([np.cos(angle), np.sin(angle)])
    n = 1 + n_tails
    positions = [np.array([x, y]) + k * r0 * direction for k in range(n)]
    tail = TYPE_NAMES.index(tail) if isinstance(tail, str) else tail
    types = [TYPE_NAMES.index(head) if isinstance(head, str) else head] + [tail] * n_tails
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


# ---- 3D (World(dim=3)) ------------------------------------------------------
def add_lipid3(world, base, direction, n_tails=2, head="head", tail=TAIL):
    """3D lipid: head at `base`, tail beads along the unit vector `direction`."""
    r0 = world.rules.bond_r0
    d = np.asarray(direction, float) / np.linalg.norm(direction)
    n = 1 + n_tails
    tail = TYPE_NAMES.index(tail) if isinstance(tail, str) else tail
    types = [TYPE_NAMES.index(head) if isinstance(head, str) else head] + [tail] * n_tails
    world.add_molecule([np.asarray(base, float) + k * r0 * d for k in range(n)], types,
                       [(k, k + 1) for k in range(n - 1)], [(k, k + 2) for k in range(n - 2)])


def bilayer_disc3(world, center, radius, n_tails=2, spacing=1.1, head="head", flat_box=False):
    """Flat round bilayer patch lying in the x-y plane (two leaflets, tails inside).

    Lipids sit on a hexagonal lattice, `spacing` apart. flat_box=True fills the
    whole x-y face of the box instead (a bilayer with no edge, through the
    periodic walls). Returns the number of lipids.
    """
    center = np.asarray(center, float)
    r0 = world.rules.bond_r0
    rows = int(world.L / (spacing * np.sqrt(3) / 2)) if flat_box else int(2 * radius / spacing) + 2
    cols = int(world.L / spacing) if flat_box else rows
    count = 0
    for a in range(rows):
        for b in range(cols):
            x = (b + 0.5 * (a % 2)) * spacing
            y = a * spacing * np.sqrt(3) / 2
            if flat_box:
                xy = np.array([x, y])
            else:
                xy = center[:2] + np.array([x, y]) - (rows - 1) * spacing / 2
                if np.hypot(*(xy - center[:2])) > radius:
                    continue
            for side in (1, -1):   # upper leaflet points down, lower points up
                top = center[2] + side * (0.5 + n_tails * r0)
                add_lipid3(world, [xy[0], xy[1], top], [0, 0, -side], n_tails, head)
                count += 1
    return count


def gas3(world, n_lipids, n_tails=2, head="head"):
    """Lipids scattered through the box on a cubic lattice, random directions."""
    length = n_tails * world.rules.bond_r0
    per = int(world.L / (length + 1.2))
    if per ** 3 < n_lipids:
        raise ValueError(f"box {world.L} fits only {per ** 3} lipids")
    spacing = world.L / per
    for c in world.rng.permutation(per ** 3)[:n_lipids]:
        mid = (np.array([c % per, c // per % per, c // per ** 2]) + 0.5) * spacing
        d = world.rng.normal(size=3)
        d /= np.linalg.norm(d)
        add_lipid3(world, mid - 0.5 * length * d, d, n_tails, head)


def vesicle3(world, center, radius, n_tails=2, area=1.2, head="head"):
    """Ready-made spherical vesicle: two leaflets around the bilayer middle at `radius`.

    Lipids spread evenly (Fibonacci sphere), `area` sigma^2 each. Returns the count.
    """
    center = np.asarray(center, float)
    r0 = world.rules.bond_r0
    count = 0
    for side in (1, -1):   # outer leaflet: heads out; inner: heads in
        r_head = radius + side * (0.5 + n_tails * r0)
        n = int(4 * np.pi * r_head ** 2 / area)
        k = np.arange(n) + 0.5
        z = 1 - 2 * k / n
        phi = np.pi * (1 + 5 ** 0.5) * k
        u = np.stack([np.sqrt(1 - z ** 2) * np.cos(phi), np.sqrt(1 - z ** 2) * np.sin(phi), z], 1)
        for d in u:
            add_lipid3(world, center + r_head * d, -side * d, n_tails, head)
        count += n
    return count
