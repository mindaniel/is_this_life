"""Measure what self-assembled."""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from .world import TAIL


def cluster_labels(world, contact=1.5):
    """Cluster id per molecule: molecules whose tails touch share an id."""
    tails = np.flatnonzero(world.type == TAIL)
    tree = cKDTree(world._wrapped()[tails], boxsize=world.L)
    pairs = tree.query_pairs(contact, output_type="ndarray")
    m = world.mol[tails]
    rows = np.concatenate([m[pairs[:, 0]], m])   # second part links each molecule to itself
    cols = np.concatenate([m[pairs[:, 1]], m])
    graph = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(world.n_mol, world.n_mol))
    return connected_components(graph, directed=False)[1]


def tail_clusters(world, contact=1.5):
    """Cluster sizes (in molecules with tails; ions don't count), largest first."""
    if not np.any(world.type == TAIL):
        return np.array([], dtype=int)
    lipids = np.unique(world.mol[world.type == TAIL])
    sizes = np.bincount(cluster_labels(world, contact)[lipids])
    return np.sort(sizes[sizes > 0])[::-1]


def water_pockets(world, cell=1.0, reach=1.0, min_area=20, one_cluster=0.8, centers=False, masks=False):
    """Find water fully enclosed by a membrane: the inside of a vesicle (a ring in 2D).

    Paint a grid: cells within `reach` of any bead are 'wall', the rest is water.
    The biggest water region is the open sea. Another water region counts as a
    vesicle if it has at least `min_area` cells and one tail cluster makes at
    least `one_cluster` of its wall (a gap between jammed clumps is not a vesicle).
    Returns pocket areas, largest first. centers=True also returns the
    (x, y) middle of each pocket, for marking them on pictures.
    masks=True also returns each pocket's grid mask (for `lining`).
    """
    from scipy import ndimage

    m = int(np.ceil(world.L / cell))
    wall = np.zeros((m, m), dtype=bool)
    owner = np.full((m, m), -1)
    if not np.any(world.type == TAIL):
        empty = np.array([], dtype=int)
        return (empty, [], []) if masks else (empty, []) if centers else empty
    # only lipids build walls: a loose ion floating in water seals nothing
    has_tail = np.zeros(world.n_mol, dtype=bool)
    has_tail[world.mol[world.type == TAIL]] = True
    lipid = has_tail[world.mol]
    idx = (world._wrapped()[lipid] / cell).astype(int) % m
    wall[idx[:, 0], idx[:, 1]] = True
    owner[idx[:, 0], idx[:, 1]] = cluster_labels(world)[world.mol[lipid]]
    r = int(np.ceil(reach / cell))
    disk = np.hypot(*np.mgrid[-r:r + 1, -r:r + 1]) <= reach / cell
    wall = ndimage.binary_dilation(np.pad(wall, r, mode="wrap"), disk)[r:-r, r:-r]

    labels, n = ndimage.label(~wall)
    # periodic box: water touching across opposite edges is the same region
    parent = np.arange(n + 1)

    def root(a):
        while parent[a] != a:
            a = parent[a]
        return a

    for a_edge, b_edge in ((labels[0], labels[-1]), (labels[:, 0], labels[:, -1])):
        for a, b in zip(a_edge, b_edge):
            if a and b:
                parent[root(a)] = root(b)
    roots = np.array([root(a) for a in range(n + 1)])
    areas = np.bincount(roots[labels.ravel()], minlength=n + 1)
    areas[0] = 0   # label 0 is wall
    sea = np.argmax(areas)
    pockets, where, shapes = [], [], []
    region = roots[labels]
    ring = int(np.ceil(reach / cell)) + 1
    for k, a in enumerate(areas):
        if k == sea or a < min_area:
            continue
        mask = np.pad(region == k, ring, mode="wrap")
        border = ndimage.binary_dilation(mask, iterations=ring)[ring:-ring, ring:-ring]
        walls = owner[border & (owner >= 0)]
        if len(walls) and np.bincount(walls).max() >= one_cluster * len(walls):
            pockets.append(a)
            shapes.append(region == k)
            # grid axis 0 is x, axis 1 is y; circular mean handles wrap-around
            ang = 2 * np.pi * np.argwhere(region == k) / m
            mid = np.angle(np.mean(np.exp(1j * ang), axis=0)) % (2 * np.pi)
            where.append(tuple((mid / (2 * np.pi) * m + 0.5) * cell))
    order = np.argsort(pockets)[::-1]
    pockets = np.array(pockets, dtype=int)[order]
    if masks:
        return pockets, [where[k] for k in order], [shapes[k] for k in order]
    if centers:
        return pockets, [where[k] for k in order]
    return pockets


def lining(world, mask, cell=1.0, near=2.5):
    """How membrane-like the wall around one pocket is: -1..1.

    For each lipid whose head is within `near` of the pocket, compare its
    direction (last tail bead -> head) with the direction to the nearest
    pocket cell. A real vesicle's inner leaflet points heads into the water
    inside: score near +1. A gap in a tangle has lipids pointing every way:
    near 0. Tails poking into the hole make it negative.
    Also returns coverage: lining lipids per unit of pocket edge. A packed
    inner leaflet is ~1 (one lipid per bead width); a hole whose sides are
    tail ends and gaps is much lower.
    Returns (score, coverage).
    """
    from scipy import ndimage

    m = mask.shape[0]
    pad = int(np.ceil(near / cell)) + 2
    big = np.pad(mask, pad, mode="wrap")
    dist, (ix, iy) = ndimage.distance_transform_edt(~big, return_indices=True)
    p = world._wrapped()
    has_tail = np.zeros(world.n_mol, dtype=bool)
    has_tail[world.mol[world.type == TAIL]] = True
    first = np.r_[0, np.flatnonzero(np.diff(world.mol)) + 1]   # first bead of each molecule
    last = np.r_[first[1:] - 1, world.n - 1]
    lip = np.flatnonzero(has_tail[world.mol[first]])
    heads, ends = first[lip], last[lip]
    g = (p[heads] / cell).astype(int) % m + pad
    close = dist[g[:, 0], g[:, 1]] * cell <= near
    edge = np.sum(mask & ~(np.roll(mask, 1, 0) & np.roll(mask, -1, 0)
                           & np.roll(mask, 1, 1) & np.roll(mask, -1, 1))) * cell
    if not np.any(close):
        return 0.0, 0.0
    g, heads, ends = g[close], heads[close], ends[close]
    target = (np.stack([ix[g[:, 0], g[:, 1]], iy[g[:, 0], g[:, 1]]], axis=1) - pad + 0.5) * cell
    to_pocket = target - p[heads]
    to_pocket -= world.L * np.round(to_pocket / world.L)
    axis = world.delta(ends, heads)
    nt = np.linalg.norm(to_pocket, axis=1)
    inside = nt < 0.5 * cell   # head already sits on a pocket cell: count it as facing in
    cos = np.where(inside, 1.0, np.einsum("ij,ij->i", axis, to_pocket)
                   / np.maximum(np.linalg.norm(axis, axis=1) * nt, 1e-9))
    return float(np.mean(cos)), len(cos) / max(edge, 1.0)


def summary(world):
    sizes = tail_clusters(world)
    alone = int(np.sum(sizes == 1))
    pockets, where = water_pockets(world, centers=True)
    return {
        "step": world.t,
        "clusters": int(np.sum(sizes > 1)),
        "largest": int(sizes[0]) if len(sizes) else 0,
        "free_lipids": alone,
        "vesicles": len(pockets),
        "vesicle_xy": where,
        "vesicle_area": [int(a) for a in pockets],
        "energy_per_lipid": float(world.energy / max(int(np.sum(sizes)), 1)),
    }
