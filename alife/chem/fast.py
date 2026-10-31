"""Compiled force loops (numba). Same math as World.forces in numpy.

Work in 2D or 3D (pos with 2 or 3 columns); 2D gives the same numbers as before.

One pass over each pair list instead of many temporary numpy arrays:
~5x faster for big worlds. If numba is missing, World falls back to numpy.
"""
import math

import numpy as np

try:
    from numba import njit
    HAVE_NUMBA = True
except ImportError:   # pragma: no cover
    HAVE_NUMBA = False

    def njit(*a, **k):
        return (lambda f: f) if not a or not callable(a[0]) else a[0]


@njit(cache=True)
def pair_forces(pos, L, pairs, typ, sigma, rc, eps, wc, F):
    """WCA repulsion + cos^2 oily attraction. Adds to F, returns energy."""
    three = pos.shape[1] == 3
    energy = 0.0
    for p in range(pairs.shape[0]):
        i, j = pairs[p, 0], pairs[p, 1]
        dx = pos[j, 0] - pos[i, 0]
        dy = pos[j, 1] - pos[i, 1]
        dx -= L * np.rint(dx / L)
        dy -= L * np.rint(dy / L)
        dz = 0.0
        if three:
            dz = pos[j, 2] - pos[i, 2]
            dz -= L * np.rint(dz / L)
        r2 = dx * dx + dy * dy + dz * dz
        ti, tj = typ[i], typ[j]
        c = rc[ti, tj]
        e = eps[ti, tj]
        reach = c + wc
        if r2 >= reach * reach:
            continue
        if r2 < c * c:
            b = sigma[ti, tj]
            sr2 = b * b / r2
            sr6 = sr2 * sr2 * sr2
            s = 24.0 * (2.0 * sr6 * sr6 - sr6) / r2   # = fmag / dist
            energy += 4.0 * (sr6 * sr6 - sr6) + 1.0 - e
        elif e > 0.0:
            dist = math.sqrt(r2)
            y = math.pi * (dist - c) / wc   # = 2x;  cos(x)^2 = (1 + cos 2x) / 2
            s = -e * math.pi / (2.0 * wc) * math.sin(y) / dist
            energy -= e * 0.5 * (1.0 + math.cos(y))
        else:
            continue
        F[j, 0] += s * dx
        F[j, 1] += s * dy
        F[i, 0] -= s * dx
        F[i, 1] -= s * dy
        if three:
            F[j, 2] += s * dz
            F[i, 2] -= s * dz
    return energy


@njit(cache=True)
def charge_forces(pos, L, pairs, typ, charge, coulomb, debye, cut, F):
    """Screened Coulomb, shifted to zero energy at the cutoff."""
    three = pos.shape[1] == 3
    energy = 0.0
    shift = math.exp(-cut / debye) / cut
    for p in range(pairs.shape[0]):
        i, j = pairs[p, 0], pairs[p, 1]
        dx = pos[j, 0] - pos[i, 0]
        dy = pos[j, 1] - pos[i, 1]
        dx -= L * np.rint(dx / L)
        dy -= L * np.rint(dy / L)
        dz = 0.0
        if three:
            dz = pos[j, 2] - pos[i, 2]
            dz -= L * np.rint(dz / L)
        dist = max(math.sqrt(dx * dx + dy * dy + dz * dz), 0.5)
        if dist >= cut:
            continue
        qq = coulomb * charge[typ[i]] * charge[typ[j]]
        screen = math.exp(-dist / debye)
        fmag = qq * screen * (1.0 / dist ** 2 + 1.0 / (debye * dist))
        energy += qq * (screen / dist - shift)
        s = fmag / dist
        F[j, 0] += s * dx
        F[j, 1] += s * dy
        F[i, 0] -= s * dx
        F[i, 1] -= s * dy
        if three:
            F[j, 2] += s * dz
            F[i, 2] -= s * dz
    return energy


@njit(cache=True)
def spring_forces(pos, L, idx, k, r0, F):
    """Harmonic springs (bonds and bends)."""
    three = pos.shape[1] == 3
    energy = 0.0
    for p in range(idx.shape[0]):
        i, j = idx[p, 0], idx[p, 1]
        dx = pos[j, 0] - pos[i, 0]
        dy = pos[j, 1] - pos[i, 1]
        dx -= L * np.rint(dx / L)
        dy -= L * np.rint(dy / L)
        dz = 0.0
        if three:
            dz = pos[j, 2] - pos[i, 2]
            dz -= L * np.rint(dz / L)
        dist = math.sqrt(dx * dx + dy * dy + dz * dz)
        stretch = dist - r0
        s = -k * stretch / max(dist, 1e-9)
        F[j, 0] += s * dx
        F[j, 1] += s * dy
        F[i, 0] -= s * dx
        F[i, 1] -= s * dy
        if three:
            F[j, 2] += s * dz
            F[i, 2] -= s * dz
        energy += 0.5 * k * stretch * stretch
    return energy


@njit(cache=True)
def cell_pairs(p, L, cutoff, mol, frozen, sel):
    """All pairs (i < j by list order) of beads in `sel` closer than cutoff.

    Cell list: beads sorted into square cells at least `cutoff` wide, each
    cell checked against itself and 4 neighbours (half the stencil, so every
    pair is found once). Same-molecule and rock-rock pairs are skipped.
    Needs at least 3 cells per side. p must be wrapped into [0, L).
    """
    nc = int(L // cutoff)
    size = L / nc
    n = sel.shape[0]
    cell = np.empty(n, np.int64)
    count = np.zeros(nc * nc + 1, np.int64)
    for a in range(n):
        cx = min(int(p[sel[a], 0] / size), nc - 1)
        cy = min(int(p[sel[a], 1] / size), nc - 1)
        cell[a] = cx * nc + cy
        count[cell[a] + 1] += 1
    start = np.cumsum(count)
    fill = start[:-1].copy()
    order = np.empty(n, np.int64)
    for a in range(n):
        order[fill[cell[a]]] = sel[a]
        fill[cell[a]] += 1
    cap = 16 * n + 1024
    out = np.empty((cap, 2), np.int64)
    m = 0
    c2 = cutoff * cutoff
    for cx in range(nc):
        for cy in range(nc):
            here = cx * nc + cy
            for k in range(5):
                if k == 0:
                    ox, oy = 0, 0
                elif k == 1:
                    ox, oy = 1, 0
                elif k == 2:
                    ox, oy = 1, 1
                elif k == 3:
                    ox, oy = 0, 1
                else:
                    ox, oy = -1, 1
                other = ((cx + ox) % nc) * nc + (cy + oy) % nc
                for a in range(start[here], start[here + 1]):
                    i = order[a]
                    b0 = a + 1 if k == 0 else start[other]
                    for b in range(b0, start[other + 1]):
                        j = order[b]
                        if mol[i] == mol[j] or (frozen[i] and frozen[j]):
                            continue
                        dx = p[j, 0] - p[i, 0]
                        dy = p[j, 1] - p[i, 1]
                        dx -= L * np.rint(dx / L)
                        dy -= L * np.rint(dy / L)
                        if dx * dx + dy * dy < c2:
                            if m == cap:
                                cap *= 2
                                bigger = np.empty((cap, 2), np.int64)
                                bigger[:m] = out[:m]
                                out = bigger
                            out[m, 0] = i
                            out[m, 1] = j
                            m += 1
    return out[:m].copy()


# the 13 neighbour cells "after" a cell in 3D (with the cell itself: every pair found once)
HALF3 = np.array([(x, y, z) for x in (-1, 0, 1) for y in (-1, 0, 1) for z in (-1, 0, 1)
                  if (x, y, z) > (0, 0, 0)], np.int64)


@njit(cache=True)
def cell_pairs3(p, L, cutoff, mol, frozen, sel, half):
    """3D version of cell_pairs: cubic cells, each checked against itself and 13 neighbours."""
    nc = int(L // cutoff)
    size = L / nc
    n = sel.shape[0]
    cell = np.empty(n, np.int64)
    count = np.zeros(nc * nc * nc + 1, np.int64)
    for a in range(n):
        cx = min(int(p[sel[a], 0] / size), nc - 1)
        cy = min(int(p[sel[a], 1] / size), nc - 1)
        cz = min(int(p[sel[a], 2] / size), nc - 1)
        cell[a] = (cx * nc + cy) * nc + cz
        count[cell[a] + 1] += 1
    start = np.cumsum(count)
    fill = start[:-1].copy()
    order = np.empty(n, np.int64)
    for a in range(n):
        order[fill[cell[a]]] = sel[a]
        fill[cell[a]] += 1
    cap = 32 * n + 1024
    out = np.empty((cap, 2), np.int64)
    m = 0
    c2 = cutoff * cutoff
    for cx in range(nc):
        for cy in range(nc):
            for cz in range(nc):
                here = (cx * nc + cy) * nc + cz
                for k in range(half.shape[0] + 1):
                    if k == 0:
                        other = here
                    else:
                        other = (((cx + half[k - 1, 0]) % nc) * nc + (cy + half[k - 1, 1]) % nc) * nc                             + (cz + half[k - 1, 2]) % nc
                    for a in range(start[here], start[here + 1]):
                        i = order[a]
                        b0 = a + 1 if k == 0 else start[other]
                        for b in range(b0, start[other + 1]):
                            j = order[b]
                            if mol[i] == mol[j] or (frozen[i] and frozen[j]):
                                continue
                            dx = p[j, 0] - p[i, 0]
                            dy = p[j, 1] - p[i, 1]
                            dz = p[j, 2] - p[i, 2]
                            dx -= L * np.rint(dx / L)
                            dy -= L * np.rint(dy / L)
                            dz -= L * np.rint(dz / L)
                            if dx * dx + dy * dy + dz * dz < c2:
                                if m == cap:
                                    cap *= 2
                                    bigger = np.empty((cap, 2), np.int64)
                                    bigger[:m] = out[:m]
                                    out = bigger
                                out[m, 0] = i
                                out[m, 1] = j
                                m += 1
    return out[:m].copy()
