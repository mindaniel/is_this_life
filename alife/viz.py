"""Drawing with matplotlib: live window or saved snapshots."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgb

from .chem.elements import COLOR
from .chem.world import TAIL

WATER = "#E6F1FB"
BOND_COLOR = "#888780"


def _bond_segments(world):
    """Bond lines, skipping ones that wrap around the periodic edge."""
    if not len(world.bonds):
        return np.zeros((0, 2, 2))
    p = world._wrapped()
    a, b = p[world.bonds[:, 0]], p[world.bonds[:, 1]]
    keep = np.all(np.abs(a - b) < world.L / 2, axis=1)
    return np.stack([a[keep], b[keep]], axis=1)


class View:
    def __init__(self, world, title="", size=7):
        self.world = world
        self.fig, self.ax = plt.subplots(figsize=(size, size))
        self.ax.set_facecolor(WATER)
        self.ax.set_xlim(0, world.L)
        self.ax.set_ylim(0, world.L)
        self.ax.set_aspect("equal")
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.title = title
        self.bonds = LineCollection([], colors=BOND_COLOR, linewidths=1)
        self.ax.add_collection(self.bonds)
        colors = np.array(COLOR)[world.type]
        sizes = np.where(world.type == TAIL, 22, 40)
        sizes = np.where(world.rules.charge[world.type] != 0, 30 + 15 * np.abs(world.rules.charge[world.type]), sizes)
        sizes = sizes * (64 / world.L) ** 2 * (size / 7) ** 2   # same look in bigger boxes
        self.beads = self.ax.scatter(*world._wrapped().T, c=colors, s=sizes, zorder=3)
        self.text = self.ax.text(0.01, 1.01, "", transform=self.ax.transAxes, fontsize=10)
        self.update()

    def show_field(self, field, vmin, vmax, n=96, cmap="coolwarm"):
        """Paint a map under the beads, e.g. temperature: field(pos) -> value."""
        if field is None:
            return
        xs = (np.arange(n) + 0.5) * self.world.L / n
        gx, gy = np.meshgrid(xs, xs)
        values = field(np.stack([gx.ravel(), gy.ravel()], axis=1)).reshape(n, n)
        if getattr(self, "field_img", None) is None:
            self.field_img = self.ax.imshow(values, origin="lower", extent=(0, self.world.L, 0, self.world.L),
                                            cmap=cmap, vmin=vmin, vmax=vmax, alpha=0.35, zorder=0)
        else:
            self.field_img.set_data(values)

    def update(self, info=None):
        self.bonds.set_segments(_bond_segments(self.world))
        self.beads.set_offsets(self.world._wrapped())
        label = f"{self.title}  step {self.world.t}"
        if info:
            label += f"  clusters {info['clusters']}  largest {info['largest']}  free {info['free_lipids']}"
            if "vesicles" in info:
                label += f"  vesicles {info['vesicles']}"
            if "kT" in info:
                label += f"  kT {info['kT']:.2f}"
            if "light" in info:
                # night water is dark, noon water is bright
                k = 0.35 + 0.65 * info["light"]
                self.ax.set_facecolor(tuple(k * c for c in to_rgb(WATER)))
        self.text.set_text(label)
        return self.bonds, self.beads, self.text

    def save(self, path, info=None):
        self.update(info)
        self.fig.savefig(path, dpi=100, bbox_inches="tight")
