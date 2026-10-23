"""Drawing with matplotlib: live window or saved snapshots."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgb, to_rgba

from .chem.elements import COLOR, INDEX, NAMES
from .chem.world import TAIL

WATER = "#E6F1FB"
BOND_COLOR = "#888780"
# sugar carbons coloured by chain length, so formaldehyde -> sugar is visible
CARB_COLORS = ((1, "#C9B3E6", "formaldehyde C1"), (2, "#8E5BB5", "small sugar C2-4"), (5, "#4B1C7A", "big sugar C5+"))


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
        self.size = size
        self.beads = self.ax.scatter(*world._wrapped().T, zorder=3)
        self._style()
        # circles around detected vesicles, so they can be checked by eye
        self.rings = self.ax.scatter([], [], s=600 * (64 / world.L) ** 2 * (size / 7) ** 2 * 4,
                                     facecolors="none", edgecolors="#E24B4A", linewidths=2, zorder=4)
        self.text = self.ax.text(0.01, 1.01, "", transform=self.ax.transAxes, fontsize=10)
        self.update()

    def _style(self):
        """Colour and size of each bead (redone when reactions add or join beads)."""
        w = self.world
        colors = np.array(COLOR, dtype=object)[w.type]
        carb = np.flatnonzero(w.type == INDEX["carb"])
        if len(carb):
            chain = np.bincount(w.mol[carb], minlength=w.n_mol)[w.mol[carb]]
            for least, color, _ in CARB_COLORS:
                colors[carb[chain >= least]] = color
        sizes = np.where(w.type == TAIL, 22, 40)
        sizes = np.where(w.rules.charge[w.type] != 0, 30 + 15 * np.abs(w.rules.charge[w.type]), sizes)
        sizes = sizes * (64 / w.L) ** 2 * (self.size / 7) ** 2   # same look in bigger boxes
        rgba = np.array([to_rgba(c) for c in colors]).reshape(-1, 4)
        crust = w.frozen & ~np.isin(w.type, [INDEX["rock"], INDEX["mineral"]])
        rgba[crust, 3] = 0.3   # dried onto the shore: faded
        self.beads.set_facecolors(rgba)
        self.beads.set_sizes(sizes)

    def legend(self):
        """Colour key for the bead types in this world, below the picture."""
        from matplotlib.lines import Line2D
        present = set(NAMES[t] for t in np.unique(self.world.type))
        items = [(COLOR[INDEX[n]], n) for n in ("head", "tail", "fa_head", "Na", "K", "Cl", "Mg", "Ca", "sugar")
                 if n in present]
        items += [(c, label) for _, c, label in CARB_COLORS]
        handles = [Line2D([], [], marker="o", ls="", color=c, markersize=7, label=label) for c, label in items]
        handles.append(Line2D([], [], marker="o", ls="", mfc="none", mec="#E24B4A", markersize=11, label="vesicle"))
        self.ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.01), ncol=5,
                       fontsize=8, frameon=False)

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
        if self.world.n != len(self.beads.get_sizes()) or (self.world.type == INDEX["carb"]).any() or self.world.frozen.any():
            self._style()
        xy = (info or {}).get("vesicle_xy") or []
        self.rings.set_offsets(np.array(xy).reshape(-1, 2))
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
        return self.bonds, self.beads, self.rings, self.text

    def save(self, path, info=None):
        self.update(info)
        self.fig.savefig(path, dpi=100, bbox_inches="tight")


class Movie:
    """Save a View as an mp4 (needs ffmpeg) or .gif, one grab() per frame."""

    def __init__(self, view, path, fps=30, dpi=100):
        from matplotlib import animation
        if path.endswith(".gif"):
            self.writer = animation.PillowWriter(fps=fps)
        else:
            import shutil
            if not shutil.which("ffmpeg"):
                raise SystemExit("ffmpeg not found: install it, or save a .gif instead")
            self.writer = animation.FFMpegWriter(fps=fps, codec="libx264", extra_args=["-pix_fmt", "yuv420p"])
        self.fig = view.fig
        self.fig.tight_layout()   # keep the title and colour key inside the frame
        self.writer.setup(self.fig, path, dpi=dpi)
        self.frames = 0

    def grab(self):
        self.writer.grab_frame()
        self.frames += 1

    def close(self):
        self.writer.finish()
