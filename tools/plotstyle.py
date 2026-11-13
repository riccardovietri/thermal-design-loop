"""plotstyle: one look for every figure in the repo.

Figures are drawn 8 in wide at 200 dpi (1600 px). Read on a page at about 700 px wide, 1 pt is about 1.2 px, so the smallest
text used here (10 pt) stays at about 12 px. White background, light grid, no top or right spine.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedLocator, MaxNLocator, NullFormatter, NullLocator

WIDTH_IN = 8.0
DPI = 200
SMALL = 10.0  # smallest font used anywhere (annotations, colour-bar ticks)

INK = "#1f2937"
MUTED = "#6b7280"
FAINT = "#d1d5db"
ACCENT = "#c2410c"  # pressure budget, recommended design
BLUE = "#1d4ed8"
THRESHOLD = "#15803d"
OBJECTIVE = "#166534"
COPPER = "#b87333"
FIN_GREY = "#9ca3af"

RC = {
    "font.size": 11,
    "axes.titlesize": 12.5,
    "axes.labelsize": 11.5,
    "xtick.labelsize": 10.5,
    "ytick.labelsize": 10.5,
    "legend.fontsize": 10.5,
    "figure.titlesize": 13,
    "axes.titlepad": 8,
    "axes.edgecolor": "#4b5563",
    "axes.labelcolor": INK,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.axisbelow": True,
    "grid.color": "#e5e7eb",
    "grid.linewidth": 0.8,
    "xtick.color": "#374151",
    "ytick.color": "#374151",
    "lines.linewidth": 1.8,
    "lines.markersize": 6,
    "legend.frameon": False,
    "legend.handlelength": 1.6,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.dpi": DPI,
    "font.family": "DejaVu Sans",
}


def apply() -> None:
    plt.rcParams.update(RC)


def figure(height: float, nrows: int = 1, ncols: int = 1, width: float = WIDTH_IN, **kw):
    apply()
    return plt.subplots(nrows, ncols, figsize=(width, height), constrained_layout=True, **kw)


def save(fig, out) -> None:
    fig.savefig(out, dpi=DPI, facecolor="white")
    plt.close(fig)


NICE_TICKS = (0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 7, 10, 12, 15, 20, 25, 30, 40, 50, 70, 100, 150, 200)


def plain_log_ticks(ax, lo: float, hi: float, axis: str = "x", ticks: tuple | None = None) -> None:
    """Readable tick labels (0.3, 1, 2, ...) on a log axis instead of matplotlib's 3×10⁻¹ style."""
    a = ax.xaxis if axis == "x" else ax.yaxis
    a.set_major_locator(FixedLocator([t for t in (ticks or NICE_TICKS) if lo <= t <= hi]))
    a.set_major_formatter(lambda v, _: f"{v:g}")
    a.set_minor_locator(NullLocator())
    a.set_minor_formatter(NullFormatter())


def field_levels(vmin: float, vmax: float, n: int = 30) -> np.ndarray:
    return np.linspace(vmin, vmax, n + 1)


def round_ticks(vmin: float, vmax: float, nbins: int = 5, min_ticks: int = 4) -> list[float]:
    """Round tick values inside [vmin, vmax], at least min_ticks of them where the range allows."""
    tol = 1e-9 * max(abs(vmin), abs(vmax), 1.0)
    for n in range(nbins, nbins + 6):
        ticks = [float(t) for t in MaxNLocator(nbins=n, steps=[1, 2, 2.5, 5, 10]).tick_values(vmin, vmax)
                 if vmin - tol <= t <= vmax + tol]
        if len(ticks) >= min_ticks:
            break
    return ticks


def colorbar(fig, mappable, ax, label: str, vmin: float, vmax: float, nbins: int = 5, min_ticks: int = 4, **kw):
    """Colour bar with ticks at round values inside [vmin, vmax] and a label that carries the unit."""
    ticks = round_ticks(vmin, vmax, nbins, min_ticks)
    cb = fig.colorbar(mappable, ax=ax, ticks=ticks, **kw)
    cb.set_label(label)
    cb.ax.tick_params(labelsize=SMALL)
    cb.ax.yaxis.set_major_formatter(lambda v, _: f"{v:g}")
    cb.outline.set_linewidth(0.6)
    return cb


def no_grid(*axes) -> None:
    for ax in axes:
        ax.grid(False)
        for s in ("top", "right"):
            ax.spines[s].set_visible(True)
