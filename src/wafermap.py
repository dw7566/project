from __future__ import annotations

import statistics
from collections import defaultdict
from math import ceil
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


WAFERMAP_PARAMETERS = [
    "extinction_ratio_db",
    "vpi_mean_v",
]

HIGHER_IS_BETTER = {
    "extinction_ratio_db": True,
    "vpi_mean_v": False,
}


def _scalar_float(value: object) -> float | None:
    if value is None:
        return None
    if isinstance(value, (int, float, np.integer, np.floating)):
        result = float(value)
    else:
        text = str(value).strip()
        if not text or ";" in text:
            return None
        try:
            result = float(text)
        except ValueError:
            return None
    return result if np.isfinite(result) else None


def _parameters(rows: list[dict[str, object]]) -> list[str]:
    parameters: list[str] = []
    for parameter in WAFERMAP_PARAMETERS:
        if any(_scalar_float(row.get(parameter)) is not None for row in rows):
            parameters.append(parameter)
    return parameters


def plot_wafermap(wafer: str, timestamp: str,
                  rows: list[dict[str, object]], path: Path) -> bool:
    parameters = _parameters(rows)
    if not parameters:
        return False

    values_by_parameter: dict[str, dict[tuple[int, int], list[float]]] = {
        parameter: defaultdict(list) for parameter in parameters
    }
    die_coords: set[tuple[int, int]] = set()

    for row in rows:
        try:
            die_col = int(str(row["die_column"]))
            die_row = int(str(row["die_row"]))
        except (KeyError, TypeError, ValueError):
            continue

        die_coords.add((die_col, die_row))
        for parameter in parameters:
            value = _scalar_float(row.get(parameter))
            if value is not None:
                values_by_parameter[parameter][(die_col, die_row)].append(value)

    if not die_coords:
        return False

    x_coords = sorted({die_col for die_col, _ in die_coords})
    y_coords = sorted({die_row for _, die_row in die_coords})
    x_index = {coord: index for index, coord in enumerate(x_coords)}
    y_index = {coord: index for index, coord in enumerate(y_coords)}

    ncols = min(3, len(parameters))
    nrows = ceil(len(parameters) / ncols)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(5.2 * ncols, 4.8 * nrows),
        squeeze=False,
        constrained_layout=True,
    )
    fig.suptitle(f"{wafer} / {timestamp}: wafermap", fontsize=15, fontweight="bold")

    cmap = mpl.colormaps["viridis"].copy()
    cmap.set_bad(color="white")

    for ax, parameter in zip(axes.flat, parameters):
        grid = np.full((len(y_coords), len(x_coords)), np.nan, dtype=float)
        for (die_col, die_row), values in values_by_parameter[parameter].items():
            if values:
                grid[y_index[die_row], x_index[die_col]] = statistics.fmean(values)

        finite_values = grid[np.isfinite(grid)]
        if finite_values.size:
            value_min = float(finite_values.min())
            value_max = float(finite_values.max())
            if value_min == value_max:
                score_grid = np.full_like(grid, 1.0)
            else:
                score_grid = (grid - value_min) / (value_max - value_min)
                if not HIGHER_IS_BETTER.get(parameter, True):
                    score_grid = 1.0 - score_grid
            score_grid[~np.isfinite(grid)] = np.nan
        else:
            score_grid = grid

        image = ax.imshow(
            np.ma.masked_invalid(score_grid),
            cmap=cmap,
            origin="lower",
            aspect="equal",
            vmin=0.0,
            vmax=1.0,
        )

        for (die_col, die_row), values in values_by_parameter[parameter].items():
            if not values:
                continue
            value = statistics.fmean(values)
            score = score_grid[y_index[die_row], x_index[die_col]]
            color_value = cmap(score) if np.isfinite(score) else (0.0, 0.0, 0.0, 1.0)
            luminance = (
                0.2126 * color_value[0] + 0.7152 * color_value[1] + 0.0722 * color_value[2]
            )
            text_color = "black" if luminance > 0.62 else "white"
            ax.text(
                x_index[die_col], y_index[die_row], f"{value:.3g}",
                ha="center", va="center", fontsize=7, color=text_color,
            )

        ax.set_title(parameter, fontsize=10)
        ax.set_xlabel("Die Column")
        ax.set_ylabel("Die Row")
        ax.set_xticks(range(len(x_coords)), x_coords)
        ax.set_yticks(range(len(y_coords)), y_coords)
        ax.set_xticks(np.arange(-0.5, len(x_coords), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(y_coords), 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=1.0)
        ax.tick_params(which="minor", bottom=False, left=False)
        colorbar = fig.colorbar(image, ax=ax, shrink=0.82)
        colorbar.set_ticks([0.0, 1.0])
        colorbar.set_ticklabels(["Bad", "Good"])

    for ax in axes.flat[len(parameters):]:
        ax.set_axis_off()

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return True
