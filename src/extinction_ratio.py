from __future__ import annotations

import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelextrema

from .xml_parser import (
    find_mzm_modulators, parse_modulation_sweeps, interpolate_sweeps,
)


def analyze_extinction_ratio(modulator: ET.Element) -> dict[str, object]:
    """Deep null 기반 소광비(ER) 추출. ER < 10 dB인 가짜 null은 제외."""
    empty: dict[str, object] = {
        "er_results": [],
        "biases": np.array([], dtype=float),
    }

    interpolated = interpolate_sweeps(parse_modulation_sweeps(modulator))
    if interpolated is None:
        return empty
    biases, wavelength, il_matrix = interpolated

    MIN_ER_DB = 10
    order = 40
    er_results: list[dict[str, object]] = []

    for i, bias in enumerate(biases):
        il = il_matrix[i]
        maxima_idx = argrelextrema(il, np.greater, order=order)[0]
        minima_idx = argrelextrema(il, np.less, order=order)[0]

        er_list: list[float] = []
        pair_info: list[dict[str, float]] = []

        for mi in minima_idx:
            null_wl = float(wavelength[mi])
            null_il = float(il[mi])

            left_peaks = maxima_idx[maxima_idx < mi]
            right_peaks = maxima_idx[maxima_idx > mi]

            if len(left_peaks) > 0 and len(right_peaks) > 0:
                lp = left_peaks[-1]
                rp = right_peaks[0]
                peak_idx = lp if il[lp] >= il[rp] else rp
                peak_wl = float(wavelength[peak_idx])
                peak_il = float(il[peak_idx])
                er = peak_il - null_il

                if er < MIN_ER_DB:
                    continue

                er_list.append(er)
                pair_info.append({
                    "null_wl": null_wl, "null_il": null_il,
                    "peak_wl": peak_wl, "peak_il": peak_il,
                    "er": er,
                })

        er_results.append({
            "bias": float(bias),
            "er_list": er_list,
            "er_mean": float(np.mean(er_list)) if er_list else 0.0,
            "er_max": float(np.max(er_list)) if er_list else 0.0,
            "er_min": float(np.min(er_list)) if er_list else 0.0,
            "pairs": pair_info,
        })

    return {"er_results": er_results, "biases": biases}


def plot_extinction_ratio_panels(axes, root: ET.Element) -> None:
    ax_bias, ax_fringe, ax_summary = axes

    modulators = find_mzm_modulators(root)
    if not modulators:
        for ax in axes:
            ax.set_axis_off()
        ax_summary.text(0.5, 0.5, "Extinction ratio\nMZM modulator not found",
                        transform=ax_summary.transAxes, ha="center", va="center", color="red")
        return

    analysis = analyze_extinction_ratio(modulators[0])
    er_results = analysis["er_results"]
    biases = analysis["biases"]

    if not er_results or all(len(r["er_list"]) == 0 for r in er_results):
        for ax in axes:
            ax.set_axis_off()
        ax_summary.text(0.5, 0.5, "Extinction ratio\nNo valid fringes found",
                        transform=ax_summary.transAxes, ha="center", va="center", color="red")
        return

    bias_arr = np.array([r["bias"] for r in er_results])
    er_means = [r["er_mean"] for r in er_results]
    er_maxs = [r["er_max"] for r in er_results]
    er_mins = [r["er_min"] for r in er_results]

    ax_bias.fill_between(bias_arr, er_mins, er_maxs, alpha=0.2,
                         color="steelblue", label="Min\u2013Max range")
    ax_bias.plot(bias_arr, er_means, "bo-", ms=7, lw=2, label="Mean ER")
    ax_bias.plot(bias_arr, er_maxs, "g^--", ms=5, lw=1, alpha=0.7, label="Max ER")
    ax_bias.plot(bias_arr, er_mins, "rv--", ms=5, lw=1, alpha=0.7, label="Min ER")
    for b, m in zip(bias_arr, er_means):
        ax_bias.annotate(f"{m:.1f}", (b, m), textcoords="offset points",
                         xytext=(0, 10), ha="center", fontsize=8)
    ax_bias.set_title("ER vs DC bias")
    ax_bias.set_xlabel("DC bias [V]")
    ax_bias.set_ylabel("Extinction ratio [dB]")
    ax_bias.legend(fontsize="x-small", loc="best")
    ax_bias.grid(True, ls="--", alpha=0.35)

    idx_0v = int(np.argmin(np.abs(bias_arr - 0.0)))
    res_0v = er_results[idx_0v]

    if res_0v["pairs"]:
        labels = [f"@{p['null_wl']:.0f}nm" for p in res_0v["pairs"]]
        er_vals = [p["er"] for p in res_0v["pairs"]]
        positions = np.arange(len(er_vals))
        bar_colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(er_vals)))
        bars = ax_fringe.bar(positions, er_vals, color=bar_colors,
                             edgecolor="black", alpha=0.85)
        for bar, val in zip(bars, er_vals):
            ax_fringe.text(bar.get_x() + bar.get_width() / 2,
                           bar.get_height() + 0.3, f"{val:.1f}",
                           ha="center", va="bottom", fontsize=8)
        ax_fringe.axhline(np.mean(er_vals), color="red", ls="--", lw=1,
                          label=f"Mean={np.mean(er_vals):.1f} dB")
        ax_fringe.set_xticks(positions)
        ax_fringe.set_xticklabels(labels, rotation=25, ha="right", fontsize=8)
    ax_fringe.set_title("ER by fringe position (@ 0 V)")
    ax_fringe.set_xlabel("Null position")
    ax_fringe.set_ylabel("Extinction ratio [dB]")
    ax_fringe.legend(fontsize="x-small")
    ax_fringe.grid(True, axis="y", ls="--", alpha=0.35)

    def _fmt(value: float, spec: str) -> str:
        return format(value, spec) if np.isfinite(value) else "n/a"

    summary_lines = ["Extinction ratio analysis", "", "Per bias:"]
    for r in er_results:
        n = len(r["er_list"])
        summary_lines.append(
            f"  V={r['bias']:+5.1f}V  ER={_fmt(r['er_mean'], '.1f'):>6s} dB  ({n} fringes)")
    summary_lines += [
        "", "Overall:",
        f"  ER range: {_fmt(min(er_mins), '.1f')} ~ {_fmt(max(er_maxs), '.1f')} dB",
        f"  Best ER:  {_fmt(max(er_maxs), '.1f')} dB",
        f"  Bias range: {bias_arr[0]:.1f} ~ {bias_arr[-1]:.1f} V",
        f"  Fringe count (0 V): {len(res_0v['pairs'])}",
    ]

    ax_summary.set_axis_off()
    ax_summary.text(0.03, 0.97, "\n".join(summary_lines),
                    transform=ax_summary.transAxes, va="top", ha="left",
                    fontsize=9, family="monospace",
                    bbox=dict(boxstyle="round,pad=0.45", fc="lightyellow",
                              ec="0.5", lw=0.8, alpha=0.95))
