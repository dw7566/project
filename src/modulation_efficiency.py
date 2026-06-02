from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelextrema

from . import config
from .xml_parser import (
    attr_any, csv_float, find_mzm_modulators,
    parse_modulation_sweeps, interpolate_sweeps,
)


def estimate_modulation_fsr(wavelength: np.ndarray, il: np.ndarray) -> float:
    band_mask = (wavelength >= 1535.0) & (wavelength <= 1575.0)
    if np.count_nonzero(band_mask) < 3:
        span = float(wavelength[-1] - wavelength[0])
        band_mask = ((wavelength >= wavelength[0] + 0.1 * span) &
                     (wavelength <= wavelength[-1] - 0.1 * span))
    maxima_idx = argrelextrema(il[band_mask], np.greater, order=40)[0]
    wl_masked = wavelength[band_mask]
    if maxima_idx.size < 2:
        return float("nan")
    spacing = np.diff(wl_masked[maxima_idx])
    spacing = spacing[np.isfinite(spacing) & (spacing > 0)]
    return float(np.mean(spacing)) if spacing.size else float("nan")


def empty_modulation_analysis() -> dict[str, object]:
    return {
        "biases": np.array([], dtype=float),
        "wavelength": np.array([], dtype=float),
        "il_matrix": np.empty((0, 0), dtype=float),
        "track_results": [],
        "fsr_nm": float("nan"),
        "mean_abs_dlambda_dv": float("nan"),
        "mean_dlambda_dv": float("nan"),
    }


def analyze_modulation_efficiency(modulator: ET.Element) -> dict[str, object]:
    result = empty_modulation_analysis()

    interpolated = interpolate_sweeps(parse_modulation_sweeps(modulator))
    if interpolated is None:
        return result
    biases, wavelength, il_matrix = interpolated
    result.update({"biases": biases, "wavelength": wavelength, "il_matrix": il_matrix})

    null_tracks: dict[int, dict[float, float]] = {}
    for index, bias in enumerate(biases):
        minima_idx = argrelextrema(il_matrix[index], np.less, order=50)[0]
        deep_minima = [idx for idx in minima_idx if il_matrix[index][idx] < -30.0]
        for minimum_idx in deep_minima:
            null_wavelength = float(wavelength[minimum_idx])
            matched = False
            for track in null_tracks.values():
                if any(abs(null_wavelength - existing) < 2.0 for existing in track.values()):
                    track[bias] = null_wavelength
                    matched = True
                    break
            if not matched:
                null_tracks[len(null_tracks)] = {bias: null_wavelength}

    full_tracks = {
        null_id: track
        for null_id, track in null_tracks.items()
        if len(track) == biases.size
    }

    track_results = []
    for null_id, track in sorted(full_tracks.items()):
        v_arr = np.asarray(sorted(track), dtype=float)
        wl_arr = np.asarray([track[bias] for bias in v_arr], dtype=float)
        if v_arr.size < 2:
            continue
        coeffs = np.polyfit(v_arr, wl_arr, 1)
        wl_fit = np.polyval(coeffs, v_arr)
        ss_res = float(np.sum((wl_arr - wl_fit) ** 2))
        ss_tot = float(np.sum((wl_arr - np.mean(wl_arr)) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else 0.0
        track_results.append(
            {
                "null_id": null_id, "v": v_arr, "wl": wl_arr,
                "coeffs": coeffs, "dlambda_dv": float(coeffs[0]),
                "wl_0v": float(np.polyval(coeffs, 0.0)), "r2": r2,
            }
        )

    if not track_results:
        return result

    dlambda_values = np.asarray([item["dlambda_dv"] for item in track_results], dtype=float)
    fsr = estimate_modulation_fsr(wavelength, il_matrix[0])
    result.update(
        {
            "track_results": track_results, "fsr_nm": fsr,
            "mean_abs_dlambda_dv": float(np.mean(np.abs(dlambda_values))),
            "mean_dlambda_dv": float(np.mean(dlambda_values)),
        }
    )
    return result


def extract_modulation_efficiency(modulator: ET.Element) -> dict[str, object]:
    empty = {
        "modulation_null_count": 0, "modulation_fsr_nm": "",
        "modulation_mean_abs_dlambda_dv_nm_per_v": "",
        "modulation_mean_dlambda_dv_nm_per_v": "",
        "modulation_dlambda_dv_by_null_nm_per_v": "",
        "modulation_null_wavelengths_0v_nm": "",
        "modulation_r2_by_null": "",
    }

    analysis = analyze_modulation_efficiency(modulator)
    track_results = analysis["track_results"]
    if not track_results:
        return empty

    return {
        "modulation_null_count": len(track_results),
        "modulation_fsr_nm": csv_float(float(analysis["fsr_nm"])),
        "modulation_mean_abs_dlambda_dv_nm_per_v": csv_float(float(analysis["mean_abs_dlambda_dv"])),
        "modulation_mean_dlambda_dv_nm_per_v": csv_float(float(analysis["mean_dlambda_dv"])),
        "modulation_dlambda_dv_by_null_nm_per_v": ";".join(
            f"{item['dlambda_dv']:.6f}" for item in track_results),
        "modulation_null_wavelengths_0v_nm": ";".join(
            f"{item['wl_0v']:.4f}" for item in track_results),
        "modulation_r2_by_null": ";".join(
            f"{item['r2']:.6f}" for item in track_results),
    }


def plot_modulation_efficiency_panels(axes, root: ET.Element) -> None:
    ax_shift, ax_bar, ax_summary = axes
    modulators = find_mzm_modulators(root)
    if not modulators:
        for ax in axes:
            ax.set_axis_off()
        ax_summary.text(0.5, 0.5, "Wavelength modulation\nMZM modulator not found",
                        transform=ax_summary.transAxes, ha="center", va="center", color="red")
        return

    analysis = analyze_modulation_efficiency(modulators[0])
    track_results = analysis["track_results"]
    if not track_results:
        for ax in axes:
            ax.set_axis_off()
        ax_summary.text(0.5, 0.5, "Wavelength modulation\nNo full deep-null tracks",
                        transform=ax_summary.transAxes, ha="center", va="center", color="red")
        return

    markers = ["o", "s", "D", "^", "v", "P", "X", "*"]
    for index, track in enumerate(track_results):
        marker = markers[index % len(markers)]
        ax_shift.scatter(track["v"], track["wl"], s=45, marker=marker, zorder=5,
                         label=f"Null @{float(track['wl_0v']):.1f} nm")
        v_fine = np.linspace(float(np.min(track["v"])) - 0.3,
                             float(np.max(track["v"])) + 0.3, 80)
        ax_shift.plot(v_fine, np.polyval(track["coeffs"], v_fine), "--", linewidth=1.2)

    ax_shift.set_title("Null wavelength vs voltage")
    ax_shift.set_xlabel("DC bias [V]")
    ax_shift.set_ylabel("Null wavelength [nm]")
    ax_shift.legend(fontsize="x-small", loc="best")
    ax_shift.grid(True, ls="--", alpha=0.35)

    labels = [f"@{float(t['wl_0v']):.0f}nm" for t in track_results]
    dlambda_values = [float(t["dlambda_dv"]) for t in track_results]
    positions = np.arange(len(dlambda_values))
    bar_colors = ["steelblue" if v < 0 else "coral" for v in dlambda_values]
    bars = ax_bar.bar(positions, dlambda_values, color=bar_colors, edgecolor="black", alpha=0.85)
    for bar, value in zip(bars, dlambda_values):
        va = "bottom" if value >= 0 else "top"
        offset = 0.003 if value >= 0 else -0.003
        ax_bar.text(bar.get_x() + bar.get_width() / 2, value + offset,
                    f"{value:.4f}", ha="center", va=va, fontsize=8)
    ax_bar.axhline(0.0, color="gray", linewidth=0.8)
    ax_bar.set_xticks(positions)
    ax_bar.set_xticklabels(labels, rotation=25, ha="right")
    ax_bar.set_title("Wavelength modulation efficiency")
    ax_bar.set_xlabel("Null position")
    ax_bar.set_ylabel("dLambda/dV [nm/V]")
    ax_bar.grid(True, axis="y", ls="--", alpha=0.35)

    def _fmt(value: object, spec: str) -> str:
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            return "n/a"
        return format(numeric, spec) if np.isfinite(numeric) else "n/a"

    summary = [
        "Wavelength modulation efficiency",
        f"Deep-null tracks: {len(track_results)}",
        f"FSR: {_fmt(analysis['fsr_nm'], '.3f')} nm",
        f"Mean |dLambda/dV|: {_fmt(analysis['mean_abs_dlambda_dv'], '.4f')} nm/V",
        f"Mean dLambda/dV: {_fmt(analysis['mean_dlambda_dv'], '.4f')} nm/V",
        "", "Per null:",
    ]
    for track in track_results:
        summary.append(
            f"@{float(track['wl_0v']):7.2f} nm  "
            f"{float(track['dlambda_dv']): .4f} nm/V  "
            f"R2={float(track['r2']):.4f}")

    ax_summary.set_axis_off()
    ax_summary.text(0.03, 0.97, "\n".join(summary), transform=ax_summary.transAxes,
                    va="top", ha="left", fontsize=9, family="monospace",
                    bbox=dict(boxstyle="round,pad=0.45", fc="lightyellow",
                              ec="0.5", lw=0.8, alpha=0.95))


def modulation_efficiency_png_path(xml_path: Path, root: ET.Element) -> Path:
    test_site_info = root.find("./TestSiteInfo")
    batch = attr_any(test_site_info, "Batch", default=xml_path.parents[2].name)
    wafer = attr_any(test_site_info, "Wafer", default=xml_path.parent.parent.name)
    timestamp = xml_path.parent.name
    return config.PNG_DIR / config.MODULATION_EFFICIENCY_PNG_DIR / batch / wafer / timestamp / f"{xml_path.stem}.png"


def analyze_modulation_efficiency_figure(xml_path: Path, root: ET.Element | None = None) -> bool:
    if root is None:
        root = ET.parse(xml_path).getroot()

    out_path = modulation_efficiency_png_path(xml_path, root)
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.subplots_adjust(left=0.055, right=0.985, bottom=0.16, top=0.80, wspace=0.32)
    plot_modulation_efficiency_panels(axes, root)

    test_site_info = root.find("./TestSiteInfo")
    batch = attr_any(test_site_info, "Batch", default="?")
    wafer = attr_any(test_site_info, "Wafer", default="?")
    device = attr_any(test_site_info, "TestSite", default="?")
    die = f"({attr_any(test_site_info, 'DieColumn', default='?')},{attr_any(test_site_info, 'DieRow', default='?')})"
    fig.suptitle(f"Wavelength Modulation Efficiency for {wafer} {die} {device}",
                 fontsize=14, fontweight="bold", y=0.97)
    fig.text(0.5, 0.91, f"Batch: {batch}  |  Date: {root.attrib.get('CreationDate', '?')}",
             ha="center", fontsize=10, color="dimgray")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True
