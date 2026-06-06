from __future__ import annotations

import re
import statistics
import xml.etree.ElementTree as ET

import numpy as np
from scipy.signal import argrelextrema

from .xml_parser import (
    attr_any,
    csv_float,
    find_mzm_modulators,
    interpolate_sweeps,
    parse_float_list,
    parse_modulation_sweeps,
)


MIN_VPI_TRACK_R2 = 0.5
MIN_ABS_D_LAMBDA_DV_NM_PER_V = 0.02


def estimate_modulation_fsr(wavelength: np.ndarray, il: np.ndarray) -> float:
    band_mask = (wavelength >= 1535.0) & (wavelength <= 1575.0)
    if np.count_nonzero(band_mask) < 3:
        span = float(wavelength[-1] - wavelength[0])
        band_mask = (
            (wavelength >= wavelength[0] + 0.1 * span)
            & (wavelength <= wavelength[-1] - 0.1 * span)
        )
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
                "null_id": null_id,
                "v": v_arr,
                "wl": wl_arr,
                "coeffs": coeffs,
                "dlambda_dv": float(coeffs[0]),
                "wl_0v": float(np.polyval(coeffs, 0.0)),
                "r2": r2,
            }
        )

    if not track_results:
        return result

    dlambda_values = np.asarray([item["dlambda_dv"] for item in track_results], dtype=float)
    fsr = estimate_modulation_fsr(wavelength, il_matrix[0])
    result.update(
        {
            "track_results": track_results,
            "fsr_nm": fsr,
            "mean_abs_dlambda_dv": float(np.mean(np.abs(dlambda_values))),
            "mean_dlambda_dv": float(np.mean(dlambda_values)),
        }
    )
    return result


def valid_vpi_track(track: dict[str, object]) -> bool:
    try:
        dlambda_dv = abs(float(track.get("dlambda_dv", float("nan"))))
        r2 = float(track.get("r2", float("nan")))
    except (TypeError, ValueError):
        return False
    return (
        np.isfinite(dlambda_dv)
        and np.isfinite(r2)
        and dlambda_dv >= MIN_ABS_D_LAMBDA_DV_NM_PER_V
        and r2 >= MIN_VPI_TRACK_R2
    )


def vpi_value_for_track(track: dict[str, object], fsr: float) -> float | None:
    if not np.isfinite(fsr) or fsr <= 0.0 or not valid_vpi_track(track):
        return None
    dlambda_dv = abs(float(track["dlambda_dv"]))
    return fsr / (2.0 * dlambda_dv)


def vpi_values_from_analysis(analysis: dict[str, object]) -> list[float]:
    fsr = float(analysis.get("fsr_nm", float("nan")))
    if not np.isfinite(fsr) or fsr <= 0.0:
        return []

    return [
        vpi
        for track in analysis.get("track_results", [])
        if (vpi := vpi_value_for_track(track, fsr)) is not None
    ]


def extract_modulation_efficiency(modulator: ET.Element) -> dict[str, object]:
    empty = {
        "modulation_null_count": 0,
        "modulation_fsr_nm": "",
        "modulation_mean_abs_dlambda_dv_nm_per_v": "",
        "modulation_mean_dlambda_dv_nm_per_v": "",
        "modulation_dlambda_dv_by_null_nm_per_v": "",
        "modulation_null_wavelengths_0v_nm": "",
        "modulation_r2_by_null": "",
        "vpi_mean_v": "",
        "vpi_min_v": "",
        "vpi_max_v": "",
        "vpi_by_null_v": "",
    }

    analysis = analyze_modulation_efficiency(modulator)
    track_results = analysis["track_results"]
    if not track_results:
        return empty

    vpi_values = vpi_values_from_analysis(analysis)
    return {
        "modulation_null_count": len(track_results),
        "modulation_fsr_nm": csv_float(float(analysis["fsr_nm"])),
        "modulation_mean_abs_dlambda_dv_nm_per_v": csv_float(float(analysis["mean_abs_dlambda_dv"])),
        "modulation_mean_dlambda_dv_nm_per_v": csv_float(float(analysis["mean_dlambda_dv"])),
        "modulation_dlambda_dv_by_null_nm_per_v": ";".join(
            f"{item['dlambda_dv']:.6f}" for item in track_results
        ),
        "modulation_null_wavelengths_0v_nm": ";".join(
            f"{item['wl_0v']:.4f}" for item in track_results
        ),
        "modulation_r2_by_null": ";".join(f"{item['r2']:.6f}" for item in track_results),
        "vpi_mean_v": csv_float(float(np.mean(vpi_values))) if vpi_values else "",
        "vpi_min_v": csv_float(float(np.min(vpi_values))) if vpi_values else "",
        "vpi_max_v": csv_float(float(np.max(vpi_values))) if vpi_values else "",
        "vpi_by_null_v": ";".join(f"{value:.6f}" for value in vpi_values),
    }


def _local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _normalized_name(value: object) -> str:
    return re.sub(r"[\s_\-./()]+", "", str(value).lower())


def _is_vpi_name(value: object) -> bool:
    return _normalized_name(value) in {"vpi", "vpivoltage", "halfwavevoltage"}


def _is_voltage_name(value: object) -> bool:
    return _normalized_name(value) in {"voltage", "bias", "dcbias"}


def _series_from_text(text: str | None) -> list[float]:
    try:
        return parse_float_list(text)
    except ValueError:
        return []


def _series_from_element(element: ET.Element) -> list[float]:
    for attr_name in ("Value", "MeasuredValue", "Result", "Data"):
        values = _series_from_text(attr_any(element, attr_name))
        if len(values) >= 2:
            return values
    return _series_from_text(element.text)


def _element_name(element: ET.Element) -> str:
    return attr_any(element, "Name", "Parameter", "Symbol", default=_local_name(element.tag))


def _voltage_series_in(element: ET.Element, expected_count: int) -> list[float] | None:
    for candidate in element.iter():
        if not (_is_voltage_name(_local_name(candidate.tag)) or _is_voltage_name(_element_name(candidate))):
            continue
        values = _series_from_element(candidate)
        if len(values) == expected_count:
            return values
    return None


def explicit_vpi_voltage_curves(root: ET.Element) -> list[dict[str, object]]:
    parent_by_child = {child: parent for parent in root.iter() for child in parent}
    curves: list[dict[str, object]] = []

    for element in root.iter():
        if not (_is_vpi_name(_local_name(element.tag)) or _is_vpi_name(_element_name(element))):
            continue

        vpi = _series_from_element(element)
        if len(vpi) < 2:
            continue

        voltage = None
        parent = parent_by_child.get(element)
        if parent is not None:
            voltage = _voltage_series_in(parent, len(vpi))
        if voltage is None:
            voltage = _voltage_series_in(root, len(vpi))
        if voltage is None:
            continue

        curves.append(
            {
                "voltage": np.asarray(voltage, dtype=float),
                "vpi": np.asarray(vpi, dtype=float),
                "label": "XML Vpi",
                "source": "xml",
            }
        )

    return curves


def fitted_vpi_voltage_curves(modulator: ET.Element) -> tuple[list[dict[str, object]], dict[str, object]]:
    analysis = analyze_modulation_efficiency(modulator)
    fsr = float(analysis.get("fsr_nm", float("nan")))
    if not np.isfinite(fsr) or fsr <= 0.0:
        return [], analysis

    curves: list[dict[str, object]] = []
    for track in analysis.get("track_results", []):
        if not valid_vpi_track(track):
            continue

        voltage = np.asarray(track["v"], dtype=float)
        wavelength = np.asarray(track["wl"], dtype=float)
        if voltage.size < 2:
            continue

        order = np.argsort(voltage)
        voltage = voltage[order]
        wavelength = wavelength[order]

        # 2차 다항식 피팅: λ(V) = aV² + bV + c
        # dλ/dV(V) = 2aV + b → 전압 의존성 반영
        deg = min(2, voltage.size - 1)
        coeffs = np.polyfit(voltage, wavelength, deg)

        if deg == 2:
            dlambda_dv = 2.0 * coeffs[0] * voltage + coeffs[1]
        else:
            dlambda_dv = np.full(voltage.shape, coeffs[0], dtype=float)

        with np.errstate(divide="ignore", invalid="ignore"):
            vpi = fsr / (2.0 * np.abs(dlambda_dv))

        keep = (
            np.isfinite(voltage)
            & np.isfinite(vpi)
            & (vpi > 0.0)
            & (np.abs(dlambda_dv) >= MIN_ABS_D_LAMBDA_DV_NM_PER_V)
        )
        if np.count_nonzero(keep) < 2:
            continue

        curves.append(
            {
                "voltage": voltage[keep],
                "vpi": vpi[keep],
                "label": f"Null @{float(track['wl_0v']):.1f} nm",
                "source": "fitted",
            }
        )

    return curves, analysis


def vpi_voltage_curves(root: ET.Element) -> tuple[list[dict[str, object]], str, dict[str, object]]:
    explicit_curves = explicit_vpi_voltage_curves(root)
    if explicit_curves:
        return explicit_curves, "xml", {}

    for modulator in find_mzm_modulators(root):
        curves, analysis = fitted_vpi_voltage_curves(modulator)
        if curves:
            return curves, "fitted", analysis

    return [], "", {}


def _mean_vpi_by_voltage(curves: list[dict[str, object]]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    grouped: dict[float, list[float]] = {}
    for curve in curves:
        voltage = np.asarray(curve["voltage"], dtype=float)
        vpi = np.asarray(curve["vpi"], dtype=float)
        for v, value in zip(voltage, vpi):
            grouped.setdefault(round(float(v), 6), []).append(float(value))

    xs = np.asarray(sorted(grouped), dtype=float)
    means = np.asarray([statistics.fmean(grouped[float(x)]) for x in xs], dtype=float)
    mins = np.asarray([min(grouped[float(x)]) for x in xs], dtype=float)
    maxs = np.asarray([max(grouped[float(x)]) for x in xs], dtype=float)
    return xs, means, mins, maxs


def plot_vpi_voltage_panels(axes, root: ET.Element) -> None:
    ax_curve = axes[0]
    ax_mean = axes[1] if len(axes) > 1 else None
    curves, _source, _analysis = vpi_voltage_curves(root)
    if not curves:
        for ax in axes:
            ax.set_axis_off()
        ax_curve.text(
            0.5,
            0.5,
            "V_pi vs voltage\nNo V_pi array or valid null tracks found",
            transform=ax_curve.transAxes,
            ha="center",
            va="center",
            color="red",
        )
        return

    for curve in curves:
        ax_curve.plot(
            curve["voltage"],
            curve["vpi"],
            marker="o",
            linestyle="-",
            linewidth=1.5,
            markersize=5,
            label=str(curve["label"]),
        )
    ax_curve.set_title("V_pi vs DC bias")
    ax_curve.set_xlabel("DC bias [V]")
    ax_curve.set_ylabel("V_pi [V]")
    ax_curve.grid(True, linestyle="--", alpha=0.35)
    ax_curve.legend(fontsize="x-small", loc="best")

    if ax_mean is None:
        return

    voltage, means, mins, maxs = _mean_vpi_by_voltage(curves)
    yerr = np.vstack([means - mins, maxs - means])
    ax_mean.errorbar(voltage, means, yerr=yerr, marker="o", capsize=4, linewidth=1.5)
    ax_mean.set_xlabel("DC bias [V]")
    ax_mean.set_ylabel("Mean V_pi [V]")
    ax_mean.grid(True, linestyle="--", alpha=0.35)
