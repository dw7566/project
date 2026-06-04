from __future__ import annotations

import re
import statistics
import xml.etree.ElementTree as ET

import numpy as np

from .modulation_efficiency import analyze_modulation_efficiency
from .xml_parser import attr_any, find_mzm_modulators, parse_float_list


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

        keep = np.isfinite(voltage) & np.isfinite(vpi) & (vpi > 0.0)
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