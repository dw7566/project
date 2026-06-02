from __future__ import annotations

import re
import statistics
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from . import config
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
        if voltage.size == 2:
            slope = (wavelength[1] - wavelength[0]) / (voltage[1] - voltage[0])
            dlambda_dv = np.full(voltage.shape, slope, dtype=float)
        else:
            dlambda_dv = np.gradient(wavelength, voltage)

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


def _format(value: object, spec: str) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return "n/a"
    return format(numeric, spec) if np.isfinite(numeric) else "n/a"


def _all_vpi(curves: list[dict[str, object]]) -> np.ndarray:
    if not curves:
        return np.array([], dtype=float)
    return np.concatenate([np.asarray(curve["vpi"], dtype=float) for curve in curves])


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
    ax_curve, ax_mean, ax_summary = axes
    curves, source, analysis = vpi_voltage_curves(root)
    if not curves:
        for ax in axes:
            ax.set_axis_off()
        ax_summary.text(
            0.5,
            0.5,
            "V_pi vs voltage\nNo V_pi array or valid null tracks found",
            transform=ax_summary.transAxes,
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

    voltage, means, mins, maxs = _mean_vpi_by_voltage(curves)
    yerr = np.vstack([means - mins, maxs - means])
    ax_mean.errorbar(voltage, means, yerr=yerr, marker="o", capsize=4, linewidth=1.5)
    ax_mean.set_title("Mean V_pi by voltage")
    ax_mean.set_xlabel("DC bias [V]")
    ax_mean.set_ylabel("Mean V_pi [V]")
    ax_mean.grid(True, linestyle="--", alpha=0.35)

    values = _all_vpi(curves)
    summary = [
        "V_pi vs voltage",
        f"Source: {source}",
        f"Curves: {len(curves)}",
        f"Voltage range: {_format(np.min(voltage), '.2f')} ~ {_format(np.max(voltage), '.2f')} V",
        f"Mean V_pi: {_format(np.mean(values), '.3f')} V",
        f"V_pi range: {_format(np.min(values), '.3f')} ~ {_format(np.max(values), '.3f')} V",
    ]
    if source == "fitted":
        summary.extend(
            [
                f"FSR: {_format(analysis.get('fsr_nm'), '.3f')} nm",
                "Formula: FSR / (2*abs(dLambda/dV))",
            ]
        )
    summary.extend(["", "Per curve:"])
    for curve in curves:
        curve_values = np.asarray(curve["vpi"], dtype=float)
        summary.append(f"{curve['label']}: mean={np.mean(curve_values):.3f} V")

    ax_summary.set_axis_off()
    ax_summary.text(
        0.03,
        0.97,
        "\n".join(summary),
        transform=ax_summary.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.45", fc="lightyellow", ec="0.5", lw=0.8, alpha=0.95),
    )


def vpi_voltage_png_path(xml_path: Path, root: ET.Element) -> Path:
    test_site_info = root.find("./TestSiteInfo")
    batch = attr_any(test_site_info, "Batch", default=xml_path.parents[2].name)
    wafer = attr_any(test_site_info, "Wafer", default=xml_path.parent.parent.name)
    timestamp = xml_path.parent.name
    return config.PNG_DIR / config.VPI_VOLTAGE_PNG_DIR / batch / wafer / timestamp / f"{xml_path.stem}.png"


def analyze_vpi_voltage_figure(xml_path: Path, root: ET.Element | None = None) -> bool:
    if root is None:
        root = ET.parse(xml_path).getroot()

    out_path = vpi_voltage_png_path(xml_path, root)
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.subplots_adjust(left=0.055, right=0.985, bottom=0.16, top=0.80, wspace=0.32)
    plot_vpi_voltage_panels(axes, root)

    test_site_info = root.find("./TestSiteInfo")
    batch = attr_any(test_site_info, "Batch", default="?")
    wafer = attr_any(test_site_info, "Wafer", default="?")
    device = attr_any(test_site_info, "TestSite", default="?")
    die = f"({attr_any(test_site_info, 'DieColumn', default='?')},{attr_any(test_site_info, 'DieRow', default='?')})"
    fig.suptitle(f"V_pi vs Voltage for {wafer} {die} {device}",
                 fontsize=14, fontweight="bold", y=0.97)
    fig.text(0.5, 0.91, f"Batch: {batch}  |  Date: {root.attrib.get('CreationDate', '?')}",
             ha="center", fontsize=10, color="dimgray")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True
