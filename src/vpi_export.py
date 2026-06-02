from __future__ import annotations

import argparse
import csv
import re
import statistics
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from src.modulation_efficiency import analyze_modulation_efficiency
    from src.xml_parser import attr_any, find_mzm_modulators
else:
    from .modulation_efficiency import analyze_modulation_efficiency
    from .xml_parser import attr_any, find_mzm_modulators


DATA_DIR = Path("data")
CSV_PATH = Path("res") / "csv" / "vpi_summary.csv"
PNG_PATH = Path("res") / "png" / "vpi_summary.png"

CSV_COLUMNS = [
    "lot",
    "wafer",
    "test_site",
    "die_column",
    "die_row",
    "timestamp",
    "device_name",
    "vpi_v",
    "vpi_source",
    "explicit_vpi_v",
    "fitted_vpi_mean_v",
    "fitted_vpi_min_v",
    "fitted_vpi_max_v",
    "fitted_vpi_by_null_v",
    "modulation_fsr_nm",
    "modulation_mean_abs_dlambda_dv_nm_per_v",
    "modulation_null_count",
    "modulation_r2_by_null",
    "filename",
    "source_file",
]


def _local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _parse_first_float(value: object) -> float | None:
    if value is None:
        return None
    match = re.search(
        r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?",
        str(value),
    )
    if match is None:
        return None
    return float(match.group(0))


def _normalized_name(value: object) -> str:
    return re.sub(r"[\s_\-./()]+", "", str(value).lower())


def _is_vpi_name(value: object) -> bool:
    return _normalized_name(value) in {"vpi", "vpivoltage", "halfwavevoltage"}


def _format_float(value: float | None, digits: int = 6) -> str:
    if value is None or not np.isfinite(value):
        return ""
    return f"{value:.{digits}g}"


def extract_explicit_vpi(element: ET.Element) -> float | None:
    """Return a V_pi value stored directly in XML, if one exists."""
    for elem in element.iter():
        if _is_vpi_name(_local_name(elem.tag)):
            value = _parse_first_float(elem.text)
            if value is not None:
                return value

        name = attr_any(elem, "Name", "Parameter", "Symbol")
        if not _is_vpi_name(name):
            continue

        for attr_name in ("Value", "MeasuredValue", "Result", "Data"):
            value = _parse_first_float(attr_any(elem, attr_name))
            if value is not None:
                return value

        value = _parse_first_float(elem.text)
        if value is not None:
            return value

    return None


def _metadata(root: ET.Element, modulator: ET.Element | None, xml_path: Path) -> dict[str, str]:
    test_site_info = root.find("./TestSiteInfo")
    device_info = modulator.find("./DeviceInfo") if modulator is not None else None
    device_name = (
        (device_info.get("Name") if device_info is not None else None)
        or (modulator.get("Name") if modulator is not None else None)
        or ""
    )

    return {
        "lot": attr_any(test_site_info, "Batch"),
        "wafer": attr_any(test_site_info, "Wafer", default="Unknown"),
        "test_site": attr_any(test_site_info, "TestSite"),
        "die_column": attr_any(test_site_info, "DieColumn", "Diecolumn", default="Unknown"),
        "die_row": attr_any(test_site_info, "DieRow", "Dierow", default="Unknown"),
        "timestamp": xml_path.parent.name,
        "device_name": device_name,
        "filename": xml_path.name,
        "source_file": str(xml_path),
    }


def _fitted_vpi_values(modulator: ET.Element | None) -> tuple[list[float], dict[str, Any]]:
    if modulator is None:
        return [], {}

    analysis = analyze_modulation_efficiency(modulator)
    fsr = float(analysis.get("fsr_nm", float("nan")))
    if not np.isfinite(fsr) or fsr <= 0.0:
        return [], analysis

    vpi_values = []
    for track in analysis.get("track_results", []):
        dlambda_dv = float(track.get("dlambda_dv", float("nan")))
        if np.isfinite(dlambda_dv) and dlambda_dv != 0.0:
            vpi_values.append(fsr / (2.0 * abs(dlambda_dv)))

    return vpi_values, analysis


def extract_vpi_rows(data_dir: Path = DATA_DIR) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    for xml_path in sorted(data_dir.rglob("*LMZ*.xml")):
        try:
            root = ET.parse(xml_path).getroot()
        except ET.ParseError as exc:
            rows.append(
                {
                    "filename": xml_path.name,
                    "source_file": str(xml_path),
                    "vpi_source": f"parse_error: {exc}",
                }
            )
            continue

        modulators = find_mzm_modulators(root)
        if not modulators:
            modulators = [None]

        for modulator in modulators:
            explicit_vpi = (
                extract_explicit_vpi(modulator)
                if modulator is not None
                else None
            )
            if explicit_vpi is None:
                explicit_vpi = extract_explicit_vpi(root)

            fitted_values, analysis = _fitted_vpi_values(modulator)
            fitted_mean = statistics.fmean(fitted_values) if fitted_values else None
            chosen_vpi = explicit_vpi if explicit_vpi is not None else fitted_mean
            source = "xml" if explicit_vpi is not None else ("fitted" if fitted_mean is not None else "")

            track_results = analysis.get("track_results", []) if analysis else []
            row = {
                **_metadata(root, modulator, xml_path),
                "vpi_v": _format_float(chosen_vpi),
                "vpi_source": source,
                "explicit_vpi_v": _format_float(explicit_vpi),
                "fitted_vpi_mean_v": _format_float(fitted_mean),
                "fitted_vpi_min_v": _format_float(min(fitted_values) if fitted_values else None),
                "fitted_vpi_max_v": _format_float(max(fitted_values) if fitted_values else None),
                "fitted_vpi_by_null_v": ";".join(_format_float(value) for value in fitted_values),
                "modulation_fsr_nm": _format_float(float(analysis.get("fsr_nm", float("nan")))) if analysis else "",
                "modulation_mean_abs_dlambda_dv_nm_per_v": (
                    _format_float(float(analysis.get("mean_abs_dlambda_dv", float("nan"))))
                    if analysis else ""
                ),
                "modulation_null_count": len(track_results),
                "modulation_r2_by_null": ";".join(
                    _format_float(float(track.get("r2", float("nan"))))
                    for track in track_results
                ),
            }
            rows.append(row)

    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _row_vpi(row: dict[str, object]) -> float | None:
    value = _parse_first_float(row.get("vpi_v"))
    return value if value is not None and np.isfinite(value) else None


def write_png(path: Path, rows: list[dict[str, object]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    valid_rows = [row for row in rows if _row_vpi(row) is not None]
    path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(max(10, min(24, len(valid_rows) * 0.35)), 6))
    if not valid_rows:
        ax.set_axis_off()
        ax.text(
            0.5,
            0.5,
            "No V_pi data found in LMZ XML files",
            ha="center",
            va="center",
            fontsize=14,
        )
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return

    values = np.asarray([_row_vpi(row) for row in valid_rows], dtype=float)
    sources = [str(row.get("vpi_source", "")) for row in valid_rows]
    colors = ["tab:blue" if source == "xml" else "tab:orange" for source in sources]
    labels = [
        f"{row.get('wafer', '?')} ({row.get('die_column', '?')},{row.get('die_row', '?')})"
        for row in valid_rows
    ]
    positions = np.arange(len(valid_rows))

    ax.scatter(positions, values, c=colors, s=45, edgecolors="black", linewidths=0.5)
    ax.axhline(float(np.mean(values)), color="black", linestyle="--", linewidth=1.0,
               label=f"Mean = {np.mean(values):.3g} V")
    ax.set_title("V_pi summary for LMZ XML files")
    ax.set_xlabel("LMZ measurement")
    ax.set_ylabel("V_pi [V]")
    ax.grid(True, axis="y", linestyle="--", alpha=0.35)
    ax.legend(loc="best")

    if len(labels) <= 35:
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=8)
    else:
        ax.set_xticks([])

    fig.text(
        0.01,
        0.01,
        f"Files with V_pi: {len(valid_rows)} / {len(rows)}   "
        "source: blue=XML, orange=fitted",
        fontsize=9,
        color="dimgray",
    )
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main(
    data_dir: str | Path = DATA_DIR,
    csv_path: str | Path = CSV_PATH,
    png_path: str | Path = PNG_PATH,
) -> list[dict[str, object]]:
    rows = extract_vpi_rows(Path(data_dir))
    write_csv(Path(csv_path), rows)
    write_png(Path(png_path), rows)

    found = sum(1 for row in rows if _row_vpi(row) is not None)
    print(f"Scanned {len(rows)} LMZ modulator entries.")
    print(f"Found V_pi for {found} entries.")
    print(f"CSV output: {csv_path}")
    print(f"PNG output: {png_path}")
    return rows


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract or estimate V_pi for LMZ XML files.")
    parser.add_argument("--data-dir", default=str(DATA_DIR), help="Root folder containing XML data.")
    parser.add_argument("--csv", default=str(CSV_PATH), help="Output CSV path.")
    parser.add_argument("--png", default=str(PNG_PATH), help="Output PNG path.")
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    main(args.data_dir, args.csv, args.png)
