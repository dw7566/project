from __future__ import annotations

import csv
import xml.etree.ElementTree as ET
from pathlib import Path

from .config import CSV_COLUMNS
from .xml_parser import (
    attr_any, die_coordinates, find_mzm_modulators,
    parse_float_list, nearest_value,
)
from .vpi_analysis import extract_modulation_efficiency, vpi_by_bias_from_modulator


def _bias_key(value: object) -> float | None:
    try:
        return round(float(str(value)), 6)
    except (TypeError, ValueError):
        return None


def summarize_xml(xml_path: Path) -> list[dict[str, object]]:
    root = ET.parse(xml_path).getroot()
    test_site_info = root.find("./TestSiteInfo")

    lot = attr_any(test_site_info, "Batch")
    wafer = attr_any(test_site_info, "Wafer")
    test_site = attr_any(test_site_info, "TestSite")
    die_column, die_row = die_coordinates(xml_path, test_site_info)
    timestamp = xml_path.parent.name

    rows: list[dict[str, object]] = []
    for modulator in find_mzm_modulators(root):
        device_info = modulator.find("./DeviceInfo")
        device_name = (
            (device_info.get("Name") if device_info is not None else None)
            or modulator.get("Name") or ""
        )

        port_combo = modulator.find("./PortCombo")
        voltage = parse_float_list(
            port_combo.findtext("./IVMeasurement/Voltage") if port_combo is not None else None
        )
        current = parse_float_list(
            port_combo.findtext("./IVMeasurement/Current") if port_combo is not None else None
        )

        modulation = extract_modulation_efficiency(modulator)
        vpi_by_bias = vpi_by_bias_from_modulator(modulator)

        for sweep in modulator.findall("./PortCombo/WavelengthSweep"):
            dc_bias = sweep.get("DCBias", "")
            wavelength = parse_float_list(sweep.findtext("./L"))
            il = parse_float_list(sweep.findtext("./IL"))
            if not wavelength or not il:
                continue

            # 해당 DC bias에서의 전류
            try:
                bias_float = float(dc_bias)
            except (TypeError, ValueError):
                bias_float = 0.0
            current_at_bias = nearest_value(voltage, current, bias_float)

            count = min(len(wavelength), len(il))
            wavelength = wavelength[:count]
            il = il[:count]
            min_index = min(range(count), key=il.__getitem__)
            max_index = max(range(count), key=il.__getitem__)
            il_min = il[min_index]
            il_max = il[max_index]

            rows.append(
                {
                    "lot": lot, "wafer": wafer, "test_site": test_site,
                    "die_column": die_column, "die_row": die_row,
                    "timestamp": timestamp, "device_name": device_name,
                    "dc_bias_v": dc_bias,
                    "current_a": current_at_bias,
                    "wavelength_start_nm": wavelength[0],
                    "wavelength_stop_nm": wavelength[-1],
                    "extinction_ratio_db": il_max - il_min,
                    **modulation,
                    "vpi_at_dc_bias_v": vpi_by_bias.get(_bias_key(dc_bias), ""),
                    "source_file": str(xml_path),
                }
            )

    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

def rows_to_dataframe(rows: list[dict[str, object]]) -> pd.DataFrame:
    """Return rows in the exact column order used by write_csv()."""
    return pd.DataFrame(rows).reindex(columns=CSV_COLUMNS)