from __future__ import annotations

import csv
import xml.etree.ElementTree as ET
from pathlib import Path

from .config import CSV_COLUMNS
from .xml_parser import (
    attr_any, find_mzm_modulators,
    parse_float_list, nearest_value,
)
from .vpi_analysis import extract_modulation_efficiency


def summarize_xml(xml_path: Path) -> list[dict[str, object]]:
    root = ET.parse(xml_path).getroot()
    test_site_info = root.find("./TestSiteInfo")

    lot = attr_any(test_site_info, "Batch")
    wafer = attr_any(test_site_info, "Wafer")
    test_site = attr_any(test_site_info, "TestSite")
    die_column = attr_any(test_site_info, "DieColumn", "Diecolumn")
    die_row = attr_any(test_site_info, "DieRow", "Dierow")
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

        current_minus_1v = nearest_value(voltage, current, -1.0)
        current_0v = nearest_value(voltage, current, 0.0)
        current_plus_1v = nearest_value(voltage, current, 1.0)
        modulation = extract_modulation_efficiency(modulator)

        for sweep in modulator.findall("./PortCombo/WavelengthSweep"):
            wavelength = parse_float_list(sweep.findtext("./L"))
            il = parse_float_list(sweep.findtext("./IL"))
            if not wavelength or not il:
                continue

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
                    "dc_bias_v": sweep.get("DCBias", ""),
                    "current_at_minus_1v_a": current_minus_1v,
                    "current_at_0v_a": current_0v,
                    "current_at_plus_1v_a": current_plus_1v,
                    "wavelength_start_nm": wavelength[0],
                    "wavelength_stop_nm": wavelength[-1],
                    "extinction_ratio_db": il_max - il_min,
                    "wavelength_at_min_il_nm": wavelength[min_index],
                    "wavelength_at_max_il_nm": wavelength[max_index],
                    **modulation,
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
