from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np


def parse_float_list(text: str | None) -> list[float]:
    if not text:
        return []
    return [float(value) for value in text.split(",") if value.strip()]


def parse_float_array(text: str | None) -> np.ndarray:
    return np.asarray(parse_float_list(text), dtype=float)


def attr_any(element: ET.Element | None, *names: str, default: str = "") -> str:
    if element is None:
        return default
    lower_map = {key.lower(): value for key, value in element.attrib.items()}
    for name in names:
        value = lower_map.get(name.lower())
        if value is not None:
            return value
    return default


def nearest_value(xs: list[float], ys: list[float], target: float) -> float | None:
    if not xs or not ys:
        return None
    limit = min(len(xs), len(ys))
    best_index = min(range(limit), key=lambda idx: abs(xs[idx] - target))
    return abs(ys[best_index])


def r2_score(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    if mask.sum() < 2:
        return float("nan")
    y_true = y_true[mask]
    y_pred = y_pred[mask]
    ss_res = float(np.sum((y_true - y_pred) ** 2))
    ss_tot = float(np.sum((y_true - np.mean(y_true)) ** 2))
    if ss_tot == 0.0:
        return 1.0 if ss_res == 0.0 else float("nan")
    return 1.0 - ss_res / ss_tot


def csv_float(value: float, digits: int = 6) -> str:
    if not np.isfinite(value):
        return ""
    return f"{value:.{digits}g}"


def find_mzm_modulators(root: ET.Element) -> list[ET.Element]:
    modulators = []
    for modulator in root.findall(".//Modulator"):
        device_info = modulator.find("./DeviceInfo")
        names = [
            (modulator.get("Name") or "").upper(),
            (device_info.get("Name") if device_info is not None else "").upper(),
        ]
        if any(name.startswith("MZM") for name in names):
            modulators.append(modulator)
    return modulators


def load_xml(file_path: Path) -> tuple[ET.Element, list[dict[str, object]], dict[str, np.ndarray]]:
    root = ET.parse(file_path).getroot()

    sweeps: list[dict[str, object]] = []
    for sweep in root.findall(".//WavelengthSweep"):
        wavelength = parse_float_array(sweep.findtext("./L"))
        il = parse_float_array(sweep.findtext("./IL"))
        count = min(wavelength.size, il.size)
        if count == 0:
            continue
        sweeps.append(
            {"L": wavelength[:count], "IL": il[:count],
             "Bias": sweep.attrib.get("DCBias", "0.0")}
        )

    iv = {"V": np.array([], dtype=float), "I": np.array([], dtype=float)}
    iv_node = root.find(".//IVMeasurement")
    if iv_node is not None:
        voltage = parse_float_array(iv_node.findtext("./Voltage"))
        current = parse_float_array(iv_node.findtext("./Current"))
        count = min(voltage.size, current.size)
        iv["V"] = voltage[:count]
        iv["I"] = current[:count]

    return root, sweeps, iv


def pick_sweep(sweeps: list[dict[str, object]], bias: str) -> dict[str, object] | None:
    for sweep in sweeps:
        if sweep["Bias"] == bias:
            return sweep
    return None


def sweep_label(sweep: dict[str, object], index: int, total: int) -> str:
    bias = sweep["Bias"]
    if index == total - 1:
        return f"Reference ({bias}V)"
    return f"Bias {bias}V"


def parse_modulation_sweeps(modulator: ET.Element) -> list[tuple[float, np.ndarray, np.ndarray]]:
    sweeps = []
    for sweep in modulator.findall("./PortCombo/WavelengthSweep"):
        try:
            bias = float(sweep.get("DCBias", "nan"))
        except ValueError:
            continue
        wavelength = parse_float_array(sweep.findtext("./L"))
        il = parse_float_array(sweep.findtext("./IL"))
        count = min(wavelength.size, il.size)
        if count == 0 or not np.isfinite(bias):
            continue
        wavelength = wavelength[:count]
        il = il[:count]
        order = np.argsort(wavelength)
        sweeps.append((bias, wavelength[order], il[order]))
    return sweeps


def interpolate_sweeps(
    sweeps: list[tuple[float, np.ndarray, np.ndarray]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    by_bias: dict[float, tuple[np.ndarray, np.ndarray]] = {}
    for bias, wavelength, il in sweeps:
        by_bias.setdefault(bias, (wavelength, il))
    if len(by_bias) < 2:
        return None

    items = sorted(by_bias.items())
    low = max(float(wavelength[0]) for _, (wavelength, _) in items)
    high = min(float(wavelength[-1]) for _, (wavelength, _) in items)
    if not low < high:
        return None

    first_wavelength = items[0][1][0]
    ref_wl = first_wavelength[(first_wavelength >= low) & (first_wavelength <= high)]
    if ref_wl.size < 3:
        return None

    il_matrix = []
    for _, (wavelength, il) in items:
        il_matrix.append(np.interp(ref_wl, wavelength, il))
    biases = np.asarray([bias for bias, _ in items], dtype=float)
    return biases, ref_wl, np.asarray(il_matrix, dtype=float)
