from __future__ import annotations

import csv
import statistics
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema, find_peaks


DATA_DIR = Path("data")
CSV_DIR = Path("res") / "csv"
PNG_DIR = Path("res") / "png"

THERMAL_VOLTAGE = 0.02585
MOD_BIAS = "-1.0"

plt.rcParams["font.family"] = ["DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


CSV_COLUMNS = [
    "lot", "wafer", "test_site", "die_column", "die_row", "timestamp",
    "device_name", "dc_bias_v", "current_at_minus_1v_a", "current_at_0v_a",
    "current_at_plus_1v_a", "wavelength_start_nm", "wavelength_stop_nm",
    "point_count", "il_min_db", "il_max_db", "il_mean_db",
    "extinction_ratio_db", "wavelength_at_min_il_nm", "wavelength_at_max_il_nm",
    "modulation_null_count", "modulation_fsr_nm",
    "modulation_mean_abs_dlambda_dv_nm_per_v",
    "modulation_mean_dlambda_dv_nm_per_v",
    "modulation_dlambda_dv_by_null_nm_per_v",
    "modulation_null_wavelengths_0v_nm", "modulation_r2_by_null",
    "source_file",
]


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


def mzi_model(wavelength: np.ndarray, A: float, B: float, wl0: float, FSR: float, phi: float, slope: float) -> np.ndarray:
    x = wavelength - wavelength.mean()
    return A + slope * x + B * np.cos(np.pi * (wavelength - wl0) / FSR + phi) ** 2


def measure_fsr(wavelength: np.ndarray, transmission_db: np.ndarray, fallback: float) -> float:
    step = float(np.median(np.diff(wavelength))) if wavelength.size > 2 else 0.01
    min_distance = max(1, int(round(2.0 / max(abs(step), 1e-6))))
    for prominence in (1.0, 0.5, 0.25):
        notches, _ = find_peaks(-transmission_db, prominence=prominence, distance=min_distance)
        if len(notches) >= 2:
            spacing = np.diff(wavelength[notches])
            spacing = spacing[np.isfinite(spacing) & (spacing > 0)]
            if spacing.size:
                return float(np.median(spacing))
    return fallback


def crest_points(wavelength: np.ndarray, values_db: np.ndarray, fsr: float) -> tuple[np.ndarray, np.ndarray]:
    """Locate the fringe crests (transmission maxima) the envelope should pass through.

    Three guards keep only genuine crests:
      * ``prominence`` (scaled to fringe depth) rejects shallow bumps -- e.g. the
        small upturn where a band edge falls mid-notch, which is otherwise picked
        up as a spurious low "crest" near 1580 nm and drags the envelope down,
        pushing the real crest near it several dB above 0.
      * band-edge anchors recover crests that ``find_peaks`` drops at the very ends.
      * a robust pass removes any candidate sitting well below the crest trend.
    """
    step = float(np.median(np.diff(wavelength))) if wavelength.size > 2 else 0.01
    distance = max(1, int(round(0.6 * fsr / max(abs(step), 1e-6))))
    span = float(np.percentile(values_db, 95) - np.percentile(values_db, 5))
    prominence = float(np.clip(0.15 * span, 2.0, 8.0))
    peaks, _ = find_peaks(values_db, distance=distance, prominence=prominence)
    xs = list(wavelength[peaks])
    ys = list(values_db[peaks])
    for lo, hi in ((wavelength.min(), wavelength.min() + fsr),
                   (wavelength.max() - fsr, wavelength.max())):
        mask = (wavelength >= lo) & (wavelength <= hi)
        if np.any(mask):
            idx = int(np.argmax(values_db[mask]))
            xs.append(float(wavelength[mask][idx]))
            ys.append(float(values_db[mask][idx]))
    xs_arr = np.asarray(xs)
    ys_arr = np.asarray(ys)
    order = np.argsort(xs_arr)
    xs_arr, ys_arr = xs_arr[order], ys_arr[order]
    _, unique_idx = np.unique(np.round(xs_arr, 3), return_index=True)
    xs_arr, ys_arr = xs_arr[unique_idx], ys_arr[unique_idx]
    # robust rejection: drop crests > 2.5 dB below the linear crest trend
    for _ in range(4):
        if xs_arr.size < 4:
            break
        trend = np.poly1d(np.polyfit(xs_arr, ys_arr, 1))(xs_arr)
        keep = (ys_arr - trend) > -2.5
        if keep.all():
            break
        xs_arr, ys_arr = xs_arr[keep], ys_arr[keep]
    return xs_arr, ys_arr


def top_envelope(wavelength: np.ndarray, values_db: np.ndarray, fsr: float, degree: int = 2) -> np.ndarray:
    """Upper (peak) envelope: monotone PCHIP interpolation through the detected
    fringe crests, with the ends held flat so there is no runaway extrapolation.

    A global polynomial cannot follow an asymmetric crest trend (it is forced
    symmetric), so the crest at one band edge floats above 0 after subtraction.
    A PCHIP curve passes exactly through every crest and does not overshoot
    between them, so all crests land on the 0 dB baseline. Falls back to a
    polynomial only when too few crests are found for interpolation.
    """
    xs, ys = crest_points(wavelength, values_db, fsr)
    if xs.size >= 2:
        envelope = PchipInterpolator(xs, ys, extrapolate=True)(wavelength)
        envelope[wavelength < xs[0]] = ys[0]
        envelope[wavelength > xs[-1]] = ys[-1]
        return envelope
    # fallback: one window per FSR so each window max is a genuine crest
    n_windows = int(np.clip(round((wavelength.max() - wavelength.min()) / fsr), 3, 12))
    edges = np.linspace(float(wavelength.min()), float(wavelength.max()), n_windows + 1)
    xs_l: list[float] = []
    ys_l: list[float] = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (wavelength >= lo) & (wavelength <= hi)
        if np.any(mask):
            idx = int(np.argmax(values_db[mask]))
            xs_l.append(float(wavelength[mask][idx]))
            ys_l.append(float(values_db[mask][idx]))
    if len(xs_l) < 2:
        return np.zeros_like(wavelength)
    degree = min(degree, len(xs_l) - 1)
    return np.poly1d(np.polyfit(xs_l, ys_l, degree))(wavelength)


def flatten_to_envelope(wavelength: np.ndarray, values_db: np.ndarray,
                        fsr: float | None = None, degree: int = 2) -> np.ndarray:
    """Second-stage flattening: subtract the residual device envelope so the
    fringe maxima land on the 0 dB baseline.

    Reference subtraction (IL - cubic_fit(reference)) only removes the *reference*
    waveguide envelope. Devices whose grating-coupler / insertion-loss envelope
    differs from the reference port -- notably the O-band LMZO devices -- keep a
    dome-shaped residual after step 1, so their fringe peaks droop several dB at
    the band edges instead of sitting flat at 0 dB. This step removes that dome.
    """
    if wavelength.size < 4:
        return values_db
    band = float(wavelength.max() - wavelength.min())
    if not fsr or fsr <= 0:
        fsr = band / 6.0
    return values_db - top_envelope(wavelength, values_db, fsr, degree)


def diode_eq(voltage: np.ndarray, Is: float, n: float) -> np.ndarray:
    exponent = np.clip(voltage / (n * THERMAL_VOLTAGE), -700, 700)
    return Is * (np.exp(exponent) - 1.0)


def device_fsr_fallback(root: ET.Element) -> float:
    site = attr_any(root.find(".//TestSiteInfo"), "TestSite").upper()
    if "LMZO" in site:
        return 10.0
    if "LMZC" in site:
        return 14.0
    return 12.0


def sweep_label(sweep: dict[str, object], index: int, total: int) -> str:
    bias = sweep["Bias"]
    if index == total - 1:
        return f"Reference ({bias}V)"
    return f"Bias {bias}V"


def csv_float(value: float, digits: int = 6) -> str:
    if not np.isfinite(value):
        return ""
    return f"{value:.{digits}g}"


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


def extract_modulation_efficiency(modulator: ET.Element) -> dict[str, object]:
    empty = {
        "modulation_null_count": 0,
        "modulation_fsr_nm": "",
        "modulation_mean_abs_dlambda_dv_nm_per_v": "",
        "modulation_mean_dlambda_dv_nm_per_v": "",
        "modulation_dlambda_dv_by_null_nm_per_v": "",
        "modulation_null_wavelengths_0v_nm": "",
        "modulation_r2_by_null": "",
    }

    interpolated = interpolate_sweeps(parse_modulation_sweeps(modulator))
    if interpolated is None:
        return empty
    biases, wavelength, il_matrix = interpolated

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
                "dlambda_dv": float(coeffs[0]),
                "wl_0v": float(np.polyval(coeffs, 0.0)),
                "r2": r2,
            }
        )

    if not track_results:
        return empty

    dlambda_values = np.asarray([item["dlambda_dv"] for item in track_results], dtype=float)
    fsr = estimate_modulation_fsr(wavelength, il_matrix[0])
    return {
        "modulation_null_count": len(track_results),
        "modulation_fsr_nm": csv_float(fsr),
        "modulation_mean_abs_dlambda_dv_nm_per_v": csv_float(float(np.mean(np.abs(dlambda_values)))),
        "modulation_mean_dlambda_dv_nm_per_v": csv_float(float(np.mean(dlambda_values))),
        "modulation_dlambda_dv_by_null_nm_per_v": ";".join(
            f"{item['dlambda_dv']:.6f}" for item in track_results
        ),
        "modulation_null_wavelengths_0v_nm": ";".join(
            f"{item['wl_0v']:.4f}" for item in track_results
        ),
        "modulation_r2_by_null": ";".join(
            f"{item['r2']:.6f}" for item in track_results
        ),
    }


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
                    "point_count": count, "il_min_db": il_min, "il_max_db": il_max,
                    "il_mean_db": statistics.fmean(il),
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
        writer = csv.DictWriter(csv_file, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def mean(values: list[float]) -> float:
    return statistics.fmean(values) if values else 0.0


def plot_wafer_summary(wafer: str, rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    best_by_die: dict[tuple[int, int], list[float]] = defaultdict(list)
    er_by_bias: dict[float, list[float]] = defaultdict(list)

    for row in rows:
        try:
            die_col = int(str(row["die_column"]))
            die_row = int(str(row["die_row"]))
            er = float(row["extinction_ratio_db"])
            bias = float(str(row["dc_bias_v"]))
        except (TypeError, ValueError):
            continue
        best_by_die[(die_col, die_row)].append(er)
        er_by_bias[bias].append(er)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    fig.suptitle(f"{wafer} MZM summary", fontsize=14)

    ax = axes[0]
    if best_by_die:
        keys = sorted(best_by_die)
        xs = [key[0] for key in keys]
        ys = [key[1] for key in keys]
        values = [max(best_by_die[key]) for key in keys]
        scatter = ax.scatter(xs, ys, c=values, cmap="viridis", marker="s",
                             s=650, edgecolors="black", linewidths=0.8)
        for x, y, value in zip(xs, ys, values):
            ax.text(x, y, f"{value:.1f}", ha="center", va="center", fontsize=8, color="white")
        fig.colorbar(scatter, ax=ax, label="Best extinction ratio [dB]")
    ax.set_title("Best ER by die")
    ax.set_xlabel("Die column")
    ax.set_ylabel("Die row")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.25)

    ax = axes[1]
    if er_by_bias:
        biases = sorted(er_by_bias)
        means = [mean(er_by_bias[bias]) for bias in biases]
        ax.plot(biases, means, marker="o", linewidth=2)
    ax.set_title("Mean ER by DC bias")
    ax.set_xlabel("DC bias [V]")
    ax.set_ylabel("Extinction ratio [dB]")
    ax.grid(True, alpha=0.3)

    fig.savefig(path, dpi=180)
    plt.close(fig)


def clean_legacy_outputs() -> None:
    """Remove files from the old flat output layout before writing the dated layout."""
    if PNG_DIR.exists():
        for wafer_dir in PNG_DIR.iterdir():
            if not wafer_dir.is_dir():
                continue
            for png_path in wafer_dir.glob("*LMZ*.png"):
                png_path.unlink()

    if CSV_DIR.exists():
        for csv_path in CSV_DIR.glob("D*_mzm_summary.csv"):
            csv_path.unlink()


def measurement_folders(data_dir: Path) -> list[tuple[str, str]]:
    folders = {
        (xml_path.parent.parent.name, xml_path.parent.name)
        for xml_path in data_dir.rglob("*.xml")
        if xml_path.parent.parent != data_dir
    }
    return sorted(folders)


def mirror_result_folders(data_dir: Path) -> list[tuple[str, str]]:
    folders = measurement_folders(data_dir)
    for wafer, timestamp in folders:
        (PNG_DIR / wafer / timestamp).mkdir(parents=True, exist_ok=True)
        (CSV_DIR / wafer / timestamp).mkdir(parents=True, exist_ok=True)
    return folders


def unique_png_path(wafer: str, timestamp: str, xml_path: Path,
                    used_names_by_folder: dict[tuple[str, str], set[str]]) -> Path:
    output_dir = PNG_DIR / wafer / timestamp
    base_name = f"{xml_path.stem}.png"
    used_names = used_names_by_folder[(wafer, timestamp)]

    if base_name not in used_names:
        used_names.add(base_name)
        return output_dir / base_name

    index = 2
    while True:
        fallback_name = f"{xml_path.stem}__{index}.png"
        if fallback_name not in used_names:
            used_names.add(fallback_name)
            return output_dir / fallback_name
        index += 1


def analyze_figure(xml_path: Path, out_path: Path) -> bool:
    root, sweeps, iv = load_xml(xml_path)
    if not sweeps:
        return False

    ref = sweeps[-1]
    ref_l = ref["L"]
    ref_il = ref["IL"]
    assert isinstance(ref_l, np.ndarray)
    assert isinstance(ref_il, np.ndarray)
    poly_func = np.poly1d(np.polyfit(ref_l, ref_il, 3))
    fsr_fallback = device_fsr_fallback(root)

    fig, axes = plt.subplots(2, 3, figsize=(20, 11))
    plt.subplots_adjust(hspace=0.35, wspace=0.3)

    ax1 = axes[0, 0]
    for index, sweep in enumerate(sweeps):
        wavelength = sweep["L"]
        il = sweep["IL"]
        assert isinstance(wavelength, np.ndarray)
        assert isinstance(il, np.ndarray)
        ax1.plot(wavelength, il, label=sweep_label(sweep, index, len(sweeps)), linewidth=1.0)
    ax1.set_title("Transmission spectra - as measured")
    ax1.set_xlabel("Wavelength [nm]")
    ax1.set_ylabel("Measured transmission [dB]")
    ax1.legend(ncol=2, fontsize="small", loc="lower center")
    ax1.grid(True, ls="--", alpha=0.5)

    ax2 = axes[0, 1]
    ref_fit = poly_func(ref_l)
    r2_ref = r2_score(ref_il, ref_fit)
    ax2.plot(ref_l, ref_il, "b", label="Raw Data", linewidth=1.0)
    ax2.plot(ref_l, ref_fit, "r--", label=f"3rd Order Fit (R^2={r2_ref:.4f})")
    ax2.set_title("Reference Fit")
    ax2.set_xlabel("Wavelength [nm]")
    ax2.set_ylabel("Transmission [dB]")
    ax2.legend()
    ax2.grid(True, ls="--", alpha=0.5)

    ax3 = axes[0, 2]
    for index, sweep in enumerate(sweeps):
        wavelength = sweep["L"]
        il = sweep["IL"]
        assert isinstance(wavelength, np.ndarray)
        assert isinstance(il, np.ndarray)
        processed = il - poly_func(wavelength)
        # The reference sweep has no fringes; cubic subtraction already flattens
        # it to ~0 dB. Running the second-stage (peak-envelope) flatten on it would
        # subtract its noise ceiling and push the whole reference ~1 dB below 0,
        # making the modulator crests appear to rise ABOVE the reference (which is
        # physically impossible). So flatten only the modulator sweeps.
        is_reference = index == len(sweeps) - 1
        if is_reference:
            flattened = processed
        else:
            fsr = measure_fsr(wavelength, processed, fallback=fsr_fallback)
            flattened = flatten_to_envelope(wavelength, processed, fsr)
        ax3.plot(wavelength, flattened,
                 label=sweep_label(sweep, index, len(sweeps)), linewidth=1.0)
    ax3.axhline(0.0, color="black", ls=":", lw=1.0, alpha=0.7)
    ax3.set_title("Transmission spectra - processed (flattened)")
    ax3.set_xlabel("Wavelength [nm]")
    ax3.set_ylabel("Normalized transmission [dB]")
    ax3.legend(ncol=2, fontsize="small", loc="lower center")
    ax3.grid(True, ls="--", alpha=0.5)

    ax4 = axes[1, 0]
    mzm = pick_sweep(sweeps, MOD_BIAS)
    if mzm is not None:
        wavelength = mzm["L"]
        il = mzm["IL"]
        assert isinstance(wavelength, np.ndarray)
        assert isinstance(il, np.ndarray)
        proc_db = il - poly_func(wavelength)
        fsr0 = measure_fsr(wavelength, proc_db, fallback=fsr_fallback)
        # apply the same second-stage flatten before converting to linear scale,
        # so the cos^2 model is not fighting a residual dome (raises LMZO R^2)
        proc_db = flatten_to_envelope(wavelength, proc_db, fsr0)
        transmission = 10 ** (proc_db / 10.0)

        t_low = float(np.nanpercentile(transmission, 5))
        t_high = float(np.nanpercentile(transmission, 95))
        span = max(t_high - t_low, 1e-3)
        p0 = [t_low, span, float(wavelength[np.nanargmax(transmission)]), fsr0, 0.0, 0.0]
        bounds = (
            [max(-0.5, t_low - span), span * 0.05, float(wavelength[0] - fsr0),
             fsr0 * 0.65, -np.pi, -0.05],
            [max(2.0, t_high + span), max(2.0, span * 3.0), float(wavelength[-1] + fsr0),
             fsr0 * 1.45, np.pi, 0.05],
        )

        ax4.plot(wavelength, transmission, color="C0", linewidth=1.1, label=f"Flat Raw ({MOD_BIAS}V)")
        try:
            popt, _ = curve_fit(mzi_model, wavelength, transmission, p0=p0, bounds=bounds, maxfev=20000)
            fit = mzi_model(wavelength, *popt)
            r2 = r2_score(transmission, fit)
            ax4.plot(wavelength, fit, "--", color="black", lw=2,
                     label=f"MZI fit (R^2={r2:.4f}, FSR={popt[3]:.2f}nm)")
            y_values = np.concatenate([transmission[np.isfinite(transmission)], fit[np.isfinite(fit)]])
        except Exception as exc:
            ax4.text(0.5, 0.5, f"Fitting Failed\n{exc}", transform=ax4.transAxes,
                     ha="center", va="center", color="red", fontsize="small")
            y_values = transmission[np.isfinite(transmission)]

        ax4.set_xlim(float(wavelength.min()), float(wavelength.max()))
        if y_values.size:
            y_min = min(-0.05, float(np.nanmin(y_values)) - 0.05)
            y_max = max(1.15, float(np.nanmax(y_values)) * 1.08)
            ax4.set_ylim(y_min, y_max)
    else:
        ax4.text(0.5, 0.5, f"Bias {MOD_BIAS}V\nNot Found", transform=ax4.transAxes,
                 ha="center", va="center", color="red")
    ax4.set_title(f"MZM Fitting ({MOD_BIAS}V) - Linear")
    ax4.set_xlabel("Wavelength [nm]")
    ax4.set_ylabel("Normalized transmission [Linear]")
    ax4.legend(loc="upper right", fontsize="small")
    ax4.grid(True, ls="--", alpha=0.5)

    ax5 = axes[1, 1]
    if iv["V"].size:
        current_abs = np.abs(iv["I"])
        positive = current_abs > 0
        ax5.semilogy(iv["V"][positive], current_abs[positive], "bo", ms=4)
    ax5.set_title("IV-curve (Log Scale)")
    ax5.set_xlabel("Voltage [V]")
    ax5.set_ylabel("Current [A]")
    ax5.grid(True, which="both", ls="--", alpha=0.5)

    ax6 = axes[1, 2]
    if iv["V"].size:
        voltage = iv["V"]
        current_abs = np.abs(iv["I"])
        positive = current_abs > 0
        ax6.semilogy(voltage[positive], current_abs[positive], "o",
                     color="tab:blue", ms=5, label="Measured IV")

        # --- reverse / low-bias region: polynomial fit on log10(|I|) ---
        reverse = (voltage < 0.5) & positive
        r2_rev = float("nan")
        if np.count_nonzero(reverse) >= 4:
            v_rev_pts = voltage[reverse]
            log_i_rev = np.log10(current_abs[reverse])
            # drop sharp single-point notches (diode zero-crossing, where |I| collapses
            # to the noise floor) so the polynomial doesn't oscillate toward that dip
            keep = np.ones(v_rev_pts.size, dtype=bool)
            for i in range(1, v_rev_pts.size - 1):
                if log_i_rev[i] < 0.5 * (log_i_rev[i - 1] + log_i_rev[i + 1]) - 1.0:
                    keep[i] = False
            v_keep, log_keep = v_rev_pts[keep], log_i_rev[keep]
            deg = min(4, v_keep.size - 1)  # cap degree -> smooth, no Runge wiggle
            rev_poly = np.poly1d(np.polyfit(v_keep, log_keep, deg))
            v_rev = np.linspace(float(v_keep.min()), float(v_keep.max()), 200)
            ax6.semilogy(v_rev, 10 ** rev_poly(v_rev), "-", color="tab:orange",
                         label="Reverse polynomial fit")
            r2_rev = r2_score(log_keep, rev_poly(v_keep))

        # --- forward region: diode fit only if the device actually turns on ---
        forward = (voltage >= 0.5) & positive
        diode_is = diode_n = r2_fwd = float("nan")
        turn_on = (np.count_nonzero(forward) >= 2 and
                   (np.log10(current_abs[forward].max()) -
                    np.log10(current_abs[forward].min())) > 1.0)  # >= ~1 decade rise
        if turn_on:
            try:
                popt, _ = curve_fit(
                    diode_eq, voltage[forward], current_abs[forward],
                    p0=[1e-15, 1.5], bounds=([1e-30, 0.5], [1e-3, 10.0]), maxfev=10000,
                )
                diode_is, diode_n = float(popt[0]), float(popt[1])
                v_fwd = np.linspace(float(voltage[forward].min()), float(voltage.max()), 100)
                ax6.semilogy(v_fwd, diode_eq(v_fwd, *popt), "-", color="tab:green",
                             label="Forward diode fit")
                r2_fwd = r2_score(current_abs[forward], diode_eq(voltage[forward], *popt))
            except Exception as exc:
                ax6.text(0.5, 0.2, f"Diode fit failed\n{exc}", transform=ax6.transAxes,
                         ha="center", color="red", fontsize="small")

        # --- stats box (n/a when the diode fit is skipped) ---
        def _fmt(value, spec):
            return "n/a" if not np.isfinite(value) else format(value, spec)

        note = "" if turn_on else "\n(flat IV: no turn-on)"
        stats = (f"Is = {_fmt(diode_is, '.3e')} A\nn = {_fmt(diode_n, '.3f')}\n"
                 f"$R^2_{{fwd}}$ = {_fmt(r2_fwd, '.4f')}\n$R^2_{{rev}}$ = {_fmt(r2_rev, '.4f')}{note}")
        ax6.text(0.03, 0.97, stats, transform=ax6.transAxes, va="top", ha="left",
                 fontsize=9, family="monospace",
                 bbox=dict(boxstyle="square,pad=0.4", fc="white", ec="0.4", lw=0.8))
    ax6.set_title("IV analysis")
    ax6.set_xlabel("Voltage [V]")
    ax6.set_ylabel("Current [A]")
    ax6.grid(True, which="both", ls="--", alpha=0.5)
    ax6.legend(fontsize="small", loc="lower left")

    test_site_info = root.find(".//TestSiteInfo")
    batch = attr_any(test_site_info, "Batch", default="?")
    wafer = attr_any(test_site_info, "Wafer", default="?")
    device = attr_any(test_site_info, "TestSite", default="?")
    die = f"({attr_any(test_site_info, 'DieColumn', default='?')},{attr_any(test_site_info, 'DieRow', default='?')})"
    title = f"Analysis for {wafer} {die} {device}"
    plt.suptitle(title, fontsize=16, y=0.98, fontweight="bold")
    sub = f"Batch: {batch}  |  Wafer: {wafer}  |  Date: {root.attrib.get('CreationDate', '?')}"
    fig.text(0.5, 0.93, sub, ha="center", fontsize=11, color="dimgray")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return True


def analyze_mzm(data_dir: Path = DATA_DIR) -> list[dict[str, object]]:
    xml_files = sorted(data_dir.rglob("*LMZ*.xml"))
    all_rows: list[dict[str, object]] = []
    used_png_names_by_folder: dict[tuple[str, str], set[str]] = defaultdict(set)

    clean_legacy_outputs()
    all_measurement_folders = mirror_result_folders(data_dir)

    total = len(xml_files)
    print(f"Found {total} MZM files", flush=True)
    for index, xml_path in enumerate(xml_files, start=1):
        # in-place progress line: [6/98] updates on the same line
        print(f"\r[{index}/{total}] {xml_path.name:<55}", end="", flush=True)
        try:
            rows = summarize_xml(xml_path)
            all_rows.extend(rows)
            root = ET.parse(xml_path).getroot()
            wafer = attr_any(root.find("./TestSiteInfo"), "Wafer", default="unknown")
            timestamp = xml_path.parent.name
            output_path = unique_png_path(wafer, timestamp, xml_path, used_png_names_by_folder)
            analyze_figure(xml_path, output_path)
        except Exception as exc:
            print(f"\n  ERROR {xml_path}: {exc}", flush=True)
    print(flush=True)  # finish the progress line with a newline

    rows_by_wafer: dict[str, list[dict[str, object]]] = defaultdict(list)
    rows_by_wafer_date: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in all_rows:
        wafer = str(row["wafer"])
        timestamp = str(row["timestamp"])
        rows_by_wafer[wafer].append(row)
        rows_by_wafer_date[(wafer, timestamp)].append(row)

    write_csv(CSV_DIR / "mzm_all_summary.csv", all_rows)
    for wafer, rows in sorted(rows_by_wafer.items()):
        write_csv(CSV_DIR / wafer / f"{wafer}_mzm_summary.csv", rows)
        # plot_wafer_summary(wafer, rows, PNG_DIR / wafer / f"{wafer}_mzm_summary.png")
    for wafer, timestamp in all_measurement_folders:
        rows = rows_by_wafer_date.get((wafer, timestamp), [])
        write_csv(CSV_DIR / wafer / timestamp / f"{wafer}_{timestamp}_mzm_summary.csv", rows)
        # if rows:
        #     plot_wafer_summary(f"{wafer} {timestamp}", rows,
        #                        PNG_DIR / wafer / timestamp / f"{wafer}_{timestamp}_mzm_summary.png")

    return all_rows


def main(data_root: str | Path = DATA_DIR, out_dir: str | Path = PNG_DIR) -> None:
    global PNG_DIR
    PNG_DIR = Path(out_dir)
    rows = analyze_mzm(Path(data_root))
    wafer_count = len({row["wafer"] for row in rows})
    file_count = len({row["source_file"] for row in rows})
    print(f"Analyzed {file_count} MZM XML files across {wafer_count} wafers.")
    print(f"CSV output: {CSV_DIR}")
    print(f"PNG output: {PNG_DIR}")


if __name__ == "__main__":
    main()
