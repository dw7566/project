from __future__ import annotations

import statistics
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from . import config
from .config import DATA_DIR, CSV_DIR, MOD_BIAS
from .xml_parser import (
    attr_any, r2_score, load_xml, pick_sweep, sweep_label,
)
from .spectrum import (
    mzi_model, measure_fsr, flatten_to_envelope, device_fsr_fallback,
)
from .iv_analysis import plot_iv_log, plot_iv_analysis
from .vpi_analysis import plot_vpi_voltage_panels, analyze_vpi_voltage_figure
from .extinction_ratio import (
    plot_extinction_ratio_panels, analyze_extinction_ratio_figure,
)
from .csv_export import summarize_xml, write_csv


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
    if config.PNG_DIR.exists():
        for wafer_dir in config.PNG_DIR.iterdir():
            if not wafer_dir.is_dir():
                continue
            for png_path in wafer_dir.glob("*LMZ*.png"):
                png_path.unlink()

    if CSV_DIR.exists():
        for csv_path in CSV_DIR.glob("D*_mzm_summary.csv"):
            csv_path.unlink()

    old_summary_png = config.PNG_DIR / "vpi_summary.png"
    if old_summary_png.exists():
        old_summary_png.unlink()

    legacy_modulation_dir = config.PNG_DIR / config.MODULATION_EFFICIENCY_PNG_DIR
    if legacy_modulation_dir.exists():
        for png_path in legacy_modulation_dir.rglob("*LMZ*.png"):
            png_path.unlink()


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
        (config.PNG_DIR / wafer / timestamp).mkdir(parents=True, exist_ok=True)
        (CSV_DIR / wafer / timestamp).mkdir(parents=True, exist_ok=True)
    return folders


def unique_png_path(wafer: str, timestamp: str, xml_path: Path,
                    used_names_by_folder: dict[tuple[str, str], set[str]]) -> Path:
    output_dir = config.PNG_DIR / wafer / timestamp
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

    fig, axes = plt.subplots(4, 3, figsize=(20, 20))
    plt.subplots_adjust(hspace=0.45, wspace=0.3)

    # Row 0, Col 0
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

    # Row 0, Col 1
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

    # Row 0, Col 2
    ax3 = axes[0, 2]
    for index, sweep in enumerate(sweeps):
        wavelength = sweep["L"]
        il = sweep["IL"]
        assert isinstance(wavelength, np.ndarray)
        assert isinstance(il, np.ndarray)
        processed = il - poly_func(wavelength)
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

    # Row 1, Col 0
    ax4 = axes[1, 0]
    mzm = pick_sweep(sweeps, MOD_BIAS)
    if mzm is not None:
        wavelength = mzm["L"]
        il = mzm["IL"]
        assert isinstance(wavelength, np.ndarray)
        assert isinstance(il, np.ndarray)
        proc_db = il - poly_func(wavelength)
        fsr0 = measure_fsr(wavelength, proc_db, fallback=fsr_fallback)
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

    # Row 1, Col 1 & 2
    plot_iv_log(axes[1, 1], iv)
    plot_iv_analysis(axes[1, 2], iv)

    # Row 2: V_pi vs voltage
    plot_vpi_voltage_panels(axes[2, :], root)

    # Row 3: 소광비
    plot_extinction_ratio_panels(axes[3, :], root)

    # 타이틀
    test_site_info = root.find(".//TestSiteInfo")
    batch = attr_any(test_site_info, "Batch", default="?")
    wafer = attr_any(test_site_info, "Wafer", default="?")
    device = attr_any(test_site_info, "TestSite", default="?")
    die = f"({attr_any(test_site_info, 'DieColumn', default='?')},{attr_any(test_site_info, 'DieRow', default='?')})"
    title = f"Analysis for {wafer} {die} {device}"
    plt.suptitle(title, fontsize=16, y=0.98, fontweight="bold")
    sub = f"Batch: {batch}  |  Wafer: {wafer}  |  Date: {root.attrib.get('CreationDate', '?')}"
    fig.text(0.5, 0.95, sub, ha="center", fontsize=11, color="dimgray")

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
        print(f"\r[{index}/{total}] {xml_path.name:<55}", end="", flush=True)
        try:
            rows = summarize_xml(xml_path)
            all_rows.extend(rows)
            root = ET.parse(xml_path).getroot()
            wafer = attr_any(root.find("./TestSiteInfo"), "Wafer", default="unknown")
            timestamp = xml_path.parent.name
            output_path = unique_png_path(wafer, timestamp, xml_path, used_png_names_by_folder)
            analyze_figure(xml_path, output_path)
            analyze_vpi_voltage_figure(xml_path, root)
            analyze_extinction_ratio_figure(xml_path, root)
        except Exception as exc:
            print(f"\n  ERROR {xml_path}: {exc}", flush=True)
    print(flush=True)

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
    for wafer, timestamp in all_measurement_folders:
        rows = rows_by_wafer_date.get((wafer, timestamp), [])
        write_csv(CSV_DIR / wafer / timestamp / f"{wafer}_{timestamp}_mzm_summary.csv", rows)

    return all_rows


def main(data_root: str | Path = DATA_DIR, out_dir: str | Path = None) -> None:
    if out_dir is not None:
        config.PNG_DIR = Path(out_dir)
    rows = analyze_mzm(Path(data_root))
    wafer_count = len({row["wafer"] for row in rows})
    file_count = len({row["source_file"] for row in rows})
    print(f"Analyzed {file_count} MZM XML files across {wafer_count} wafers.")
    print(f"CSV output: {CSV_DIR}")
    print(f"PNG output: {config.PNG_DIR}")
