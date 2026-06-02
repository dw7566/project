from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_DIR = Path("data")
CSV_DIR = Path("res") / "csv"
PNG_DIR = Path("res") / "png"
MODULATION_EFFICIENCY_PNG_DIR = "07_wavelength_modulation_efficiency"
VPI_VOLTAGE_PNG_DIR = "07_vpi_vs_voltage"
EXTINCTION_RATIO_PNG_DIR = "08_extinction_ratio"

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
    "vpi_mean_v", "vpi_min_v", "vpi_max_v", "vpi_by_null_v",
    "source_file",
]
