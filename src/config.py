from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_DIR = Path("data")
CSV_DIR = Path("res") / "csv"
PNG_DIR = Path("res") / "png"

THERMAL_VOLTAGE = 0.02585
MOD_BIAS = "-1.0"

plt.rcParams["font.family"] = ["DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

CSV_COLUMNS = [
    "lot", "wafer", "test_site", "die_column", "die_row", "timestamp",
    "device_name", "dc_bias_v", "current_a",
    "extinction_ratio_db",
    "vpi_at_dc_bias_v",
    "source_file",
]