![header](https://capsule-render.vercel.app/api?type=wave&color=auto&height=200&section=header&text=SPDAP&fontSize=70)

###### Silicon Photonics Data Analysis Automation Pipeline
### Contents

1. [Introduction](#1-introduction)
2. [Project Information](#2-project-information)
3. [Project Structure](#3project-structure)
4. [Install and Run](#4-install-and-run)
5. [Module Descriptions](#5-description-of-the-module-file-feature)
6. [Output Examples](#6-output-examples)
7. [CSV Column Definitions](#7-csv-column-description)
8. [Project Structure](#8-project-structure)

---

# SPDAP :
##### Hi !
##### Thank you for looking at our project. SPDAP stands for "Silicon Photonics Data Analysis Automation Pipeline."
##### This project automates silicon photonics wafer-scale MZM measurement analysis and generates organized CSV and PNG results.

---

## 1. Introduction

We aim to develop a Python-based automation pipeline for silicon photonics data analysis.
Wafer-scale MZM (Mach-Zehnder Modulator) measurement data analysis has been error-prone when done manually.

This project solves that by automating the entire process.

The goal is to process wafer-scale XML measurement data, extract MZM-related device information, key optical/electrical characteristics, and save comprehensive analysis outputs.

#### - Main Features
- **Parsing**: Extract wavelength sweep, insertion loss, voltage, and current data from MZM XML files
- **Fitting**: Perform optical spectrum normalization, MZM parameter fitting via FSR detection, and R² quality metrics
- **Visualization**: Generate die-level 9-panel analysis figures and wafer-level extinction ratio + bias analysis heatmaps
- **CSV Output**: Organize and save measurement data and fitting parameters with wafer/device/timestamp hierarchy

#### - Contributors

If you have any questions, please contact us at the following email.

|     Name      |         E-mail          |
|:-------------:|:-----------------------:|
| Dong-Min Kim  |  dm1656@hanyang.ac.kr   |
|  Sang-Uk Kim  | tkdlek850@hanyang.ac.kr |
| Jae-Hyeok Lee |  dw7566@hanyang.ac.kr   |

---

## 2. Project information

####
+ **Detailed project**

    Main analysis targets
    + Wafer
    + Die row & column
    + MZM XML measurement files
    + Optical spectrum (wavelength sweep) and IV characteristics
    + Reference sweep, Extinction Ratio

####

---


## 3. Project Structure

```
project/
├── run.py                      # Execution entry point
├── requirements.txt            # Python dependencies
├── README.md                   # This file
├── project_explanation.ipynb   # Jupyter notebook with detailed explanations
├── src/
│   ├── config.py               # Configuration: paths, bias, constants
│   ├── main.py                 # Main pipeline orchestrator
│   ├── xml_parser.py           # XML parsing & data extraction
│   ├── spectrum.py             # Spectral processing & MZI model
│   ├── iv_analysis.py          # IV curve plotting & analysis
│   ├── extinction_ratio.py     # Extinction ratio calculation
│   ├── vpi_analysis.py         # Vπ extraction & analysis
│   ├── csv_export.py           # CSV row summarization & export
│   ├── wafermap.py             # Wafer-level heatmap generation
│
├── data/                       # Input directory for XML measurement files
│   └── (organized by wafer_id/timestamp/)
└── res/                        # Output directory for results
    ├── csv/                    # CSV analysis results
    │   ├── mzm_all_summary.csv # Global summary
    │   └── {wafer_id}/         # Per-wafer CSVs
    │       └── {timestamp}.csv
    └── png/                    # Generated analysis figures
        └── {wafer_id}/
            └── {timestamp}/
                ├── *.png       # Die-level analysis figures (3x3 panels)
                └── wafermap.png # Wafer summary heatmap
```


## 4. Install and Run

####
* **Getting Started**
   + Enter the Terminal and install required packages from `requirements.txt`:
   
```bash
pip install -r requirements.txt
```

* **How to Run**

  + **Data Preparation**: Place the raw XML measurement data under the `data` directory.

```
data/
├── {wafer_id}/
│   ├── {timestamp}/
│   │   ├── *LMZ*.xml
│   │   └── ...
│   └── ...
└── ...
```

  + **Run the main script**. SPDAP will analyze available MZM XML files and save results automatically.

```bash
python run.py
```

  + **Results**: 
    - CSV outputs are saved in `res/csv/`
    - PNG die-level figures are saved in `res/png/{wafer_id}/{timestamp}/`
    - Wafer-map figures are saved in `res/png/{wafer_id}/{timestamp}/wafermap.png`


---

## 5. Description of the module file feature


* **xml_parser.py** — `XML parsing & shared utilities`
  + Parses raw XML measurement files and provides shared utility functions used across all modules.
  + Extracts wavelength sweep data, IV measurements, and device metadata from the nested XML structure (`PortCombo > IVMeasurement`, `PortCombo > WavelengthSweep`).
  + Key functions: `load_xml()` for sweep/IV extraction, `find_mzm_modulators()` for device filtering, `interpolate_sweeps()` for bias-aligned wavelength grids, and `r2_score()` for fit quality.
```python
# Locate the IVMeasurement nested under PortCombo
port_combo = modulator.find("./PortCombo")
voltage = parse_float_list(port_combo.findtext("./IVMeasurement/Voltage"))
current = parse_float_list(port_combo.findtext("./IVMeasurement/Current"))

# Align multi-bias sweeps onto a common wavelength grid
biases, ref_wl, il_matrix = interpolate_sweeps(parse_modulation_sweeps(modulator))
```

* **iv_analysis.py** — `IV curve plotting & diode fitting`
  + Plots the measured IV curve on a log-scale axis and fits it with a diode model using `scipy.optimize.curve_fit`.
  + The reverse-bias region is fitted with a polynomial, and the forward-bias region uses the Shockley diode equation to extract saturation current `Is` and ideality factor `n`.
```python
def diode_eq(voltage, Is, n):
    exponent = np.clip(voltage / (n * THERMAL_VOLTAGE), -700, 700)
    return Is * (np.exp(exponent) - 1.0)

popt, _ = curve_fit(
    diode_eq, voltage[forward], current_abs[forward],
    p0=[1e-15, 1.5], bounds=([1e-30, 0.5], [1e-3, 10.0])
)
# Result: Is, n, R²_fwd, R²_rev annotated on the plot
```

* **spectrum.py** — `Spectral flattening & MZI model fitting`
  + Normalizes transmission spectra by removing the insertion-loss envelope, then fits a cosine-squared MZI model.
  + The envelope is estimated by finding crest points across FSR windows and interpolating with `PchipInterpolator`. FSR is detected automatically via `scipy.signal.find_peaks`.
```python
def mzi_model(wavelength, A, B, wl0, FSR, phi, slope):
    x = wavelength - wavelength.mean()
    return A + slope * x + B * np.cos(np.pi * (wavelength - wl0) / FSR + phi) ** 2

# Flatten raw IL by subtracting the top envelope before fitting
flattened = flatten_to_envelope(wavelength, il_db, fsr=device_fsr_fallback(root))
```

* **extinction_ratio.py** — `Extinction ratio analysis`
  + Extracts peak-to-null extinction ratio (ER) for each DC bias by detecting local maxima and minima in the interpolated IL spectra. Pairs with ER < 10 dB are filtered as false nulls.
  + Plots ER vs. DC bias (mean/min/max range) and a per-fringe bar chart at 0 V.
```python
for mi in minima_idx:
    peak_il = max(il[left_peak], il[right_peak])   # nearest adjacent peak
    er = peak_il - il[mi]                           # ER in dB
    if er < MIN_ER_DB:                              # discard shallow fringes
        continue
    er_list.append(er)

```

* **vpi_analysis.py** — `Vπ & dλ/dV extraction`
  + Tracks deep-null wavelength positions (< −30 dB) across bias sweeps and fits a linear `dλ/dV` slope per null. Vπ is computed from FSR and the fitted slope.
  + Filters unreliable tracks with R² < 0.5 or |dλ/dV| < 0.02 nm/V to ensure only physically meaningful results are exported.
```python
# Track how each deep null shifts with voltage
coeffs = np.polyfit(v_arr, wl_arr, 1)   # linear fit: wl = coeffs[0]*V + coeffs[1]
dlambda_dv = coeffs[0]                  # nm/V

# Vπ from FSR and modulation efficiency
Vpi = fsr / (2.0 * abs(dlambda_dv))
```

* **csv_export.py** — `CSV summary generation`
  + Iterates over all MZM modulators in each XML file and assembles one row per wavelength sweep, combining device metadata, IV spot values, IL statistics, extinction ratio, and modulation efficiency.
  + Saves the aggregated table to `res/csv/` using Python's built-in `csv.DictWriter`.
```python
rows.append({
    "lot": lot, "wafer": wafer, "die_column": die_column, "die_row": die_row,
    "dc_bias_v": sweep.get("DCBias"),
    "current_at_minus_1v_a": nearest_value(voltage, current, -1.0),
    "il_min_db": il_min, "il_max_db": il_max,
    "extinction_ratio_db": il_max - il_min,
    **modulation,   # dλ/dV, Vπ, FSR from modulation_efficiency
})
write_csv(path, rows)
```

* **src/wafermap.py** - `Wafer Map module`
  + Generates 2D heatmaps of die performance metrics (extinction ratio, bias response)
  + Maps die position (column, row) to measured parameters
  + Provides visual overview of wafer-scale uniformity

* **src/main.py** - `Main Pipeline module`
  + Orchestrates the entire analysis workflow
  + Iterates over all LMZ XML files in the data directory
  + Calls all analysis and visualization routines
  + Manages result organization and file naming

---

## 6. Output Examples

### Die-Level Analysis Figure (9-panel, 3x3)

<img width="937" height="873" alt="image" src="https://github.com/user-attachments/assets/713eb442-ca66-46d6-888f-f134825e9806" />

```
```

| Panel | Content |
|-------|---------|
| (0,0) | Raw transmission spectra for all bias conditions |
| (0,1) | Reference sweep with 3rd-order polynomial fit + R² |
| (0,2) | Flattened (processed) spectra after reference normalization |
| (1,0) | MZI model fit on linear-scale transmission (FSR, R² shown) |
| (1,1) | Dark current IV in log-log scale |
| (1,2) | Forward bias Shockley + reverse bias analysis |
| (2,0) | MZI model fit on dB-scale residuals (FSR detection) |
| (2,1) | Vπ vs wavelength + statistics |
| (2,2) | Extinction ratio vs DC bias across measurement points |

<img width="1586" height="743" alt="image" src="https://github.com/user-attachments/assets/cf059dff-7e9f-4160-9966-a3551f1a0c39" />

---

## 7. CSV Column Description

| Column Name           | Unit | Description                         |
|-----------------------|------|-------------------------------------|
| `lot`                 | — | Lot ID from XML                     |
| `wafer`               | — | Wafer ID (e.g., D08, D24)           |
| `test_site`           | — | Test site identifier                |
| `die_column`          | — | Die column index                    |
| `die_row`             | — | Die row index                       |
| `timestamp`           | — | Measurement timestamp folder name   |
| `device_name`         | — | Device descriptor from XML          |
| `dc_bias_v`           | V | DC bias voltage applied to MZM      |
| `current_a`           | A | Current measured at DC bias         |
| `extinction_ratio_db` | dB | Extinction ratio (IL_max - IL_min)  |
| `vpi_at_dc_bias_v`    | V | Vπ at DC bias                       |
| `source_file`         | — | Original XML filename               |


---

## Requirements

This project requires Python 3.7+ with the following dependencies:

- `numpy` — Numerical computing
- `scipy` — Scientific computing and curve fitting
- `matplotlib` — Data visualization
- `lxml` — XML parsing

See `requirements.txt` for pinned versions.

---

## Documentation

For a detailed walkthrough and explanation of the project, see `project_demo.ipynb`
(Jupyter Notebook).

---

## License

This project is developed by silicon photonics team at Hanyang University

---

**Last Updated**: 2026-06-09
