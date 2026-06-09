![header](https://capsule-render.vercel.app/api?type=wave&color=auto&height=200&section=header&text=SPDAP&fontSize=70)

###### Silicon Photonics Data Analysis Automation Pipeline
### Contents

1. [Introduction](#1-introduction)
2. [Project Information](#2-project-information)
3. [Installation and Usage](#3-installation-and-usage)
4. [Module Descriptions](#4-module-descriptions)
5. [Output Examples](#5-output-examples)
6. [CSV Column Definitions](#6-csv-column-definitions)
7. [Configuration](#7-configuration)
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
| Kim dong min  |  dm1656@hanyang.ac.kr   |
| Kim sang wook | tkdlek850@hanyang.ac.kr |
| Lee jae hyeok |  dw7566@hanyang.ac.kr   |

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

## 3. Install and Run

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
    - Wafer-level summary figures are saved in `res/png/{wafer_id}/{timestamp}/wafermap.png`

---

## 4. Description of the module file feature

* **XML Parser module** (`src/xml_parser.py`)
  + Extracts raw wavelength, insertion loss, DC bias, and measurement metadata from XML files
  + Organizes sweep data by bias voltage and identifies reference sweep
  + Provides utility functions for robust attribute extraction and R² calculation

* **Spectrum Processing module** (`src/spectrum.py`)
  + Implements MZI model for fitting transmission spectra
  + Detects Free Spectral Range (FSR) from processed data
  + Flattens spectral data to envelope for robust peak extraction
  + Normalizes transmission using reference sweep polynomial fit

* **IV Analysis module** (`src/iv_analysis.py`)
  + Extracts voltage and current relationships
  + Supports multiple bias conditions

* **Extinction Ratio module** (`src/extinction_ratio.py`)
  + Calculates extinction ratio from high and low transmission points
  + Generates extinction ratio vs DC bias analysis plots
  + Validates and filters data based on measurement quality

* **Vπ Analysis module** (`src/vpi_analysis.py`)
  + Extracts half-wave voltage (Vπ) from phase/modulation characteristics
  + Analyzes Vπ variation across different wavelengths and measurement conditions
  + Provides statistical summary (mean, min, max)

* **CSV Export module** (`src/csv_export.py`)
  + Summarizes key parameters from each XML file into a structured row
  + Consolidates multi-file results into wafer-level CSV tables
  + Manages timestamp-level organization of CSV outputs

* **Wafer Map module** (`src/wafermap.py`)
  + Generates 2D heatmaps of die performance metrics (extinction ratio, bias response)
  + Maps die position (column, row) to measured parameters
  + Provides visual overview of wafer-scale uniformity

* **Main Pipeline module** (`src/main.py`)
  + Orchestrates the entire analysis workflow
  + Iterates over all LMZ XML files in the data directory
  + Calls all analysis and visualization routines
  + Manages result organization and file naming

* **CSV Utilities module** (`src/tocsv.py`)
  + Provides helper functions for CSV file operations
  + Manages column ordering and data validation

---

## 5. Output Examples

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

## 6. CSV Column Description

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

## 7. Configuration

### Modify Data Directory

Edit `src/config.py` to change input/output paths:

```python
# src/config.py
DATA_DIR = Path("data")
CSV_DIR = Path("res") / "csv"
PNG_DIR = Path("res") / "png"
```

### Modify Modulation Bias Point

Change the bias voltage used for MZM FSR fitting:

```python
# src/config.py
MOD_BIAS = "-1.0"  # Modulation bias in Volts
```

### Thermal Voltage Constant

Adjust thermal voltage (affects IV curve fitting):

```python
# src/config.py
THERMAL_VOLTAGE = 0.02585  # At room temperature ~25°C
```

### Output Figure DPI

Modify resolution in `src/main.py`:

```python
fig.savefig(out_path, dpi=110, bbox_inches="tight")
```

### Figure Size and Layout

Adjust figure dimensions and spacing in `src/main.py`:

```python
# Current: 3x3 grid, 18x15 inches
fig, axes = plt.subplots(3, 3, figsize=(18, 15))
fig.subplots_adjust(hspace=0.45, wspace=0.3)
```

---

## 8. Project Structure

```
project/
├── run.py                      # Execution entry point
├── requirements.txt            # Python dependencies
├── README.md                   # This file
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
│   └── tocsv.py                # CSV utilities
├── doc/
│   └── project_explanation.ipynb # Jupyter notebook with detailed explanations
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

For a detailed walkthrough and explanation of the project, see `doc/project_explanation.ipynb` (Jupyter Notebook).

---

## License

This project is developed by the Silicon Photonics Research Team at Hanyang University.

---

**Last Updated**: 2026-06-06
