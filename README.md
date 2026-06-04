![header](https://capsule-render.vercel.app/api?type=wave&color=auto&height=200&section=header&text=SPDAP&fontSize=70)

###### Silicon Photonics Data Analysis Automation Pipeline
### Contents

[1. Introduction](#1-introduction)\
[2. Project information](#2-project-information)\
[3. Install and Run](#3-install-and-run)\
[4. Description of the module file feature](#4-description-of-the-module-file-feature)\
[5. Run file algorithm](#5-run-file-algorithm)\
[6. Output Examples](#6-output-examples)\
[7. CSV Column Description](#7-csv-column-description)\
[8. Configuration](#8-configuration)\
[9. Project Structure](#9-project-structure)

---

# SPDAP :
##### Hi !
##### Thank you for looking at our project. SPDAP stands for "Silicon Photonics Data Analysis Automation Pipeline."
##### This project automates silicon photonics wafer-scale MZM measurement analysis and generates organized CSV and PNG results.

---

## 1. Introduction

We aim to develop a Python-based automation pipeline for silicon photonics data analysis.
Wafer-scale MZM (Mach-Zehnder Modulator) measurement data analysis has been tedious and error-prone when done manually.
This project solves that by automating the entire process.

The goal is to process wafer-scale XML measurement data, extract MZM-related device information, key optical/electrical characteristics, and save comprehensive analysis outputs.

#### - Main Features
- **Parsing**: Extract wavelength sweep, insertion loss, voltage, and current data from MZM XML files
- **Fitting**: Perform optical spectrum normalization, MZM parameter fitting via FSR detection, and R² quality metrics
- **Visualization**: Generate die-level 12-panel analysis figures and wafer-level extinction ratio + bias analysis heatmaps
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
    + Reference sweep, Dark/Light curves, and Extinction Ratio

####
+ **Run file description**

   SPDAP scans the `data` directory and extracts XML files whose names include `LMZ`.\
   Then, it parses wavelength sweep, insertion loss, voltage, and current data, runs fitting/analysis logic, and saves CSV tables and PNG figures under the `res` directory.\
   The pipeline includes reference fitting, MZM FSR detection, V_pi analysis, and extinction ratio extraction across multiple DC bias points.

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
  + Plots dark current IV curves in log scale with fitting diagnostics
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

* **Modulation Efficiency module** (`src/modulation_efficiency.py`)
  + Fits transmission spectra in dB scale using residual MZI model
  + Performs envelope detection on flattened spectra
  + Extracts modulation efficiency metrics and FSR from residuals

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

---

## 5. Run file algorithm

* **Preparation**
   + Scan `data/` directory for XML files matching `*LMZ*.xml` pattern
   + Create result directory hierarchy mirroring data structure

* **Execution (Per XML File)**
   1. **Parse XML**: Load metadata, sweeps (by DC bias), and IV data
   2. **Reference Fit**: Fit reference sweep with 3rd-order polynomial
   3. **Spectrum Processing**: 
      - Flatten each biased sweep using reference and FSR detection
      - Extract transmission in both linear and dB scales
   4. **Parameter Extraction**:
      - MZM FSR fitting in dB scale with residual envelope
      - Extinction ratio calculation across bias points
      - Vπ analysis from phase response
      - IV curve analysis
   5. **Generate Figure**: Create 12-panel PNG with all analysis results
   6. **Extract CSV Row**: Summarize key parameters (IL stats, ER, Vπ, R², source file)

* **Aggregation**
   + Collect all CSV rows from individual files
   + Write wafer-level CSV: `res/csv/{wafer_id}/{timestamp}.csv`
   + Write global CSV: `res/csv/mzm_all_summary.csv`
   + Generate wafer summary figures: extinction ratio by die + by bias heatmaps

* **Output**
   + Analysis figures: `res/png/{wafer_id}/{timestamp}/*.png`
   + Wafer map: `res/png/{wafer_id}/{timestamp}/wafermap.png`
   + CSV data: `res/csv/mzm_all_summary.csv` + `res/csv/{wafer_id}/{timestamp}.csv`

---

## 6. Output Examples

### Die-Level Analysis Figure (12-panel)

```
┌──────────────────────┬──────────────────────┬──────────────────────┐
│ (0,0) Transmission   │ (0,1) Reference Fit  │ (0,2) Flattened      │
│       as measured    │       (3rd order)    │       spectra        │
├──────────────────────┼──────────────────────┼──────────────────────┤
│ (1,0) MZM Linear     │ (1,1) IV Log scale   │ (1,2) IV Analysis    │
│       FSR fit        │       fitting        │                      │
├──────────────────────┼──────────────────────┼──────────────────────┤
│ (2,0) MZM dB         │ (2,1) Vπ Voltage     │ (2,2) Extinction     │
│       Residual fit   │       curves         │       Ratio          │     │
├──────────────────────┼──────────────────────┼──────────────────────┤
```

### Wafer-Level Summary Figure

Two panels:
- **Left**: Die position (col, row) heatmap of best extinction ratio by die
- **Right**: Extinction ratio vs DC bias (mean ± range)

---

## 7. CSV Column Description

| Column Name | Unit | Description |
|-------------|------|-------------|
| `lot` | — | Lot ID from XML |
| `wafer` | — | Wafer ID (e.g., D08, D24) |
| `test_site` | — | Test site identifier |
| `die_column` | — | Die column index |
| `die_row` | — | Die row index |
| `timestamp` | — | Measurement timestamp folder name |
| `device_name` | — | Device descriptor from XML |
| `dc_bias_v` | V | DC bias voltage applied to MZM |
| `current_at_minus_1v_a` | A | Current measured at -1V |
| `current_at_0v_a` | A | Current measured at 0V |
| `current_at_plus_1v_a` | A | Current measured at +1V |
| `il_min_db` | dB | Minimum insertion loss in sweep |
| `il_max_db` | dB | Maximum insertion loss in sweep |
| `il_mean_db` | dB | Mean insertion loss in sweep |
| `extinction_ratio_db` | dB | Extinction ratio (IL_max - IL_min) |
| `wavelength_at_min_il_nm` | nm | Wavelength where IL is minimum |
| `wavelength_at_max_il_nm` | nm | Wavelength where IL is maximum |
| `modulation_r2_by_null` | — | R² of MZI fit on dB-scale residuals |
| `vpi_mean_v` | V | Mean Vπ across measurement conditions |
| `vpi_min_v` | V | Minimum Vπ measured |
| `vpi_max_v` | V | Maximum Vπ measured |
| `vpi_by_null_v` | V | Vπ extracted from fringe null method |
| `source_file` | — | Original XML filename |

---

## 8. Configuration

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

```
## 9. Project Structure

project/
├── run.py                      # Execution entry point
├── src/
│   ├── __init__.py
│   ├── config.py               # Configuration: paths, bias, constants
│   ├── main.py                 # Main pipeline orchestrator
│   ├── xml_parser.py           # XML parsing & data extraction
│   ├── spectrum.py             # Spectral processing & MZI model
│   ├── iv_analysis.py          # IV curve plotting & analysis
│   ├── extinction_ratio.py     # Extinction ratio calculation
│   ├── vpi_analysis.py         # Vπ extraction & analysis
│   ├── modulation_efficiency.py # Modulation fitting in dB scale
│   ├── csv_export.py           # CSV row summarization & export
│   ├── wafermap.py             # Wafer-level heatmap generation
│   ├── datalocation.py         # Data path utilities
│   ├── runfile.py              # Run configuration helpers
│   ├── tocsv.py                # CSV writing utilities
│   └── example.py              # Example usage
├── data/                       # Input directory for XML measurement files
│   └── (organized by wafer_id/timestamp/)
├── res/                        # Output directory for results
│   ├── csv/                    # CSV analysis results
│   │   ├── mzm_all_summary.csv # Global summary
│   │   └── {wafer_id}/         # Per-wafer CSVs
│   │       └── {timestamp}.csv
│   └── png/                    # Generated analysis figures
│       └── {wafer_id}/
│           └── {timestamp}/
│               ├── *.png       # Die-level analysis figures
│               └── wafermap.png # Wafer summary heatmap
├── requirements.txt            # Python dependencies
└── README.md                   # This file


## Requirements

This project requires Python 3.7+ with the following dependencies:

- `numpy` — Numerical computing
- `scipy` — Scientific computing and curve fitting
- `matplotlib` — Data visualization
- `lxml` — XML parsing

See `requirements.txt` for pinned versions.

---

## License

This project is developed by the Silicon Photonics Research Team at Hanyang University.

---

**Last Updated**: 2026-06-04
