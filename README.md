![header](https://capsule-render.vercel.app/api?type=wave&color=auto&height=200&section=header&text=SPDAP&fontSize=70)

###### Silicon Photonics Data Analysis Automation Pipeline
### Contents

[1. Introduction](#1-introduction)\
[2. Project information](#2-project-information)\
[3. Install and Run](#3-install-and-run)\
[4. Description of the module file feature](#4-description-of-the-module-file-feature)\
[5. Run file algorithm](#5-run-file-algorithm)

---

# SPDAP :
##### Hi !
##### Thank you for looking at our project. SPDAP stands for "Silicon Photonics Data Analysis Automation Pipeline."
##### This project automates silicon photonics wafer-scale measurement analysis and generates organized CSV and PNG results.

---

## 1. Introduction
We aim to develop a Python-based automation pipeline for silicon photonics data analysis.
The goal is to process wafer-scale XML measurement data, extract MZM-related device information, key optical/electrical characteristics, and save analysis outputs.

#### - Main Features
- **Parsing**: Extract wavelength sweep, insertion loss, voltage, and current data from MZM XML files
- **Fitting**: Perform optical spectrum normalization, MZM parameter fitting, and quality estimation (R-squared)
- **Visualization**: Generate die-level 6-panel analysis figures and wafer-level heatmaps
- **CSV Output**: Organize and save measurement data and fitting parameters

#### - contributors : If you have any questions, please contact us at the following email.

|     name      |         E-mail          |
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
    + Optical spectrum and IV characteristics

####
+ **Run file description**

   SPDAP scans the `data` directory and extracts XML files whose names include `LMZ`.\
   Then, it parses wavelength sweep, insertion loss, voltage, and current data, runs fitting/summary logic, and saves CSV tables and PNG figures under the `res` directory.

---

## 3. Install and Run

####
* Getting Started
   + Enter the Terminal and install required packages. \
```bash
pip install -r requirements.txt
```

* How to Run

  + **Data Preparation**: Place the raw XML measurement data under the `data` directory.

  + **Run the main script**. SPDAP will analyze available MZM XML files and save results automatically.

```bash
python run.py
```

  + **Results**: CSV outputs are saved in `res/csv`, and generated analysis figures are saved in `res/png`.

---

## 4. Description of the module file feature

* **Fitting module**
   + The graph is drawn by parsing raw wavelength sweep, insertion loss, current, and voltage data from XML files.
   + The module normalizes transmission spectra using a reference sweep, performs MZM fitting, estimates fitting quality such as R-squared, and visualizes IV behavior and optical response.

* **CSV module**
  + It contains a variety of measurement information, including lot, wafer, test site, die column, die row, timestamp, device name, bias, current values, wavelength range, and insertion-loss statistics.
  + Creates a dataframe-style CSV summary so that the measured information in the XML files can be viewed at a glance.
  + Saves this data frame in CSV format in the `res/csv` folder.

---

## 5. Run file algorithm

* **Preparation**
   + Place the measurement XML files under the `data` directory.

* **Execution**
   + The `main` function in `run.py` is executed.
   + `src/MZMfitting.py` searches for MZM XML files, creates analysis figures, and writes wafer-level and timestamp-level CSV summaries.

* **Output**
   + Analysis results are generated automatically and saved in the `res` directory.
   + CSV summaries are stored in `res/csv` folder.
   + Generated figures are stored in `res/png` folder.

---

## Project Structure

```
project/
├── data/                    # Input directory for XML measurement files
├── res/                     # Output directory for results
│   ├── csv/                # CSV analysis results
│   └── png/                # Generated analysis figures
├── src/
│   ├── MZMfitting.py       # Main fitting and analysis module
│   └── ...                 # Other module files
├── run.py                  # Main execution script
├── requirements.txt        # Python dependencies
└── README.md              # This file
```

---

## Requirements

This project requires Python 3.7+ with the following dependencies:
- See `requirements.txt` for a complete list of required packages.

---

## License

This project is developed by the Silicon Photonics Research Team at Hanyang University.

---

**Last Updated**: 2026-06-04
