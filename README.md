![header](https://capsule-render.vercel.app/api?type=wave&color=auto&height=200&section=header&text=SPDAP&fontSize=70)

###### Silicon Photonics Data Analysis Automation Pipeline
### Contents

[1. Introduction](#1-introduction)\
[2. Project information](#2-project-information)\
[3. Run](#3-install-and-run)\
[4. Description of the module file feature](#4-description-of-the-module-file-feature)

---

# SPDAP :
##### Hi !
##### Thank you for looking at our project. SPDAP stands for "Silicon Photonics Data Analysis Automation Pipeline."
##### This project automates silicon photonics wafer-scale measurement analysis and generates organized CSV and PNG results.

---

## 1. Introduction
We aim to develop a Python-based automation pipeline for silicon photonics data analysis.
The goal is to process wafer-scale XML measurement data, extract MZM-related device information, key optical/electrical characteristics, and save analysis outputs.

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

   First, SPDAP scans the `data` directory and extracts XML files whose names include `LMZ`.\
   Then, it parses wavelength sweep, insertion loss, voltage, and current data, runs fitting/summary logic, and saves CSV tables and PNG figures under the `res` directory.

---

## 3. Install and Run

####
* Getting Started
   + Enter the Terminal, write down `pip install -r requirements.txt` and download the required packages. \
```powershell
pip install -r requirements.txt
```

* How to Run
  + Place the raw XML measurement data under the `data` directory.

  + Run the main script. SPDAP will analyze available MZM XML files and save results automatically.

```powershell
python run.py
```

  + CSV outputs are saved in `res/csv`, and generated analysis figures are saved in `res/png`.

---

## 4. Description of the module file feature

* Fitting module
   + The graph is drawn by parsing raw wavelength sweep, insertion loss, current, and voltage data from XML files.
   + The module normalizes transmission spectra using a reference sweep, performs MZM fitting, estimates fitting quality such as R-squared, and visualizes IV behavior and optical response.

* CSV module
  + It contains a variety of measurement information, including lot, wafer, test site, die column, die row, timestamp, device name, bias, current values, wavelength range, insertion-loss statistics, and extinction ratio.
  + Create a dataframe-style CSV summary so that the measured information in the XML files can be viewed at a glance.
  + Save this data frame in CSV format in the `res/csv` folder.

 ---
## 5. Run file algorithm
* First of all
   + Prepare the measurement XML files under the `data` directory.


* And then
  + The `main` function in `run.py` is executed.
* Next
   + `src/MZMfitting.py` searches for MZM XML files, creates analysis figures, and writes wafer-level and timestamp-level CSV summaries.
---

### :warning:precautions

 1) Raw XML files should follow the existing wafer/timestamp folder structure under `data`.
 2) Files without valid MZM wavelength sweep data can be skipped during analysis.
 3) Existing result files under the old flat output layout may be cleaned before new SPDAP outputs are written.
