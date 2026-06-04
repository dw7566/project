![header](https://capsule-render.vercel.app/api?type=wave&color=auto&height=200&section=header&text=GPDO%20Wafer%20Analyzer&fontSize=70)

###### Germanium Photodetector Data Analysis Automation Pipeline
### Contents

[1. Introduction](#1-introduction)\
[2. Project information](#2-project-information)\
[3. Install and Run](#3-install-and-run)\
[4. Description of the module file feature](#4-description-of-the-module-file-feature)\
[5. Run file algorithm](#5-run-file-algorithm)

---

# GPDO Wafer Analyzer :
##### Hi !
##### Thank you for looking at our project. GPDO Wafer Analyzer automates **Germanium Photodetector XML 측정 데이터를 자동으로 파싱·피팅·시각화하는 웨이퍼 분석 파이프라인** 입니다.
##### This project automates wafer-scale GPDO measurement analysis and generates organized CSV and PNG results.

---

## 1. Introduction
We aim to develop a Python-based automation pipeline for GPDO data analysis.
광소자 공정 연구에서 웨이퍼 한 장에는 수십 개의 GPDO 다이가 있고, 다이마다 XML 측정 파일이 생성됩니다.
이 파일들을 수작업으로 열어 그래프를 그리고 파라미터를 뽑아내는 작업은 시간이 오래 걸리고 실수가 생기기 쉽습니다.

The goal is to process wafer-scale XML measurement data, extract GPDO device information including Dark/Light/Spectrum current data, and generate analysis outputs with key parameters.

#### - Main Features
- **파싱**: GPDO XML에서 Dark/Light/Spectrum 전류, Reference IL 추출
- **피팅**: Shockley 다이오드 모델, Power-law 역바이어스, 광전류 계산
- **시각화**: 다이별 6-패널 PNG + 웨이퍼 전체 히트맵
- **CSV 출력**: 다이별 핵심 파라미터(Iph, n, R, peak λ 등) 정리

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
    + Wafer (웨이퍼)
    + Die row & column (다이 위치)
    + GPDO XML measurement files (GPDO XML 측정 파일)
    + Dark/Light/Spectrum current data (암전류/명전류/스펙트럼 데이터)
    + Optical and electrical characteristics (광학·전기 특성)

####
+ **Run file description**

   GPDO Wafer Analyzer scans the `data` directory and extracts XML files whose names include `GPDO`.\
   Then, it parses dark current, light current, spectrum, and voltage data, runs fitting/analysis logic, and saves CSV tables and PNG figures under the `res` directory.

####
+ **Output Parameters**

| 파라미터 | 기호 | 설명 |
|----------|------|------|
| 광전류 | `Iph` | 역바이어스(-1.5V 이하)에서 Light − Dark 차감으로 추출한 포토커런트 |
| 이상계수 | `n_d` | Shockley 다이오드 이상계수. 1에 가까울수록 이상적인 pn 접합 |
| 응답도 | `R_resp` | 단위 광파워당 전류 [A/W]. GPDO 성능의 핵심 지표 |
| 측정 파장 | `lc_wl` | Light Current 측정에 사용된 단일 파장 [nm] |
| 스펙트럼 피크 파장 | `peak_wl` | 파장 스윕에서 전류가 최대인 파장 [nm] |
| 순방향 R² | `r2_fwd` | 순방향 Dark Current 피팅 결정계수 (1에 가까울수록 피팅 품질 우수) |
| 광전류 R² | `r2_photo` | 역바이어스 구간 광전류 포화 균일도 (1에 가까울수록 안정적인 포화) |

---

## 3. Install and Run

####
* Getting Started
   + Enter the Terminal and install required packages. \
```bash
pip install numpy scipy matplotlib lxml pandas
```

* How to Run

  + **Data Preparation**: Place the raw XML measurement data under the `data` directory.

```
data/
└── HY202103/
    ├── D08/
    │   └── 20190526_082853/
    │       ├── HY202103_D08_(-1,-1)_LION1_DCM_GPDO.xml
    │       └── ...
    └── D24/
```

  + **Run the main script**. GPDO Wafer Analyzer will analyze available GPDO XML files and save results automatically.

```bash
# 전체 처리
python run.py

# GPDO만 처리
python run.py GPDO

# 특정 디바이스 여러 개 선택
python run.py GPDO LMZC
```

  + **Results**: CSV outputs are saved in `res/csv`, and generated analysis figures are saved in `res/png` and `res/heatmap`.

```
res/
├── D08-GPDO/
│   └── 20190526_082853/
│       ├── png/                     # 다이별 6-패널 분석 그래프
│       │   ├── HY202103_D08_(-1,-1)_LION1_DCM_GPDO.png
│       │   └── ...
│       └── heatmap/                 # 웨이퍼 히트맵
│           ├── heatmap_R_resp.png
│           └── heatmap_n_d.png
└── csv/
    ├── D08_Result.csv               # 웨이퍼별 CSV
    └── Total_Result.csv             # 전체 통합 CSV
```

---

## 4. Description of the module file feature

* **Fitting module** (`src/fitting/fitting_engine.py`)
   + The analysis is performed by parsing dark current, light current, spectrum, and voltage data from XML files.
   + The module implements Shockley diode model for forward bias, Power-law model for reverse bias, and photo-current extraction.
   + Estimates fitting quality such as R-squared and visualizes IV behavior and optical response.

* **Parser module** (`src/parser/gpdo_parser.py`)
  + Extracts Dark/Light/Spectrum current, Reference IL, wavelength, fiber power from GPDO XML files.
  + Converts raw measurement data into structured numpy arrays for further analysis.

* **Plotter module** (`src/plotting/plotter.py` & `heatmap_plotter.py`)
  + Generates 6-panel analysis figures for each die: Reference Spectrum, Dark I–V, Light I–V, Photo Current, Spectrum, and Responsivity.
  + Creates wafer-level heatmaps showing spatial distribution of Responsivity and ideality factor.

* **CSV module** (`src/tocsv/gpdo_csv.py`)
  + Collects and organizes measurement information including wafer ID, die position, timestamp, device name, bias, current values, wavelength, and fitting parameters.
  + Creates a dataframe-style CSV summary for quick data review.
  + Saves results in CSV format in the `res/csv` folder.

---

## 5. Run file algorithm

* **Preparation**
   + Configure `config.py`: Set PROJECT_NAME, WAFER_IDS, DEVICE_CONFIG
   + Place the measurement XML files under the `data` directory.

* **Execution**
   + The `main` function in `run.py` is executed.
   + `src/analyzer/gpdo_analyzer.py` searches for GPDO XML files, creates analysis figures, and writes wafer-level and timestamp-level CSV summaries.
   + For each XML file: parse data → fit parameters → generate 6-panel plot → save to CSV

* **Output**
   + Analysis results are generated automatically and saved in the `res` directory.
   + CSV summaries are stored in `res/csv` folder.
   + Generated die-level figures are stored in `res/{wafer_id}-GPDO/{timestamp}/png/` folder.
   + Generated wafer-level heatmaps are stored in `res/{wafer_id}-GPDO/{timestamp}/heatmap/` folder.

---

## Project Structure

```
project/
├── run.py                      # 실행 진입점 (Main execution script)
├── config.py                   # 경로·웨이퍼·디바이스 설정 (Configuration)
├── data/                       # 원본 XML (Input directory for XML measurement files)
│   └── HY202103/
│       ├── D08/
│       └── D24/
├── res/                        # 결과 저장 (Output directory for results)
│   ├── csv/                    # CSV 분석 결과
│   ├── D08-GPDO/
│   │   └── {timestamp}/
│   │       ├── png/            # 다이별 6-패널 분석 그래프
│   │       └── heatmap/        # 웨이퍼 히트맵
│   └── D24-GPDO/
└── src/
    ├── parser/
    │   └── gpdo_parser.py      # XML 파싱 (XML Parser)
    ├── fitting/
    │   └── fitting_engine.py   # 피팅 모델 + 연산 (Fitting Models)
    ├── plotting/
    │   ├── plotter.py          # 다이 6-패널 PNG (Die-level Plotter)
    │   └── heatmap_plotter.py  # 웨이퍼 히트맵 PNG (Heatmap Plotter)
    ├── analyzer/
    │   └── gpdo_analyzer.py    # 전체 파이프라인 통합 (Pipeline Orchestrator)
    └── tocsv/
        └── gpdo_csv.py         # CSV 저장 (CSV Export)
```

---

## Configuration

### 웨이퍼 추가 / 제거 (Add/Remove Wafers)

```python
# config.py
WAFER_IDS = ["D07", "D08", "D23", "D24"]
```

### 프로젝트 데이터셋 변경 (Change Project Dataset)

```python
# config.py
PROJECT_NAME = "HY202103"   # data/ 바로 아래 폴더명과 일치
```

### 특정 디바이스에 다른 웨이퍼 셋 지정 (Configure Device-Specific Wafers)

```python
# config.py
DEVICE_CONFIG = {
    "GPDO": dict(
        wafer_ids = ["D08", "D24"],   # 이 디바이스만 D08·D24 처리
        save_root = "GPDO",
    ),
}
```

---

## Requirements

This project requires Python 3.10+ with the following dependencies:

- `numpy` - Numerical computing
- `scipy` - Scientific computing and fitting
- `matplotlib` - Data visualization
- `lxml` - XML parsing
- `pandas` - Data analysis and CSV handling

---

## License

This project is developed by the Silicon Photonics Research Team at Hanyang University.

---

**Last Updated**: 2026-06-04
