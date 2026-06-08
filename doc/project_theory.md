# SPDAP 프로젝트 설명 및 이론 문서

## 1. 프로젝트 개요

SPDAP는 Silicon Photonics Data Analysis Automation Pipeline의 약자이다. 이 프로젝트는 실리콘 포토닉스 기반 MZM, 즉 Mach-Zehnder Modulator 소자의 웨이퍼 단위 XML 측정 데이터를 자동으로 읽고, 광학 스펙트럼과 전기적 IV 특성을 분석한 뒤 CSV와 PNG 결과물로 정리하는 파이프라인이다.

기존에는 여러 die와 bias 조건에서 측정된 XML 파일을 사람이 직접 열어보고, 필요한 그래프를 만들고, 주요 성능 지표를 계산해야 했다. SPDAP는 이 과정을 자동화하여 반복 작업을 줄이고, die 단위 분석과 wafer 단위 비교를 일관된 방식으로 수행할 수 있게 한다.

이 프로젝트의 주요 분석 대상은 `data/` 디렉터리 아래에 있는 `*LMZ*.xml` 파일이다. 각 XML 파일에는 MZM 소자의 wavelength sweep 데이터, insertion loss 데이터, DC bias 정보, IV 측정 데이터가 포함된다. 분석 결과는 `res/csv/`와 `res/png/` 아래에 저장된다.

## 2. 실행 방법

필요한 Python 패키지를 설치한 뒤 `run.py`를 실행한다.

```bash
pip install -r requirements.txt
python run.py
```

`run.py`는 내부적으로 `src.main.main()`을 호출한다. 메인 함수는 `data/` 폴더를 재귀적으로 탐색하여 파일명에 `LMZ`가 포함된 XML 파일을 찾고, 각 파일에 대해 분석과 시각화를 수행한다.

## 3. 입력 및 출력 구조

입력 데이터는 일반적으로 다음과 같은 계층 구조를 가진다.

```text
data/
  HY202103/
    D08/
      20190712_113254/
        HY202103_D08_(0,0)_LION1_DCM_LMZC.xml
        ...
```

파일 경로에는 lot, wafer, timestamp, die 좌표, test site, device type 정보가 포함된다. 프로젝트는 XML 내부의 `TestSiteInfo`와 파일명 정보를 함께 사용하여 분석 결과에 metadata를 기록한다.

출력은 크게 세 가지이다.

- `res/png/{wafer}/{timestamp}/*.png`: die 단위 3x3 분석 그래프
- `res/png/{wafer}/{timestamp}/wafermap.png`: wafer 단위 성능 분포 그래프
- `res/csv/{wafer}/{timestamp}.csv`, `res/csv/mzm_all_summary.csv`: 분석 지표 요약 CSV

현재 코드 기준 CSV의 기본 컬럼은 다음과 같다.

| 컬럼 | 의미 |
|---|---|
| `lot` | 측정 lot 또는 batch ID |
| `wafer` | wafer ID |
| `test_site` | test site 또는 device group |
| `die_column`, `die_row` | wafer 내 die 좌표 |
| `timestamp` | 측정 시간 폴더명 |
| `device_name` | XML에 기록된 device 이름 |
| `dc_bias_v` | wavelength sweep에 사용된 DC bias |
| `current_at_minus_2v_a` | -2 V 근처의 전류 |
| `current_at_minus_1v_a` | -1 V 근처의 전류 |
| `current_at_0v_a` | 0 V 근처의 전류 |
| `current_at_plus_1v_a` | +1 V 근처의 전류 |
| `extinction_ratio_db` | sweep 내 최대 IL과 최소 IL의 차이 |
| `vpi_at_dc_bias_v` | 해당 bias에서 추출한 Vpi |
| `source_file` | 원본 XML 파일 경로 |

## 4. 전체 분석 파이프라인

SPDAP의 전체 처리 흐름은 다음과 같다.

1. `data/`에서 `*LMZ*.xml` 파일을 찾는다.
2. XML에서 wavelength sweep, insertion loss, DC bias, IV 데이터를 파싱한다.
3. reference sweep을 3차 다항식으로 fitting하여 baseline을 추정한다.
4. 각 bias sweep에서 baseline과 envelope trend를 제거하여 spectrum을 정규화한다.
5. MZI 모델 fitting을 통해 FSR과 fitting 품질을 확인한다.
6. IV curve를 로그 스케일로 표시하고 forward/reverse 영역의 fitting 특성을 확인한다.
7. null 위치 변화로부터 modulation efficiency와 Vpi를 계산한다.
8. peak와 null의 차이로 extinction ratio를 계산한다.
9. die 단위 3x3 PNG와 wafer 단위 wafermap을 생성한다.
10. bias별, die별 분석 결과를 CSV로 저장한다.

## 5. MZM 이론 배경

MZM은 Mach-Zehnder Interferometer 구조를 이용한 광 변조기이다. 입력 광은 두 개의 arm으로 나뉘고, 각 arm을 지나면서 서로 다른 위상 변화를 겪은 뒤 다시 합쳐진다. 두 경로의 위상 차이에 따라 출력 광은 보강 간섭 또는 상쇄 간섭을 일으킨다.

두 arm의 위상 차이를 `Delta phi`라고 하면 이상적인 출력 세기는 다음과 같은 형태로 표현할 수 있다.

```text
T = A + B * cos^2(Delta phi / 2)
```

이 프로젝트의 `src/spectrum.py`에서는 wavelength에 따른 MZI 응답을 다음과 같은 모델로 fitting한다.

```text
T(lambda) = A + slope * x + B * cos^2(pi * (lambda - lambda0) / FSR + phi)
```

여기서 `lambda`는 파장, `FSR`은 Free Spectral Range, `phi`는 위상 offset, `A`와 `B`는 transmission offset과 modulation amplitude를 의미한다. `slope` 항은 파장에 따른 완만한 기울기를 보정하기 위해 사용된다.

## 6. Insertion Loss와 Transmission Spectrum

Insertion Loss, 즉 IL은 소자를 통과하면서 광 신호가 얼마나 손실되었는지를 dB 단위로 나타내는 값이다. XML의 wavelength sweep에는 파장별 IL 값이 저장되어 있으며, 각 DC bias 조건에서 서로 다른 spectrum이 나타난다.

노트북의 raw transmission spectra 그래프는 이 측정값을 그대로 표시한다. 이 그래프를 보면 bias 변화에 따라 간섭 fringe와 null 위치가 이동하는지 확인할 수 있다. 다만 raw spectrum에는 광원 세기 변화, coupling loss, 측정 setup의 wavelength-dependent trend가 함께 포함될 수 있기 때문에 바로 정량 분석에 사용하기 어렵다.

이를 보정하기 위해 reference sweep을 3차 다항식으로 fitting하고, 각 sweep에서 이 reference baseline을 제거한다. 이후 envelope flattening을 적용해 fringe의 상대적인 변화가 더 잘 보이도록 만든다.

## 7. FSR 이론

FSR, Free Spectral Range는 인접한 peak 또는 null 사이의 파장 간격을 의미한다. MZI 구조에서는 두 arm의 optical path length 차이 때문에 주기적인 간섭 spectrum이 나타나며, 그 주기가 FSR이다.

일반적으로 FSR은 다음 관계와 연결된다.

```text
FSR ~= lambda^2 / (n_g * Delta L)
```

여기서 `lambda`는 중심 파장, `n_g`는 group index, `Delta L`은 두 arm의 길이 차이이다. 즉 arm 길이 차이가 클수록 FSR은 작아지고, 길이 차이가 작을수록 FSR은 커진다.

코드에서는 `scipy.signal.find_peaks`를 이용해 spectrum의 notch 위치를 찾고, 인접 notch 사이의 median spacing을 FSR로 추정한다. notch 검출이 충분하지 않을 경우 device type에 따라 fallback 값을 사용한다. 예를 들어 `LMZC` 계열은 14 nm, `LMZO` 계열은 10 nm를 fallback으로 사용한다.

## 8. Reference Fit과 Flattening

측정 spectrum에는 MZM 자체의 간섭 특성뿐 아니라 optical fiber coupling, grating coupler, 장비 응답, 광원 power variation 등이 함께 반영된다. 이런 완만한 trend를 제거해야 실제 modulation fringe를 안정적으로 분석할 수 있다.

SPDAP는 reference sweep을 3차 다항식으로 fitting한다.

```text
baseline(lambda) = p3 * lambda^3 + p2 * lambda^2 + p1 * lambda + p0
```

그 다음 각 sweep에서 baseline을 빼고, peak envelope을 추정하여 spectrum을 flatten한다. 이 과정의 목적은 절대 power level보다 peak-null 구조와 bias에 따른 fringe 이동을 더 명확하게 보기 위한 것이다.

## 9. Extinction Ratio 이론

Extinction Ratio, ER은 optical modulator가 켜진 상태와 꺼진 상태를 얼마나 잘 구분하는지를 나타내는 지표이다. dB 단위에서는 일반적으로 peak transmission과 null transmission의 차이로 계산할 수 있다.

```text
ER [dB] = T_peak [dB] - T_null [dB]
```

ER이 클수록 on/off contrast가 크다는 뜻이며, 변조기의 신호 구분 능력이 좋다고 해석할 수 있다. 이 프로젝트에서는 각 bias sweep에서 local maximum과 local minimum을 찾고, null 주변의 peak와 비교하여 ER을 계산한다. `src/extinction_ratio.py`에서는 10 dB 미만의 약한 null은 유효한 fringe로 보지 않고 제외한다.

노트북의 `Extinction Ratio versus DC Bias` 그래프는 bias별 ER의 평균, 최소, 최대를 함께 보여준다. 이를 통해 특정 bias에서 소자의 contrast가 좋아지는지, bias 변화에 따라 ER이 안정적인지 확인할 수 있다.

## 10. Vpi 이론

Vpi는 half-wave voltage를 의미한다. MZM에서 출력 상태를 최대 transmission에서 최소 transmission으로 바꾸려면 두 arm 사이에 `pi`만큼의 위상 차이를 만들어야 한다. 이때 필요한 전압을 Vpi라고 한다.

Vpi가 작을수록 더 낮은 전압으로 같은 광 변조를 만들 수 있으므로, 일반적으로 modulation efficiency가 좋다고 해석한다.

MZM spectrum에서 bias가 변하면 null wavelength가 이동한다. 이 이동량을 `d lambda / dV`로 나타내면, FSR과 연결하여 Vpi를 추정할 수 있다.

```text
Vpi = FSR / (2 * |d lambda / dV|)
```

여기서 `d lambda / dV`는 bias 변화에 따른 null wavelength 이동 기울기이다. FSR의 절반만큼 null이 이동하면 위상은 pi만큼 변한 것으로 볼 수 있기 때문에 위 식을 사용한다.

`src/vpi_analysis.py`는 각 bias 조건에서 깊은 null 위치를 찾고, 같은 null이 bias에 따라 어떻게 이동하는지 track한다. 모든 bias에서 추적 가능한 null에 대해 wavelength와 voltage의 관계를 fitting하고, 이 기울기를 이용해 Vpi를 계산한다. track의 품질이 너무 낮거나 기울기가 너무 작으면 유효하지 않은 Vpi로 판단하고 제외한다.

## 11. IV Curve 이론

IV curve는 전압에 따른 전류 특성을 나타낸다. MZM의 phase shifter 또는 PN junction 구조에서는 전압에 따라 leakage current, forward current, reverse current 특성이 달라진다.

SPDAP는 IV 데이터를 두 가지 방식으로 표시한다.

- 로그 스케일 IV 그래프: 전류 범위가 매우 넓기 때문에 `abs(I)`를 log scale로 표시한다.
- IV analysis 그래프: reverse 영역은 다항식으로 fitting하고, forward turn-on이 보이면 diode equation으로 fitting한다.

forward diode fitting에는 다음 형태가 사용된다.

```text
I = Is * (exp(V / (n * Vt)) - 1)
```

여기서 `Is`는 saturation current, `n`은 ideality factor, `Vt`는 thermal voltage이다. 코드에서는 `Vt = 0.02585 V`를 사용한다. IV 분석은 광학 성능과 별개로 소자의 전기적 leakage, diode turn-on, 측정 이상 여부를 확인하는 데 사용된다.

## 12. Die-Level 3x3 그래프 해석

각 XML 파일에 대해 생성되는 die-level PNG는 3x3 패널로 구성된다.

| 위치 | 그래프 | 해석 |
|---|---|---|
| (0,0) | Raw transmission spectra | bias별 원본 IL spectrum |
| (0,1) | Reference fit | reference sweep과 3차 다항식 baseline |
| (0,2) | Flattened spectra | baseline과 envelope을 제거한 정규화 spectrum |
| (1,0) | MZM linear fit | linear transmission에서 MZI model fitting |
| (1,1) | IV log scale | 전압-전류 특성의 로그 표현 |
| (1,2) | IV analysis | reverse/forward fitting과 fitting 품질 |
| (2,0) | MZM dB residual fit | dB domain에서 residual 보정 후 MZI fitting |
| (2,1) | Vpi vs DC bias | bias별 half-wave voltage |
| (2,2) | ER vs DC bias | bias별 extinction ratio |

이 그래프는 한 die의 광학적 성능과 전기적 특성을 한 번에 확인하기 위한 보고서형 결과물이다.

## 13. Wafermap 해석

Wafermap은 여러 die의 분석 결과를 wafer 좌표계 위에 표시한 summary figure이다. 현재 프로젝트의 wafermap은 die 위치별 extinction ratio와 bias별 평균 ER 경향을 보여준다.

die 위치별 ER 분포를 보면 wafer 중심과 edge 사이의 성능 차이, 특정 영역의 공정 불균일, 비정상 die를 빠르게 확인할 수 있다. bias별 평균 ER 그래프는 전체 wafer에서 어떤 DC bias 조건이 평균적으로 좋은 contrast를 만드는지 판단하는 데 사용된다.

## 14. 코드 모듈 역할

| 파일 | 역할 |
|---|---|
| `run.py` | 프로젝트 실행 진입점 |
| `src/main.py` | 전체 파이프라인 제어, PNG/CSV 생성 |
| `src/xml_parser.py` | XML parsing, metadata와 sweep/IV 데이터 추출 |
| `src/spectrum.py` | MZI model, FSR 측정, envelope flattening |
| `src/iv_analysis.py` | IV log plot 및 diode/reverse fitting |
| `src/extinction_ratio.py` | ER 계산 및 bias별 ER 그래프 |
| `src/vpi_analysis.py` | null tracking 기반 modulation efficiency와 Vpi 계산 |
| `src/csv_export.py` | XML별 요약 row 생성 및 CSV 저장 |
| `src/wafermap.py` | wafer 단위 summary map 생성 |
| `src/config.py` | 입력/출력 경로, thermal voltage, modulation bias, CSV 컬럼 설정 |

## 15. 프로젝트 결과의 의미

이 프로젝트의 결과는 단순히 그래프를 자동 생성하는 것에 그치지 않는다. 각 die의 optical spectrum, IV 특성, ER, Vpi를 같은 방식으로 계산하므로 die 간 비교가 가능해진다. 또한 wafermap을 통해 wafer 전체의 공간적 성능 분포를 확인할 수 있어 공정 균일성 분석에도 활용할 수 있다.

예를 들어 ER이 높은 die는 optical contrast가 좋다는 뜻이고, Vpi가 낮은 die는 낮은 전압으로 효율적인 변조가 가능하다는 뜻이다. IV curve에서 reverse leakage가 크거나 forward fitting이 비정상적이면 전기적 결함 또는 측정 이상 가능성을 의심할 수 있다.

## 16. 주의점 및 한계

자동 분석은 일관성과 속도 면에서 장점이 있지만, 모든 측정 데이터에 대해 완벽한 해석을 보장하지는 않는다. 다음 사항을 고려해야 한다.

- XML 구조가 예상과 다르면 일부 데이터가 파싱되지 않을 수 있다.
- spectrum의 noise가 크거나 null이 충분히 깊지 않으면 FSR, ER, Vpi 추정이 불안정할 수 있다.
- Vpi는 null tracking과 FSR 추정에 의존하므로, null 위치가 bias 전체에서 안정적으로 추적되어야 한다.
- reference sweep이 실제 baseline을 대표하지 못하면 flattening 결과가 왜곡될 수 있다.
- CSV의 `extinction_ratio_db`는 현재 코드 기준으로 sweep 내 `IL_max - IL_min` 방식의 요약값이며, 그래프의 ER 분석은 local peak-null pair 기반 분석이다.

따라서 최종 판단에서는 CSV 값뿐 아니라 die-level PNG와 wafermap을 함께 확인하는 것이 필요하다.

## 17. 요약

SPDAP는 wafer-scale MZM XML 측정 데이터를 자동으로 분석하는 Python 기반 파이프라인이다. 이 프로젝트는 XML parsing, reference fitting, spectrum flattening, MZI model fitting, FSR 추정, ER 계산, Vpi 추출, IV curve 분석, die-level graph 생성, wafermap 생성, CSV 저장까지 하나의 흐름으로 처리한다.

핵심 이론은 Mach-Zehnder 간섭, FSR, extinction ratio, half-wave voltage, diode IV 특성이다. 최종 결과를 통해 개별 die의 성능뿐 아니라 wafer 전체의 공정 균일성과 bias 조건별 성능 변화를 한눈에 평가할 수 있다.
