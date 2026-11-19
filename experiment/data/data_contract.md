# Data Contract: German County-Level Weekly Varicella Panel (2015--2019)

This directory contains the public raw datasets and deterministic processing pipeline for the empirical application in Section 4 of the manuscript (German county-level weekly Varicella surveillance counts and staggered 16-state summer school holidays, 2015--2019, with three pre-sample lag weeks from late 2014).

## Raw Inputs (`experiment/data/raw/germany/`)

| File | Source | License / Access | SHA-256 |
| --- | --- | --- | --- |
| `rki_varicella_2014_2019.json` | Robert Koch-Institut (RKI) SurvStat@RKI 2.0 SOAP service (`https://tools.rki.de/SurvStat/SurvStatWebService.svc`), pathogen `[PathogenOut].[KategorieNz].[Krankheit DE].&[Windpocken]`, weekly county counts 2014--2019 and 2017 RKI county reference populations (`Count / Incidence * 1e5`) | Open public health surveillance data (RKI) | `b55d8973975aba56de69b2f06adc4174e1f136432115efe87a94c06ba4f45d01` |
| `germany_school_holidays_2014_2019.json` | OpenHolidays API (`https://openholidaysapi.org/SchoolHolidays`), official Kultusministerkonferenz (KMK) German federal-state school holiday schedules for all 16 Länder (`DE-BW` through `DE-TH`), 2014--2019 | Open public schedule data | `fc799aa000bd362105c95e2a722ac2a5279abc8d3767b2cf18a97cd721bd008a` |
| `counties_ags.geojson` | German administrative county polygons keyed by 5-digit Amtlicher Gemeindeschlüssel (AGS), `m-ad/geofeatures-ags-germany` (`geojson/counties.json`) | Open geodata (BKG / GeoBasis-DE derived) | `8a0d0bf29ff31b6e536622613ce6e2f6198a6999bd725a631b4b1707780a6abd` |

## Processing Rules (`experiment/code/data_prep.py`)

1. **Spatial units ($N = 400$ counties across $G = 16$ federal states)**:
   - RKI reports 412 district keys including the 12 administrative boroughs of Berlin (`11001`--`11012`). Berlin's 12 boroughs are summed into a single city-state county unit (`11000`), matching the standard 400-county German Kreise/kreisfreie Städte partition.
   - Known county boundary consolidations between the GeoJSON vintage and RKI 2017+ AGS codes are reconciled explicitly: Göttingen (`03152` + `03156` $\to$ `03159`), Wartburgkreis/Eisenach (`16056` $\to$ `16063`), and Berlin boroughs (`11*` $\to$ `11000`). All 400 counties match with 0 isolated units.
   - Two counties share an undirected contiguity edge if their boundary polygons share at least one boundary vertex rounded to $10^{-3}$ degrees (~100 m), yielding 1,043 undirected edges (mean degree 5.21). The spatial weight matrix $W$ is row-normalised.
2. **Temporal window ($T = 261$ estimation weeks + $S = 3$ pre-sample lag weeks)**:
   - Estimation weeks: ISO weeks `2015-W01` through `2019-W52` (261 weeks; 2015 has 53 ISO weeks).
   - Pre-sample weeks for constructing the renewal pressure at `2015-W01` ($S = 3$): `2014-W50`, `2014-W51`, `2014-W52`.
3. **Policy exposure ($D_{it}$)**:
   - For each county $i$ in federal state $s(i)$ and ISO week $t$ (Monday--Sunday), let $f_{it} \in [0, 1]$ be the fraction of calendar days covered by state $s(i)$'s official summer school holiday (*Sommerferien*).
   - Binary indicator $D^{(0)}_{it} = \mathbf{1}\{7 f_{it} \ge 4\}$ (majority of the week in summer holiday).
   - Primary exposure aligns with Varicella's ~14-day incubation and reporting lag via a 1-week lag: $D_{it} = D^{(0)}_{i,t-1}$ (with lags 0, 2, 3, 4, continuous fraction $f_{i,t-1}$, and full event-time profiles around holiday onset $e \in \{-2, \dots, 7\}$ evaluated directly in the sensitivity suite).
4. **Background importation ($\omega_i$)**:
   - Population-proportional background $\omega_i = \text{pop}_i / 10^5$ (with $\times 0.25$ and $\times 4$ sensitivity). Four counties with zero 2017 RKI count (preventing direct `Count / Incidence` population inversion) are assigned the median county population.
5. **Derived spatial groupings**:
   - 4 German macro-regions (North: SH, HH, NI, HB, MV; West: NW, HE, RP, SL; South: BW, BY; East: BE, BB, SN, ST, TH) for region $\times$ week fixed-effects sensitivity.
   - 38 NUTS-2 government districts (*Regierungsbezirke* / state-level equivalents from the first two digits of AGS or NUTS-2 prefix) for low-zero spatial aggregation sensitivity.

## Output (`experiment/data/processed/`)

- `panel.npz`: compressed NumPy archive containing `Y`, `D_summer`, `Dshare_summer`, `D_all`, `Dshare_all`, `pop`, `state`, `cl`, `states`, `region`, `nuts2`, `ags`, `weeks`, `edges`, `W_indptr`, `W_indices`, `W_data`.
- `panel_summary.json`: machine-readable summary statistics and checksum verification record.
