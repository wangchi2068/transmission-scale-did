# Data Contract: German County-Level Weekly Varicella Panel (2015--2019)

This directory contains the public raw datasets, automated download script (`fetch_germany.py`), and deterministic processing pipeline for the empirical stress test in the manuscript (German county-level weekly Varicella surveillance counts and staggered 16-state school holidays, 2015--2019, with three pre-sample lag weeks from late 2014).

## Download & Line-Ending Stability

- Run `python experiment/data/fetch_germany.py` to download all three raw files directly from the public RKI SurvStat@RKI 2.0 SOAP 1.2 WebService, OpenHolidays API, and GitHub county GeoJSON repository, or `python experiment/data/fetch_germany.py --check` to verify local files.
- `.gitattributes` marks `experiment/data/raw/**` as `binary`, ensuring Git preserves exact LF (`\n`) bytes and identical SHA-256 checksums across Linux, macOS, and Windows checkouts.

## Raw Inputs (`experiment/data/raw/germany/`)

Archived SurvStat extraction date: **`2026-11-19`** (OpenHolidays and AGS GeoJSON fetched same day; persisted at `https://github.com/wangchi2068/transmission-scale-did`). Run `python experiment/data/fetch_germany.py` to download fresh upstream copies canonically (volatile fetch-date metadata is written separately to `fetch_metadata.json` so canonical re-downloads match the archived LF SHA-256 checksums whenever upstream records have not changed) or `--check` to verify the exact `2026-11-19` LF archive.

| File | Source & Query Specification | Access / Licence | Bytes (LF) | SHA-256 |
| --- | --- | --- | ---: | --- |
| `rki_varicella_2014_2019.json` | RKI SurvStat@RKI 2.0 SOAP 1.2 (`https://tools.rki.de/SurvStat/SurvStatWebService.svc`), extracted `2026-11-19`, Cube `SurvStat`, Pathogen `[PathogenOut].[KategorieNz].[Krankheit DE].&[Windpocken]`, Column `[ReportingDate].[Week].[Week]` (*Meldewoche*, notification week; symptom-onset week is not exposed by the public SurvStat web service), Row `[DeutschlandNodes].[Kreise71Web].[CountyKey71]`, Measure `[Measures].[FallCount_71_Web]` (all notified cases without `[ReferenzDefinition]` restriction) and `[Measures].[Inzidenz_71_Web]` for 2017 reference populations (`Count / Incidence * 1e5`), plus 2015/2018 RKI incidence-denominator fallback for the three counties reporting zero cases in 2017 (`09188`, `09763`, `16069`). | Public RKI SurvStat@RKI 2.0 open surveillance data (Robert Koch Institute) | 2,475,214 | `0735f0211c4d51389a20cb93dbcd4d5b48f8f7790d1c68176c8841c12d0b2489` |
| `germany_school_holidays_2014_2019.json` | OpenHolidays API (`https://openholidaysapi.org/SchoolHolidays`), extracted `2026-11-19`, official Kultusministerkonferenz (KMK) German federal-state school holiday schedules (Summer, Easter, Autumn, Winter, Pentecost, Christmas) for all 16 Länder (`DE-SH` through `DE-TH`), 2014--2019. | Open Data Commons Open Database License (ODbL-1.0), `openpotato/openholidaysapi.data` (attribute OpenHolidays / KMK) | 261,944 | `4f51fd93e04f0f7ae37a5d78aea0966d12f74acf15b7293a795f155826fb695e` |
| `counties_ags.geojson` | Official German 5-digit Amtlicher Gemeindeschlüssel (AGS) land county polygons (`gf=4`, 400 features) from BKG VG250 WFS (`https://sgx.geodatenzentrum.de/wfs_vg250`, layer `vg250:vg250_krs`), fetched `2026-11-20`. | Data licence Germany -- attribution -- Version 2.0 (`dl-de/by-2-0`), © GeoBasis-DE / BKG | 9,285,129 | `5003586e2c43f803828b106bb7de11b58079529cb78bf76649e3d2ec703acacc` |

## Surveillance Context & Processing Rules (`experiment/code/data_prep.py`)

1. **Notification regime & date dimension**:
   - Nationwide mandatory physician/laboratory reporting of varicella (*Windpocken*) under the German *Infektionsschutzgesetz* took effect on 29 March 2013 (replacing prior state-only mandates in several eastern Länder); our estimation window begins in January 2015 (`2015-W01`), less than two years after national rollout, which motivates checking 2016--2019 subsamples alongside state linear trends.
   - Weekly counts are indexed by RKI notification week (`ReportingDate` / *Meldewoche*), not symptom-onset date, which is unavailable in the public SurvStat web service. Notification delays combine with the ~14-day serial interval when interpreting event-time profiles.
2. **Spatial units ($N = 400$ counties across $G = 16$ federal states, post-July 2021 RKI Kreise classification)**:
   - RKI reports 412 district keys including the 12 administrative boroughs of Berlin (`11001`--`11012`). Berlin's 12 boroughs are summed into city-state county unit `11000`, yielding the standard 400-county partition.
   - Administrative mergers reconciling the earlier GeoJSON boundary vintage to the post-July 2021 RKI `Kreise71Web` key set: Göttingen (`03152` + `03156` $\to$ `03159`, effective 1 November 2016), Wartburgkreis / Eisenach (`16056` $\to$ `16063`, effective 1 July 2021), and Berlin boroughs (`11*` $\to$ `11000`). All 400 counties match with zero isolated units.
   - Two counties share an undirected contiguity edge if their boundary polygons share at least one boundary vertex rounded to $10^{-3}$ degrees ($\approx 100$ m), yielding 1,044 undirected county edges (mean degree 5.220; compared with the earlier 1,043-edge graph, official BKG VG250 includes the shared boundary across the lower Elbe between Pinneberg `01056` and Stade `03359`) and 29 state-contiguity edges (diameter 5). Weight matrix $W$ is row-normalised.
3. **Temporal window ($T = 261$ estimation weeks + $S = 3$ pre-sample lag weeks)**:
   - Estimation weeks: ISO weeks `2015-W01` through `2019-W52` ($T = 261$; 2015 has 53 ISO weeks). Pre-sample lag weeks: `2014-W50` through `2014-W52`.
4. **Policy & placebo exposures**:
   - Let $f_{it} \in [0, 1]$ be the fraction of ISO week $t$ (Monday--Sunday) covered by state $s(i)$'s official summer holiday (*Sommerferien*). Primary binary exposure sets $D^{(0)}_{it} = \mathbf{1}\{7 f_{it} \ge 4\}$ lagged one week ($D_{it} = D^{(0)}_{i,t-1}$; 11.8% treated county-weeks). Threshold variants ($\ge 3$, $\ge 5$ days), continuous $f_{it}$, lags 0--4, all school holidays, and short-holiday placebos (Easter, Autumn, Winter, Pentecost) are built identically and stored in `panel.npz`.
   - Note on summer rotation: 14 states rotate summer start dates across ISO weeks 25--31 over 2015--2019, whereas Baden-Württemberg (`08`, 44 counties) and Bavaria (`09`, 96 counties) always begin last (ranks 15--16, ISO weeks 30--31) in all five years.
5. **Background importation ($\omega_i$) & regional groupings**:
   - $\omega_i = \text{population}_i / 10^5$ using official RKI 2017 county denominators (and official RKI 2015/2018 denominators for `09188`, `09763`, `16069`; 0 median-imputed counties, national total 82.89 million).
   - 4 German macro-regions (North: SH, HH, NI, HB, MV; West: NW, HE, RP, SL; South: BW, BY; East: BE, BB, SN, ST, TH) and 38 NUTS-2 government districts (*Regierungsbezirke*).
