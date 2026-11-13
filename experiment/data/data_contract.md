# Data contract — application (US county mask orders, June–November 2020)

Scope: the single empirical illustration pre-specified in the method brief (state/territorial public mask
orders and weekly county COVID-19 cases). No other outcome, population or research question is added.
All files are public; all were obtained through the platform data library, route recorded in
`data/raw/<name>/.route`. Only the CDC dataset-metadata JSON (licence/attribution record, 23 kB) was fetched
directly with `curl`, because the library had no copy of the metadata endpoint.

Downloaded files are treated as untrusted: each source lives in its own directory under `data/raw/`, code
lives in `src/`, paths are passed as arguments, and readers run under `python3 -I`.

## 1. Sources

| dir (`data/raw/`) | content | publisher / URL | licence | route |
|---|---|---|---|---|
| `jhu_confirmed/time_series_covid19_confirmed_US.csv` | cumulative confirmed cases, county × day | JHU CSSE COVID-19 repository, `github.com/CSSEGISandData/COVID-19` (`csse_covid_19_time_series/`) | CC BY 4.0 | platform data library |
| `jhu_deaths/time_series_covid19_deaths_US.csv` | used only for its `Population` column | same | CC BY 4.0 | platform data library |
| `cdc_masks/rows.csv` | "U.S. State and Territorial Public Mask Mandates From April 10, 2020 through August 15, 2021 by County by Day", dataset `62d6-pm5i` | CDC, `data.cdc.gov/api/views/62d6-pm5i/rows.csv` | US federal work (public domain); metadata licence field empty; attribution: Mara Howard-Williams, Public Health Law Program, CDC | platform data library |
| `cdc_masks/metadata_62d6-pm5i.json` | dataset metadata (attribution, description) | `data.cdc.gov/api/views/62d6-pm5i.json` | as above | direct curl (not in library) |
| `census_adj2010/county_adjacency.txt` | 2010 county adjacency file | US Census Bureau, `www2.census.gov/geo/docs/reference/county_adjacency.txt` | public domain | platform data library |
| `acs_commuting/table1.xlsx` | ACS 2011–2015 county-to-county commuting flows | US Census Bureau, `www2.census.gov/programs-surveys/demo/tables/metro-micro/2015/commuting-flows-2015/table1.xlsx` | public domain | platform data library |

Upstream status: JHU stopped updating on 2023-03-10 and the CDC dataset is final, so the files are frozen;
a checksum mismatch means a corrupted or substituted file and `fetch_data.sh` aborts.

## 2. Integrity (recorded 2026-11-12; checked by `fetch_data.sh` against `data/checksums.tsv`)

| file | bytes | sha256 |
|---|---:|---|
| jhu_confirmed/time_series_covid19_confirmed_US.csv | 17,801,222 | e8fe5837b1a7ce09e6187ef778bf69fdcfa9fcd7015b1ba597f248b0736d8cf3 |
| jhu_deaths/time_series_covid19_deaths_US.csv | 11,865,683 | 5eab93b48daea1c31051cbcb0701b73d66105babb5f0e98f44c89173d86a37d6 |
| census_adj2010/county_adjacency.txt | 726,724 | b6df774db9ad5b94c47036379e207dd4622e46ef22c01ca2b7603453a4901014 |
| cdc_masks/rows.csv | 217,183,558 | 2354fc3169c01a5fff9297ae3a232fc2c441d4c66e191433b6dc6d86d22d70c0 |
| cdc_masks/metadata_62d6-pm5i.json | 23,334 | d7d3ad8e24714579c7fd275051941443037e843c27e23ed4f43de3fb8dc3cf5c |
| acs_commuting/table1.xlsx | 7,046,877 | f98b113197b8426fdf707bbffc219a27ae55842222bc67932210c9fe8fc47610 |

## 3. Key raw variables

* JHU: `UID`, `FIPS`, `Admin2`, `Province_State`, `Combined_Key`, daily cumulative columns `M/D/YY`;
  deaths file: `UID`, `Population`.
* CDC: `FIPS_State`, `FIPS_County`, `date` (`M/D/YYYY`), `Face_Masks_Required_in_Public`
  (`Yes`/`No`/`NA`), `order_code`; every county-day of the window is present (checked).
* Census adjacency: tab-separated, (county name, FIPS, neighbour name, neighbour FIPS), the focal county
  given on the first line of each block only; read with latin-1.
* ACS commuting: 7 header rows; residence state/county, workplace state (3-digit) / county, workers.

## 4. Processing rules (`src/data_prep.py`, fixed before any estimation)

1. Contiguous US + DC (Alaska, Hawaii, territories excluded); clusters = 49 "states" (48 + DC).
2. Window: estimation weeks 2020-06-01 … 2020-11-23 (week starts; T = 26 Monday-start weeks) plus
   S = 2 pre-sample lag weeks from 2020-05-18 (28 weeks stored).
3. Weekly cases = difference of cumulative counts at week ends; negative revisions truncated at 0
   (417 unit-weeks, total truncated mass 4,006 cases).
4. Reporting units: Utah health districts reported by JHU are merged with their member counties (6 units);
   Dukes + Nantucket (MA) merged; Kansas City (MO) assigned to Jackson County; Michigan prison rows and
   "Out of …"/"Unassigned" rows dropped. Cases in the window: 11,751,953 over all JHU rows, 11,528,182 in
   the analysis units (8,711 in dropped contiguous-US rows; the rest are AK/HI/territories).
5. Treatment: D_it = 1 if a state/territorial public mask order is in force on ≥ 4 of the 7 days of the
   week (population-weighted day share for merged units); `Dshare` keeps the continuous share. Opt-outs
   permitted under state law (notably Kansas) are **not** coded in the CDC data; see §7.
6. Networks: W_adj = row-normalised symmetric 2010 Census adjacency (two recodes: 46113→46102 Shannon→
   Oglala Lakota, 51515→51019 Bedford city→Bedford County; Connecticut uses its 8 historical counties in
   both JHU and the 2010 file); W_comm (sensitivity) = row-normalised C + Cᵀ of ACS commuter counts.
7. Exposure offset ω_i = pop_i / 10⁵ (importation term; ×0.25 and ×4 sensitivities).

## 5. Processed file `data/processed/panel.npz`

| key | shape / type | meaning |
|---|---|---|
| `Y` | (28, 3091) int | weekly new cases |
| `D` | (28, 3091) {0,1} | mask order indicator (rule 5) |
| `Dshare` | (28, 3091) [0,1] | share of days under an order |
| `pop` | (3091,) | JHU population |
| `state` | (3091,) | state FIPS (cluster) |
| `units` | (3091,) str | unit id (county FIPS or merged id) |
| `edges` | (9185, 2) | adjacency pairs i < j |
| `comm` | (131124, 3) | commuting (i, j, workers) |
| `week_start` | (28,) | Monday dates |

Verified by `data_prep.py` (`panel_summary.json`) and again at load time by `application.verify`
(and `tests/test_application_data.py`): 3,091 units, 49 clusters, 28 weeks, Y non-negative integer, D binary,
pop > 0, no isolated unit, no adjacency FIPS outside the units, no CDC FIPS missing. Descriptives:
56.8 % of unit-weeks treated; 2,291 units ever treated in the window, 571 treated at its start; 82 units
switch on→off (D is not strictly absorbing — handled by the estimator, which does not require absorbing
treatment); median weekly cases 27; 7.4 % zero unit-weeks.

## 6. Generation interval (weekly renewal weights)

Daily generation-interval distributions are discretised to weeks by
p_s = E_u[F(7(s+1) − u) − F(7s − u)], u ~ U[0, 7) (infection day uniform within the week). Lags s = 1, 2 are
kept and renormalised to w; the same-week mass (s = 0, not representable in a weekly lagged renewal model)
and the mass beyond S = 2 are reported in `application_results.json` so that the approximation is visible.
* Main: Ferretti et al. (2020, Science) Weibull(shape 2.826, scale 5.665) → w ≈ (0.964, 0.036), same-week
  mass ≈ 0.30.
* Sensitivity: Nishiura et al. (2020, IJID) serial interval lognormal(mean 4.7, sd 2.9) → w ≈ (0.913,
  0.087), same-week mass ≈ 0.39; plus a grid w₁ ∈ {0.6, …, 1.0}.

## 7. Known limitations recorded before analysis

* Kansas: county opt-outs from the state order are not in the CDC data, so D overstates treatment there;
  the pre-specified sensitivity drops Kansas and Texas (whose order also had county exemptions). A
  Kansas-only "opt-out" sub-analysis is **not** run because the opt-out status is not available in the
  data obtained.
* Confirmed cases depend on testing; testing volume is not modelled (absorbed only to the extent captured by
  unit and week effects).
* The 2010 adjacency and 2011–2015 commuting flows are pre-pandemic structural networks.
