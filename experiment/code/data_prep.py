"""Build the county-week analysis panel for the application (method_brief.md section 9).

Inputs (data/raw/, untrusted; run this script with `python3 -I`):
  jhu_confirmed/time_series_covid19_confirmed_US.csv   cumulative confirmed cases (JHU CSSE, CC BY 4.0)
  jhu_deaths/time_series_covid19_deaths_US.csv         Population column
  cdc_masks/rows.csv                                   state/territorial mask orders by county-day (CDC 62d6-pm5i)
  census_adj2010/county_adjacency.txt                  2010 county adjacency (Census)
  acs_commuting/table1.xlsx  (optional)                ACS 2011-2015 county-to-county commuting flows
Output: data/processed/panel.npz and data/processed/panel_summary.json.

Rules (pre-specified):
  * contiguous US + DC (Alaska, Hawaii and territories excluded); clusters = states;
  * window 2020-06-01 .. 2020-11-29 (T = 26 Monday-start weeks) plus S = 2 lag weeks from 2020-05-18;
  * weekly cases = difference of cumulative counts at week ends; negative revisions truncated to 0;
  * reporting units: JHU reports Utah by health district for some counties, so districts and their member
    counties are merged into one unit; Dukes and Nantucket (MA) merged; Kansas City (MO) assigned to Jackson
    County; prison rows (MI) and "Out of"/"Unassigned" rows dropped (their share of cases is reported);
  * D_it = 1 if a state order requiring masks in public applies on >= 4 of the 7 days (population-weighted
    over member counties for merged units). Opt-outs under state law are NOT coded in the CDC data.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.dirname(HERE)
RAW = os.path.join(RUN, "data", "raw")
OUT = os.path.join(RUN, "data", "processed")

START = pd.Timestamp("2020-06-01")
T_WEEKS = 26
S_LAGS = 2
EXCLUDE_STATE_FIPS = {2, 15}                      # Alaska, Hawaii (non-contiguous)
UTAH_DISTRICTS = {
    "Bear River": [49003, 49005, 49033],
    "Central Utah": [49023, 49027, 49031, 49039, 49041, 49055],
    "Southeast Utah": [49007, 49015, 49019],
    "Southwest Utah": [49001, 49017, 49021, 49025, 49053],
    "TriCounty": [49009, 49013, 49047],
    "Weber-Morgan": [49057, 49029],
}
ADJ_RECODE = {46113: 46102, 51515: 51019}         # Shannon->Oglala Lakota (2015), Bedford city->Bedford Co (2013)


def week_starts():
    first = START - pd.Timedelta(weeks=S_LAGS)
    return [first + pd.Timedelta(weeks=k) for k in range(S_LAGS + T_WEEKS)]


def load_jhu():
    c = pd.read_csv(os.path.join(RAW, "jhu_confirmed", "time_series_covid19_confirmed_US.csv"))
    d = pd.read_csv(os.path.join(RAW, "jhu_deaths", "time_series_covid19_deaths_US.csv"),
                    usecols=["UID", "Population"])
    c = c.merge(d, on="UID", how="left", validate="one_to_one")
    return c


def county_unit_map(c: pd.DataFrame):
    """Map JHU rows to reporting units; return (row -> unit key) series and unit membership dict."""
    member = {}                                    # unit key -> list of county FIPS
    row_unit = pd.Series(index=c.index, dtype=object)
    district_of = {f: nm for nm, fl in UTAH_DISTRICTS.items() for f in fl}
    dropped = []
    for i, r in c.iterrows():
        f = r["FIPS"]
        key = None
        if pd.isna(f):
            ck = str(r["Combined_Key"])
            if ck.startswith("Dukes and Nantucket"):
                key = "U25007"
            elif ck.startswith("Kansas City"):
                key = "C29095"
            else:
                nm = ck.split(",")[0].strip()
                if nm in UTAH_DISTRICTS:
                    key = "UT_" + nm
        else:
            f = int(f)
            st = f // 1000
            if f < 1000 or st >= 57 or st in EXCLUDE_STATE_FIPS or f % 1000 == 0 or f >= 80000:
                key = None
            elif f in district_of:
                key = "UT_" + district_of[f]
            elif f in (25007, 25019):
                key = "U25007"
            else:
                key = "C%05d" % f
        row_unit[i] = key
        if key is None:
            dropped.append(str(r["Combined_Key"]))
        elif not pd.isna(r["FIPS"]):
            member.setdefault(key, []).append(int(r["FIPS"]))
        else:
            member.setdefault(key, [])
    return row_unit, member, dropped


def weekly_cases(c: pd.DataFrame, row_unit: pd.Series):
    ws = week_starts()
    date_cols = {pd.Timestamp(pd.to_datetime(col, format="%m/%d/%y")): col for col in c.columns
                 if col.count("/") == 2}
    ends = [w + pd.Timedelta(days=6) for w in ws]
    prev = [w - pd.Timedelta(days=1) for w in ws]
    keep = row_unit.notna()
    cum_end = c.loc[keep, [date_cols[e] for e in ends]].to_numpy(float)
    cum_prev = c.loc[keep, [date_cols[p] for p in prev]].to_numpy(float)
    units = row_unit[keep].to_numpy()
    inc = cum_end - cum_prev                                   # rows x weeks
    df = pd.DataFrame(inc, columns=range(len(ws)))
    df["unit"] = units
    agg = df.groupby("unit").sum()
    neg = int((agg.to_numpy() < 0).sum())
    neg_mass = float(-agg.to_numpy()[agg.to_numpy() < 0].sum())
    Y = np.clip(agg.to_numpy(), 0, None).T                     # weeks x units
    # share of window cases in dropped rows (Out of / Unassigned / prisons / AK / HI excluded separately)
    all_inc = (c[[date_cols[e] for e in ends]].to_numpy(float) - c[[date_cols[p] for p in prev]].to_numpy(float))
    return agg.index.to_list(), Y, neg, neg_mass, all_inc


def mask_policy(units, member, pop_county):
    cdc = pd.read_csv(os.path.join(RAW, "cdc_masks", "rows.csv"),
                      usecols=["FIPS_State", "FIPS_County", "date", "Face_Masks_Required_in_Public"],
                      dtype={"Face_Masks_Required_in_Public": "object"})
    cdc["date"] = pd.to_datetime(cdc["date"], format="%m/%d/%Y")
    cdc["fips"] = cdc["FIPS_State"] * 1000 + cdc["FIPS_County"]
    ws = week_starts()
    lo, hi = ws[0], ws[-1] + pd.Timedelta(days=6)
    cdc = cdc[(cdc.date >= lo) & (cdc.date <= hi)].copy()
    cdc["yes"] = (cdc["Face_Masks_Required_in_Public"] == "Yes").astype(float)
    cdc["week"] = ((cdc.date - lo).dt.days // 7).astype(int)
    share = cdc.groupby(["fips", "week"]).yes.mean().unstack("week")      # county x week share of days
    ndays = cdc.groupby(["fips", "week"]).size().unstack("week")
    missing, Dshare = [], np.zeros((len(ws), len(units)))
    for j, u in enumerate(units):
        fl = member[u]
        ok = [f for f in fl if f in share.index]
        if len(ok) < len(fl):
            missing += [f for f in fl if f not in share.index]
        if not ok:
            Dshare[:, j] = np.nan
            continue
        wts = np.array([pop_county.get(f, 0.0) for f in ok], float)
        wts = wts / wts.sum() if wts.sum() > 0 else np.full(len(ok), 1.0 / len(ok))
        Dshare[:, j] = (share.loc[ok].to_numpy() * wts[:, None]).sum(0)
    full_weeks = bool((ndays.to_numpy() == 7).all())
    return Dshare, missing, full_weeks


def adjacency(units, member):
    unit_of = {f: u for u, fl in member.items() for f in fl}
    path = os.path.join(RAW, "census_adj2010", "county_adjacency.txt")
    edges = set()
    cur = None
    unknown = set()
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            if parts[1].strip():
                cur = ADJ_RECODE.get(int(parts[1]), int(parts[1]))
            nb = ADJ_RECODE.get(int(parts[3]), int(parts[3]))
            a, b = unit_of.get(cur), unit_of.get(nb)
            if a is None or b is None:
                if cur // 1000 not in EXCLUDE_STATE_FIPS and cur // 1000 < 57:
                    unknown.add(cur if a is None else nb)
                continue
            if a != b:
                edges.add((min(a, b), max(a, b)))
    idx = {u: k for k, u in enumerate(units)}
    E = np.array(sorted((idx[a], idx[b]) for a, b in edges), int)
    return E, sorted(f for f in unknown if f // 1000 < 57 and f // 1000 not in EXCLUDE_STATE_FIPS)


def commuting_edges(units, member):
    """ACS 2011-2015 county-to-county commuting flows (optional). Returns (i, j, workers) for i != j."""
    path = os.path.join(RAW, "acs_commuting", "table1.xlsx")
    if not os.path.exists(path):
        return None
    # layout (checked): 7 title/header rows; columns = residence state, county, names, workplace state
    # (3-digit, "0SS" for US states), workplace county, names, workers, margin of error
    df = pd.read_excel(path, header=None, skiprows=7, dtype=str, engine="openpyxl")
    df = df.iloc[:, [0, 1, 4, 5, 8]]
    df.columns = ["rs", "rc", "ws", "wc", "flow"]
    df = df.dropna()
    df = df[df.rs.str.fullmatch(r"[0-9]+") & df.ws.str.fullmatch(r"0[0-9][0-9]")
            & df.wc.str.fullmatch(r"[0-9][0-9][0-9]")]
    df["r"] = df.rs.str[-2:].astype(int) * 1000 + df.rc.astype(int)
    df["w"] = df.ws.str[-2:].astype(int) * 1000 + df.wc.astype(int)
    df["flow"] = pd.to_numeric(df.flow.str.replace(",", ""), errors="coerce")
    unit_of = {f: u for u, fl in member.items() for f in fl}
    unit_of.update({ADJ_RECODE.get(f, f): unit_of[f] for f in list(unit_of)})
    idx = {u: k for k, u in enumerate(units)}
    df["i"] = df.r.map(lambda f: idx.get(unit_of.get(ADJ_RECODE.get(f, f))))
    df["j"] = df.w.map(lambda f: idx.get(unit_of.get(ADJ_RECODE.get(f, f))))
    df = df.dropna(subset=["i", "j", "flow"])
    df = df[df.i != df.j]
    g = df.groupby(["i", "j"]).flow.sum().reset_index()
    return g[["i", "j", "flow"]].to_numpy(float)


def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    c = load_jhu()
    row_unit, member, dropped = county_unit_map(c)
    units, Yu, neg, neg_mass, all_inc = weekly_cases(c, row_unit)
    pop_county = {int(f): float(p) for f, p in zip(c.FIPS, c.Population) if not pd.isna(f)}
    pop = np.array([sum(pop_county.get(f, 0.0) for f in member[u]) for u in units])
    state = np.array([int(member[u][0]) // 1000 if member[u] else 29 for u in units])
    Dshare, missing_cdc, full_weeks = mask_policy(units, member, pop_county)
    # drop units without policy coding or with zero population
    keep = np.isfinite(Dshare).all(0) & (pop > 0)
    units = [u for u, k in zip(units, keep) if k]
    Yu, Dshare, pop, state = Yu[:, keep], Dshare[:, keep], pop[keep], state[keep]
    E, unknown_adj = adjacency(units, member)
    comm = commuting_edges(units, member)
    D = (Dshare >= 4.0 / 7.0 - 1e-12).astype(float)
    ws = week_starts()
    # summary of dropped mass
    keep_rows = row_unit.notna()
    win = slice(S_LAGS, S_LAGS + T_WEEKS)
    total_us = float(np.clip(all_inc[:, win], 0, None).sum())
    st_all = (c.FIPS.fillna(0) // 1000).astype(int)
    contig = ~st_all.isin(list(EXCLUDE_STATE_FIPS)) & (st_all < 57) & ~c.Province_State.isin(
        ["Alaska", "Hawaii", "American Samoa", "Guam", "Northern Mariana Islands", "Puerto Rico", "Virgin Islands",
         "Diamond Princess", "Grand Princess"])
    dropped_contig = float(np.clip(all_inc[(~keep_rows & contig).to_numpy()][:, win], 0, None).sum())
    np.savez_compressed(os.path.join(OUT, "panel.npz"), Y=Yu, D=D, Dshare=Dshare, pop=pop, state=state,
                        units=np.array(units), edges=E,
                        comm=comm if comm is not None else np.zeros((0, 3)),
                        week_start=np.array([str(w.date()) for w in ws]))
    summ = dict(
        created_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), runtime_sec=round(time.time() - t0, 1),
        n_units=len(units), n_states=int(len(np.unique(state))), weeks_total=len(ws), S=S_LAGS, T=T_WEEKS,
        first_week=str(ws[0].date()), first_estimation_week=str(ws[S_LAGS].date()), last_week=str(ws[-1].date()),
        n_edges=int(len(E)), isolated_units=[units[k] for k in range(len(units))
                                             if k not in set(E.ravel().tolist())],
        adjacency_fips_not_in_units=unknown_adj[:50],
        commuting_edges=None if comm is None else int(len(comm)),
        negative_weekly_increments_truncated=neg, negative_mass_truncated=neg_mass,
        cases_window_all_rows=total_us, cases_window_units=float(Yu[win].sum()),
        cases_window_dropped_rows_contiguous=dropped_contig,
        dropped_rows_examples=sorted(set(dropped))[:30], cdc_fips_missing=missing_cdc,
        cdc_all_unit_weeks_have_7_days=full_weeks,
        share_unit_weeks_treated=float(D[win].mean()), median_weekly_cases=float(np.median(Yu[win])),
        frac_zero_unit_weeks=float((Yu[win] == 0).mean()),
        n_units_ever_treated_in_window=int((D[win].max(0) == 1).sum()),
        n_units_treated_at_window_start=int((D[S_LAGS] == 1).sum()),
        n_units_on_to_off=int(((D[win][1:] == 0) & (D[win][:-1] == 1)).any(0).sum()),
        merged_units={u: member[u] for u in units if not u.startswith("C") or u == "C29095"},
    )
    with open(os.path.join(OUT, "panel_summary.json"), "w") as fh:
        json.dump(summ, fh, indent=1)
    print(json.dumps({k: v for k, v in summ.items() if k not in ("merged_units",)}, indent=1))


if __name__ == "__main__":
    sys.exit(main())
