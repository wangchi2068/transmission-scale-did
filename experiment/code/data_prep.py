"""Build the processed German county-level weekly Varicella panel (2015-2019 + 3 pre-sample weeks).

Inputs (experiment/data/raw/germany/):
  - rki_varicella_2014_2019.json           Weekly RKI SurvStat Varicella counts + 2017 populations
  - germany_school_holidays_2014_2019.json Official 16-state KMK school holiday dates (OpenHolidays)
  - counties_ags.geojson                   5-digit AGS German county polygons

Outputs (experiment/data/processed/):
  - panel.npz
  - panel_summary.json
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import re
from collections import defaultdict, deque
from pathlib import Path

import numpy as np
import scipy.sparse as sps

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "experiment" / "data"
RAW_DIR = DATA_DIR / "raw" / "germany"
OUT_DIR = DATA_DIR / "processed"

EXPECTED_SHA256 = {
    "rki_varicella_2014_2019.json": "0735f0211c4d51389a20cb93dbcd4d5b48f8f7790d1c68176c8841c12d0b2489",
    "germany_school_holidays_2014_2019.json": "4f51fd93e04f0f7ae37a5d78aea0966d12f74acf15b7293a795f155826fb695e",
    "counties_ags.geojson": "5003586e2c43f803828b106bb7de11b58079529cb78bf76649e3d2ec703acacc",
}

STATE_ABBR = {
    1: "SH", 2: "HH", 3: "NI", 4: "HB", 5: "NW", 6: "HE", 7: "RP", 8: "BW",
    9: "BY", 10: "SL", 11: "BE", 12: "BB", 13: "MV", 14: "SN", 15: "ST", 16: "TH",
}

# 4 German macro-regions: 0=North, 1=West, 2=South, 3=East
REGION_OF_STATE = {
    1: 0, 2: 0, 3: 0, 4: 0, 13: 0,
    5: 1, 6: 1, 7: 1, 10: 1,
    8: 2, 9: 2,
    11: 3, 12: 3, 14: 3, 15: 3, 16: 3,
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_checksums(raw_dir: Path = RAW_DIR) -> dict:
    verified = {}
    mismatches = []
    for fname, expected in EXPECTED_SHA256.items():
        p = raw_dir / fname
        if not p.exists():
            mismatches.append(f"{fname}: missing (expected {expected})")
            continue
        actual = sha256_file(p)
        if actual != expected:
            mismatches.append(f"{fname}: got {actual}, expected {expected}")
        else:
            verified[fname] = dict(bytes=p.stat().st_size, sha256=actual)
    if mismatches:
        raise RuntimeError("Checksum verification failed:\n  " + "\n  ".join(mismatches))
    return verified


def extract_ags(rki_key: str):
    m = re.findall(r"&\[([^\]]+)\]", rki_key)
    if not m:
        return None, None
    state = m[0]
    ags = m[-1]
    if state == "11":  # Merge Berlin's 12 boroughs into city-state 11000
        ags = "11000"
    return state, ags


def nuts2_code(ags: str) -> str:
    """Map 5-digit AGS to Germany's 38 NUTS-2 Regierungsbezirke (or state if single NUTS-2)."""
    st = ags[:2]
    # States with sub-state NUTS-2 regions (3-digit AGS prefix defines Regierungsbezirk):
    # NI (03: 031..034), NW (05: 051,053,055,057,059), HE (06: 064..066),
    # RP (07: 071..073), BW (08: 081..084), BY (09: 091..097), SN (14: 145..147)
    if st in {"03", "05", "06", "07", "08", "09", "14"}:
        return ags[:3]
    return st + "0"


def build_adjacency(geo_path: Path, ags_list: list[str]):
    gj = json.loads(geo_path.read_text(encoding="utf-8"))
    ags_to_idx = {a: i for i, a in enumerate(ags_list)}
    v2polys = defaultdict(set)
    matched = set()

    def iter_rings(geom):
        t = geom["type"]
        coords = geom["coordinates"]
        if t == "Polygon":
            for ring in coords:
                yield ring
        elif t == "MultiPolygon":
            for poly in coords:
                for ring in poly:
                    yield ring

    for feat in gj["features"]:
        raw = str(feat.get("properties", {}).get("ags") or feat.get("id") or "").strip()
        if raw.isdigit():
            raw = raw.zfill(5)[:5]
        # Known AGS consolidations
        if raw in ("03152", "03156") and "03159" in ags_to_idx:
            raw = "03159"
        if raw == "16056" and "16063" in ags_to_idx:
            raw = "16063"
        if raw.startswith("11") and "11000" in ags_to_idx:
            raw = "11000"
        if raw not in ags_to_idx:
            continue
        idx = ags_to_idx[raw]
        matched.add(idx)
        for ring in iter_rings(feat["geometry"]):
            for pt in ring:
                v2polys[(round(pt[0], 3), round(pt[1], 3))].add(idx)

    if len(matched) != len(ags_list):
        raise RuntimeError(f"Expected {len(ags_list)} matched counties, got {len(matched)}")

    edges_set = set()
    for polys in v2polys.values():
        if len(polys) >= 2:
            pl = sorted(polys)
            for a_i in range(len(pl)):
                for b_i in range(a_i + 1, len(pl)):
                    edges_set.add((pl[a_i], pl[b_i]))

    edges = np.array(sorted(edges_set), dtype=np.int32)
    N = len(ags_list)
    rows = np.concatenate([edges[:, 0], edges[:, 1]])
    cols = np.concatenate([edges[:, 1], edges[:, 0]])
    A = sps.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(N, N)).tocsr()
    A.setdiag(0.0)
    A.eliminate_zeros()
    deg = np.asarray(A.sum(1)).ravel()
    if (deg == 0).any():
        raise RuntimeError(f"Found {(deg == 0).sum()} isolated counties")
    W = sps.diags(1.0 / deg) @ A
    return W.tocsr(), edges, deg


def graph_diameter(adj: np.ndarray) -> int:
    n = adj.shape[0]
    diam = 0
    for s in range(n):
        dist = [-1] * n
        dist[s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for v in np.where(adj[u])[0]:
                if dist[v] < 0:
                    dist[v] = dist[u] + 1
                    q.append(v)
        if min(dist) < 0:
            return -1
        diam = max(diam, max(dist))
    return diam


def build_holiday_matrices(hol_json: dict, weeks: list[tuple[int, int]], ags_list: list[str],
                           ags_state: dict[str, str], holiday_filter: str | None, min_days: int = 4):
    hol_days = defaultdict(set)
    for st_code, srec in hol_json.items():
        for h in srec["holidays"]:
            hname = " ".join(n["text"] for n in h.get("name", [])).lower()
            if holiday_filter and holiday_filter not in hname:
                continue
            d0 = dt.date.fromisoformat(h["startDate"])
            d1 = dt.date.fromisoformat(h["endDate"])
            cur = d0
            while cur <= d1:
                hol_days[st_code].add(cur)
                cur += dt.timedelta(days=1)

    T_all = len(weeks)
    N = len(ags_list)
    D = np.zeros((T_all, N), dtype=np.float64)
    Dshare = np.zeros((T_all, N), dtype=np.float64)
    for t, (yr, wk) in enumerate(weeks):
        try:
            mon = dt.date.fromisocalendar(yr, wk, 1)
        except ValueError:
            mon = dt.date.fromisocalendar(yr, wk - 1, 1) + dt.timedelta(days=7)
        days = [mon + dt.timedelta(days=k) for k in range(7)]
        for i, a in enumerate(ags_list):
            st = ags_state[a]
            cnt = sum(1 for dd in days if dd in hol_days[st])
            Dshare[t, i] = cnt / 7.0
            D[t, i] = 1.0 if cnt >= min_days else 0.0
    return D, Dshare


def main():
    checksums = verify_checksums()
    rki = json.loads((RAW_DIR / "rki_varicella_2014_2019.json").read_text(encoding="utf-8"))
    hol = json.loads((RAW_DIR / "germany_school_holidays_2014_2019.json").read_text(encoding="utf-8"))

    years = (2014, 2015, 2016, 2017, 2018, 2019)
    ags_state = {}
    for k in rki["years"][str(years[0])]["rows"]:
        st, ags = extract_ags(k)
        if ags and ags.isdigit() and len(ags) == 5:
            ags_state[ags] = st
    ags_list = sorted(ags_state.keys())
    ags_idx = {a: i for i, a in enumerate(ags_list)}
    N = len(ags_list)

    weeks = []
    Y_rows = []
    for yr in years:
        yd = rki["years"][str(yr)]
        n_w = len(yd["cols"])
        Y_yr = np.zeros((n_w, N), dtype=np.float64)
        for k, rec in yd["rows"].items():
            _, ags = extract_ags(k)
            if ags in ags_idx:
                Y_yr[:, ags_idx[ags]] += np.asarray(rec["counts"], dtype=np.float64)
        for _, w_cap in yd["cols"]:
            weeks.append((yr, int(w_cap)))
        Y_rows.append(Y_yr)
    Y_all = np.vstack(Y_rows)

    pop = np.zeros(N, dtype=np.float64)
    for k, rec in rki["population_2017"].items():
        _, ags = extract_ags(k)
        if ags in ags_idx:
            pop[ags_idx[ags]] += float(rec["pop"])
    for k, rec in rki.get("population_fallback", {}).items():
        _, ags = extract_ags(k)
        if ags in ags_idx and pop[ags_idx[ags]] <= 0:
            pop[ags_idx[ags]] = float(rec["pop"])
    n_imputed_pop = int((pop <= 0).sum())
    med_pop = float(np.median(pop[pop > 0]))
    pop = np.where(pop > 0, pop, med_pop)

    D_sum_all, Dshare_sum_all = build_holiday_matrices(hol, weeks, ags_list, ags_state, "summer", 4)
    D_sum_3d_all, _ = build_holiday_matrices(hol, weeks, ags_list, ags_state, "summer", 3)
    D_sum_5d_all, _ = build_holiday_matrices(hol, weeks, ags_list, ags_state, "summer", 5)
    D_any_all, Dshare_any_all = build_holiday_matrices(hol, weeks, ags_list, ags_state, None, 4)
    D_east_all, _ = build_holiday_matrices(hol, weeks, ags_list, ags_state, "easter", 4)
    D_aut_all, _ = build_holiday_matrices(hol, weeks, ags_list, ags_state, "autumn", 4)
    D_win_all, _ = build_holiday_matrices(hol, weeks, ags_list, ags_state, "winter", 4)
    D_pent_all, _ = build_holiday_matrices(hol, weeks, ags_list, ags_state, "pentecost", 4)

    W, edges, deg = build_adjacency(RAW_DIR / "counties_ags.geojson", ags_list)

    state_arr = np.array([int(ags_state[a]) for a in ags_list], dtype=np.int32)
    states = np.unique(state_arr)
    cl = np.searchsorted(states, state_arr).astype(np.int32)
    region = np.array([REGION_OF_STATE[int(s)] for s in state_arr], dtype=np.int32)
    nuts2_labels = [nuts2_code(a) for a in ags_list]
    nuts2_uniq = sorted(set(nuts2_labels))
    nuts2 = np.array([nuts2_uniq.index(lbl) for lbl in nuts2_labels], dtype=np.int32)

    # State contiguity graph & diameter
    G = len(states)
    A_st = np.zeros((G, G), dtype=bool)
    for i, j in edges:
        if cl[i] != cl[j]:
            A_st[cl[i], cl[j]] = True
            A_st[cl[j], cl[i]] = True
    st_edges = int(A_st.sum() // 2)
    st_diam = graph_diameter(A_st)

    # Slice to 2014-W50 .. 2019-W52 (3 pre-sample weeks + 261 estimation weeks = 264 rows)
    idx_2015_w1 = weeks.index((2015, 1))
    start_idx = idx_2015_w1 - 3  # 2014-W50
    Y = Y_all[start_idx:]
    D_summer = D_sum_all[start_idx:]
    D_summer_3d = D_sum_3d_all[start_idx:]
    D_summer_5d = D_sum_5d_all[start_idx:]
    Dshare_summer = Dshare_sum_all[start_idx:]
    D_all = D_any_all[start_idx:]
    Dshare_all = Dshare_any_all[start_idx:]
    D_easter = D_east_all[start_idx:]
    D_autumn = D_aut_all[start_idx:]
    D_winter = D_win_all[start_idx:]
    D_pentecost = D_pent_all[start_idx:]
    weeks_win = weeks[start_idx:]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        OUT_DIR / "panel.npz",
        Y=Y,
        D_summer=D_summer,
        D_summer_3d=D_summer_3d,
        D_summer_5d=D_summer_5d,
        Dshare_summer=Dshare_summer,
        D_all=D_all,
        Dshare_all=Dshare_all,
        D_easter=D_easter,
        D_autumn=D_autumn,
        D_winter=D_winter,
        D_pentecost=D_pentecost,
        Y_full_2014_2019=Y_all,
        D_summer_full=D_sum_all,
        D_all_full=D_any_all,
        pop=pop,
        state=state_arr,
        cl=cl,
        states=states,
        region=region,
        nuts2=nuts2,
        nuts2_codes=np.array(nuts2_uniq),
        ags=np.array(ags_list),
        weeks=np.array(weeks_win, dtype=np.int32),
        edges=edges,
        W_indptr=W.indptr,
        W_indices=W.indices,
        W_data=W.data,
        W_shape=np.array(W.shape, dtype=np.int32),
    )

    Y_est = Y[3:]
    D_sum_est_lag1 = np.vstack([np.zeros((1, N)), D_summer[:-1]])[3:]
    counts_per_state = [int((cl == g).sum()) for g in range(G)]
    # Summer onset weeks across 16 states for each year 2015..2019
    stagger_by_year = {}
    for yr in range(2015, 2020):
        st_starts = {}
        for g, st_num in enumerate(states):
            col = np.where(cl == g)[0][0]
            wks_on = [wk for t, (y, wk) in enumerate(weeks_win) if y == yr and D_summer[t, col] == 1.0]
            if wks_on:
                st_starts[STATE_ABBR[int(st_num)]] = dict(start_week=min(wks_on), end_week=max(wks_on), n_weeks=len(wks_on))
        starts_vals = [v["start_week"] for v in st_starts.values()]
        stagger_by_year[str(yr)] = dict(
            earliest_start_week=min(starts_vals),
            latest_start_week=max(starts_vals),
            span_weeks=max(starts_vals) - min(starts_vals),
            states=st_starts,
        )

    summary = dict(
        N=N,
        G=G,
        T_est=int(Y_est.shape[0]),
        S_pre=3,
        first_presample_week=list(weeks_win[0]),
        first_est_week=list(weeks_win[3]),
        last_est_week=list(weeks_win[-1]),
        total_cases_2015_2019=int(Y_est.sum()),
        total_cases_2014_2019=int(Y_all.sum()),
        mean_count=float(Y_est.mean()),
        median_count=float(np.median(Y_est)),
        iqr_count=[float(np.quantile(Y_est, 0.25)), float(np.quantile(Y_est, 0.75))],
        max_count=int(Y_est.max()),
        zero_share=float((Y_est == 0).mean()),
        summer_share_lag1=float(D_sum_est_lag1.mean()),
        all_holiday_share_lag0=float(D_all[3:].mean()),
        counties_per_state=dict(min=min(counts_per_state), median=float(np.median(counts_per_state)),
                                max=max(counts_per_state), by_state={STATE_ABBR[int(s)]: c for s, c in zip(states, counts_per_state)}),
        n_nuts2=len(nuts2_uniq),
        county_edges=int(len(edges)),
        county_mean_degree=float(deg.mean()),
        county_min_degree=int(deg.min()),
        county_max_degree=int(deg.max()),
        state_edges=st_edges,
        state_mean_degree=float(A_st.sum(1).mean()),
        state_graph_diameter=st_diam,
        total_population_2017=float(pop.sum()),
        median_county_pop=med_pop,
        n_imputed_pop_counties=n_imputed_pop,
        stagger_by_year=stagger_by_year,
        checksums=checksums,
    )
    (OUT_DIR / "panel_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k != "stagger_by_year"}, indent=2))
    print("\nSummer holiday staggering span by year:")
    for yr, rec in stagger_by_year.items():
        print(f"  {yr}: weeks {rec['earliest_start_week']}..{rec['latest_start_week']} (span {rec['span_weeks']} wk)")


if __name__ == "__main__":
    main()
