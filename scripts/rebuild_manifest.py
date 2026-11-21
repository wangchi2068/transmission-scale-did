"""Rebuild tables/manifest.json registering ALL active LaTeX tables with LF SHA-256 digests
and authoritative generator / data-source provenance (resolves E8(a)).
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAB_DIR = ROOT / "tables"

TABLE_REGISTRY = [
    ("tab:sim-core-tau", "table_1_core_tau.tex", "experiment/code/make_tables.py", ["experiment/analysis/results.json"], 8000, "Main Table 1: core tau bias and coverage across comparators"),
    ("tab:sim-m0gi", "table_8_m0gi_posthoc.tex", "experiment/code/make_tables.py", ["experiment/analysis/results.json"], 8000, "Main Table 2: post-hoc M0-GI vs pre-specified M0 and M0-G"),
    ("tab:envelope", "table_envelope.tex", "experiment/code/make_tables_germany.py", ["experiment/analysis/results.json", "experiment/analysis/germany_calibrated_sim.json", "verification/same_week/results.json", "verification/root_selection/full_summary_tsummary.json"], None, "Main Table 3: operating envelope and diagnostic workflow"),
    ("tab:app", "table_7_application.tex", "experiment/code/make_tables_germany.py", ["experiment/analysis/germany_application_results.json"], 1, "Main Table 4: sparse German varicella illustrative implementation"),
    ("tab:sim-m0-full", "table_9_m0_full.tex", "experiment/code/make_tables.py", ["experiment/analysis/results.json"], 8000, "Web Table: full M0 tau performance measures with MCSEs"),
    ("tab:sim-core-spill", "table_2_core_spill.tex", "experiment/code/make_tables.py", ["experiment/analysis/results.json"], 8000, "Web Table: core spillover parameters kappa and phi"),
    ("tab:sim-null", "table_3_null.tex", "experiment/code/make_tables.py", ["experiment/analysis/results.json"], 8000, "Web Table: Type I error under null tau0=0"),
    ("tab:sim-roots", "table_6_roots.tex", "experiment/code/make_tables.py", ["experiment/analysis/results.json"], 8000, "Web Table: multiple roots and boundary solutions"),
    ("tab:sim-misspec", "table_4_misspec.tex", "experiment/code/make_tables.py", ["experiment/analysis/results.json"], 8000, "Web Table: M0 under X1-X5 misspecification"),
    ("tab:sim-psi", "table_5_psi.tex", "experiment/code/make_tables.py", ["experiment/analysis/results.json", "experiment/analysis/core_diagnostics.json"], 8000, "Web Table: cases averted psi inference"),
    ("tab:rev-rsel", "table_10_rsel.tex", "verification/revision/make_revision_tables.py", ["verification/root_selection/full_summary_tsummary.json"], 2000, "Web Table: boundary root selection validation"),
    ("tab:rev-poproots", "table_17_population_roots.tex", "verification/revision/make_revision_tables.py", ["verification/revision/population_roots/results.json"], None, "Web Table: large-N population second-root multiplicity"),
    ("tab:rev-nonnested", "table_19_nonnested_weakprior.tex", "experiment/code/make_verification_tables.py", ["verification/non_nested_sim/results.json", "verification/same_week/grid_results.json", "verification/weak_prior_m2/results.json"], None, "Web Table: non-nested daily DGP, same-week w0 grid, and weak-prior M2"),
    ("tab:rev-clustergrid", "table_20_cluster_grid.tex", "experiment/code/make_verification_tables.py", ["verification/root_selection/cluster_grid_results.json"], 250, "Web Table: controlled G=10..100 and skewness x T grid"),
    ("tab:app-all", "table_7b_application_sensitivity.tex", "experiment/code/make_tables_germany.py", ["experiment/analysis/germany_application_results.json"], 1, "Web Table: German application temporal, threshold, and subsample sensitivities"),
    ("tab:rev-gi", "table_13_gi_timing.tex", "experiment/code/make_tables_germany.py", ["experiment/analysis/germany_application_results.json"], 1, "Web Table: generation-interval and reporting-dispersion sensitivity"),
    ("tab:rev-regional", "table_11_regional.tex", "experiment/code/make_tables_germany.py", ["experiment/analysis/germany_application_results.json"], 1, "Web Table: regional FE, BW/BY exclusion, trends, and seasonal harmonics"),
    ("tab:rev-calsim", "table_12_calibrated_sim.tex", "experiment/code/make_tables_germany.py", ["experiment/analysis/germany_calibrated_sim.json"], 500, "Web Table: German design-calibrated few-cluster simulation"),
    ("tab:rev-robust", "table_16_robust_sets.tex", "experiment/code/make_tables_germany.py", ["experiment/analysis/germany_application_results.json"], 1, "Web Table: Stock-Wright confidence sets"),
    ("tab:rev-event", "table_14_event_weights.tex", "experiment/code/make_tables_germany.py", ["experiment/analysis/germany_application_results.json"], 1, "Web Table: event-time dynamic profile on total-policy scale"),
    ("tab:rev-ascert", "table_15_ascertainment.tex", "experiment/code/make_tables_germany.py", ["experiment/analysis/germany_application_results.json"], 1, "Web Table: ascertainment tipping points and omega sensitivity"),
    ("tab:rev-placebo", "table_18_placebos.tex", "experiment/code/make_tables_germany.py", ["experiment/analysis/germany_application_results.json"], 1, "Web Table: non-summer placebos, calendar shifts, and cross-state permutations"),
]


def main():
    res_meta = json.loads((ROOT / "experiment" / "analysis" / "results.json").read_text(encoding="utf-8"))
    entries = []
    for tid, fn, gen, srcs, n_seeds, note in TABLE_REGISTRY:
        p = TAB_DIR / fn
        if not p.exists():
            raise FileNotFoundError(f"Missing table file: {p}")
        dig = hashlib.sha256(p.read_bytes()).hexdigest()
        entries.append(dict(
            id=tid, path=fn, sha256=dig, generator=gen,
            data_sources=srcs, n_seeds=n_seeds, note=note,
        ))
    man = dict(
        schema_version=2,
        provenance=res_meta.get("provenance"),
        note="Complete provenance manifest covering all 22 active tables (Main Tables 1-4 and Web Tables 1-18).",
        tables=entries,
    )
    (TAB_DIR / "manifest.json").write_text(json.dumps(man, indent=2) + "\n", encoding="utf-8")
    print(f"Rebuilt tables/manifest.json with {len(entries)} verified entries.")


if __name__ == "__main__":
    main()
