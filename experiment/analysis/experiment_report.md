---
title: Transmission-scale difference-in-differences under feedback and contagion
subtitle: Monte Carlo study and US county COVID-19 application of the quasi-differenced GMM estimator
keywords: difference-in-differences; reproduction number; sequential exogeneity; spatial contagion; GMM; Monte Carlo
abstract: |
  We evaluate quasi-differenced GMM (QD-GMM) for policy effects on the reproduction number in a multiplicative renewal panel with policy feedback and network contagion. The study has 25 scenarios with 8000 replications each, comparison with standard estimators and an application to US county COVID-19 cases and mask orders (3091 counties, 49 states). When the coupling to neighbours' cases is positive, QD-GMM is close to unbiased for the direct effect (bias −0.002 to 0.001) with coverage of 94.5 % to 95.4 %, and its size under the null is 4.7 % to 5.2 %. Poisson fixed effects under-covers in all 16 core scenarios, and two-step log R and log(Y+1) regressions fail. When the coupling is zero, the estimating equations have spurious roots. QD-GMM then fails (coverage 0.6 % to 38.8 %), the pre-specified global root search only partly helps (52.9 % to 71.0 %), and a post-hoc selection rule restores coverage to 87.5 % to 95.2 %. Intervals for cases averted under-cover in some designs. In the application, the constant-effect estimate is 0.040 (95 % CI −0.217 to 0.296), the over-identification test rejects, and the event-time effects are all negative.
---

# Experiment Report: transmission-scale DiD under feedback and contagion

**Authored by:** Apevon Science (with human review strongly recommended before publication or production use)
**Mode:** Measured (Monte Carlo study run by `src/simulate.py` and `src/simulate_supp.py`; application on real data by `src/application.py`; all numbers aggregated by `src/aggregate.py` into `results.json`)
**Date:** 2026-11-12
**Hardware:** Linux pod with 16 CPU cores and 32 GiB RAM, no GPU; main Monte Carlo run 201.4 min with 16 worker processes (`logs/sim.log`); supplementary M0-GI run about 190 min at low priority on 4 workers in parallel with the main run, then 47.2 min on 16 workers for the remaining 9320 of 20000 tasks (`logs/run_supp_console.log`)
**Data contract:** `data_contract.md`

## Problem

Observational evaluations of non-pharmaceutical interventions such as mask orders usually apply two-way fixed-effects difference-in-differences to cases or log cases. That outcome scale is inconsistent with transmission dynamics. Three further problems arise. Policies react to past cases (feedback), so the regressors are predetermined rather than strictly exogenous. Infections and policies of neighbouring units contaminate the comparison (contagion spillover). And the direct effect on the reproduction number is not the number of cases averted.

The method under study places parallel trends on log transmission in a multiplicative renewal model, E[Y(i,t) | past] = exp(α(i) + γ(t) + τ D(i,t) + κ D̄(i,t)) · Λ(i,t; φ). It estimates (τ, κ, φ) by quasi-differenced GMM (QD-GMM, method M0). The quasi-difference removes the unit effect α(i) while using only sequentially exogenous instruments. Here τ is the direct log effect of a unit's own policy on its reproduction number, κ the effect of the neighbours' policy share, φ the relative infectiousness of neighbours' recent cases, and ψ the cumulative cases averted.

The experiment asks four questions, fixed in advance in `theory/simulation_spec.md`:

1. Does M0 give approximately unbiased estimates of τ, κ and φ, with valid cluster-robust standard errors and t(G−1) intervals, under feedback, spillover, overdispersion and small counts?
2. Is inference on ψ valid for the realised-path target ψ†?
3. How biased are the standard alternatives (Poisson fixed effects M1, its split-panel bias correction M1-BC, two-step log R DiD M2, TWFE on log(Y+1) M3, pooled renewal Poisson M4), and does the direction match the theory?
4. What are the size and power of the Hansen J test, and how robust is M0 to the misspecifications X1–X5?

The pre-specified support criterion is that, in feedback scenarios, M0 has |bias| ≤ 2 MCSE and coverage in [0.93, 0.97] while M2/M3 fail. The contribution would be weakened if M1 and M1-BC were also unbiased with correct coverage.

## Design

**Data-generating mechanism.** There is one generator (`src/dgp.py`, a renewal model with network pressure). It has two generation-interval lags w = (0.85, 0.15) and negative-binomial innovations with size k = 10, or Poisson in part of the null block. Counties are random points with a symmetrised 6-nearest-neighbour graph, grouped into G = 50 contiguous "states" (the clusters).

Adoption is absorbing and takes one of two forms:

- *random*, where the state adoption week is uniform;
- *feedback*, a Texas-like rule in which a state adopts once its lagged incidence crosses a propensity-dependent threshold, and a county within an adopting state enforces once its own lagged per-capita cases reach 12.

The design has 25 scenarios (`src/scenarios.py`):

- **Core factorial** (1–16): N ∈ {500, 2000} × T ∈ {12, 26} × adoption {random, feedback} × spillover {(φ0, κ0) = (0.4, −0.1), (0, 0)}, with τ0 = log 0.8.
- **Null block** (17–20): τ0 = 0, with adoption {random, feedback} × k ∈ {10, ∞}.
- **Misspecification block** (21–25):
  - X1: wrong generation interval.
  - X2: reporting drop after adoption.
  - X3: heterogeneous τ(i).
  - X4: susceptible depletion.
  - X5: dynamic effects.

**Methods.** Every replication fits the same data set with:

- M0;
- M0-ET, the event-time version with six bins;
- M0-OI, two-step GMM with two extra instruments, which provides the J test;
- M0′, with κ = φ = 0 imposed;
- M1, M1-BC, M2, M3 and M4;
- an infeasible Oracle with known α(i).

M0-G, a global root search over φ ∈ [0, 3] with selection by the over-identified M0-OI criterion, was added before the full run because a pilot showed multiple roots near φ0 = 0 (`experiment_design.md`, deviation 2). **M0-GI** was added *after* interim results had been seen and is labelled post hoc everywhere (deviation 7). It uses the same candidates as M0-G but selects the root by moments built from the φ-instrument at the other reference values 0.05 and 2.0. It was run by a separate runner, `src/simulate_supp.py`, on identical seeds and data sets, and its rows are merged by (scenario, replication) after a data-fingerprint check.

**Performance measures.** For each scenario, method and parameter the study reports:

- bias, empirical SE, mean model SE, relative SE error and RMSE;
- coverage of the 95 % t(G−1) interval, its width and the rejection rate of H0: θ = 0;
- the convergence rate, with failures counted and never replaced by another method's estimate.

Each measure has its Monte Carlo SE (Morris, White & Crowther 2019). Targets are the generating values. X3 has two targets: the per-replication switcher mean and the first-order pseudo-true value. ψ is evaluated against the realised-path ψ† and against a population ψ obtained from 5,000 common-random-number pairs per scenario.

**Replications and compute.** The target was R = 8000 per scenario, chosen in `run_plan.json` from a benchmark. Replications were run in three waves (200 → 1000 → 8000), so a run stopped early still leaves a balanced design. Seeds are 20261112 + 100000·scenario + replication.

**Application.** The application uses weekly US county cases (JHU CSSE) and state/territorial public mask orders (CDC dataset 62d6-pm5i), with clusters G = 49 states. The network is the 2010 Census county adjacency, with ACS commuting flows as a sensitivity analysis, and the period runs from June to November 2020 (`data_contract.md`).

The pre-specified sensitivity analyses are:

- the Nishiura serial interval;
- a w1 grid;
- the commuting network;
- policy lagged one week;
- ω × 0.25 and × 4;
- φ̄ ∈ {0.25, 1};
- excluding Kansas and Texas;
- a lead-indicator placebo.

## Method

M0 solves the just-identified moment conditions E[s(i,t) u(i,t; θ)] = 0, where u(i,t; θ) = Y(i,t+1) m(i,t) / m(i,t+1) − Y(i,t) is the multiplicative quasi-difference. The instruments are s · (week dummies, D, D̄, φ̄ N/Λ(φ̄)), with reference value φ̄ = 0.5. Here N(i,t) = Σ over s of w(s) Σ over j of W(i,j) Y(j,t−s) is the neighbour input to the pressure Λ(i,t; φ) = L(i,t) + φ N(i,t) + ω(i). These are functions of the past only, so they are valid under feedback (A2).

The variance is a cluster-robust sandwich with the G/(G−1) correction, and intervals use t(G−1). ψ is computed by the ν recursion (counterfactual renewal with the policy index set to zero), with a delta-method SE.

The implementation is `src/qdgmm.py`. The comparators are in `src/comparators.py`, ψ in `src/psi.py`, and the root-selection variants are `qdgmm.fit_global(selection="oi")` (M0-G) and `qdgmm.fit_global(selection="inv")` (M0-GI).

## Results

All numbers below are copied from `metrics_key_numbers.md`, which `src/report_numbers.py` writes from `results.json` in the same `run_all.sh` call that produces the tables and figures. The main study reached the target R = 8000 in all 25 scenarios (200,000 replication rows, run complete = true). The post-hoc M0-GI rows come from the supplementary runner, which also reached R = 8000 in every scenario (200,000 rows merged, 0 fingerprint mismatches; `tables/table_8_m0gi_posthoc.tex`).

Table 1. Estimation of τ in the core and null blocks: range over scenarios of the bias (mean of estimate minus generating value), coverage of the nominal 95 % t(G−1) interval, rejection rate of H0: τ = 0 and share of replications with an estimate; R = 8000 replications per scenario.

| block | method | R | bias | coverage | reject H0 | converged |
|---|---|---|---|---|---|---|
| core, spillover (1-4, 9-12) | M0 | 8000 | -0.002 to 0.001 | 94.5 % to 95.4 % | 77.1 % to 100.0 % | 100.0 % to 100.0 % |
| core, spillover (1-4, 9-12) | M0-G | 8000 | -0.002 to 0.001 | 94.5 % to 95.4 % | 77.1 % to 100.0 % | 100.0 % to 100.0 % |
| core, spillover (1-4, 9-12) | M0-GI (post hoc) | 8000 | -0.002 to 0.001 | 94.5 % to 95.4 % | 77.1 % to 100.0 % | 100.0 % to 100.0 % |
| core, spillover (1-4, 9-12) | M0′ | 8000 | -0.248 to -0.114 | 0.0 % to 22.4 % | 100.0 % to 100.0 % | 100.0 % to 100.0 % |
| core, spillover (1-4, 9-12) | M1 | 8000 | -0.027 to -0.017 | 78.1 % to 88.2 % | 99.8 % to 100.0 % | 100.0 % to 100.0 % |
| core, spillover (1-4, 9-12) | M1-BC | 8000 | -0.034 to 0.082 | 25.7 % to 82.7 % | 87.9 % to 100.0 % | 14.7 % to 100.0 % |
| core, spillover (1-4, 9-12) | M2 | 8000 | -0.110 to 0.013 | 0.0 % to 94.0 % | 100.0 % to 100.0 % | 100.0 % to 100.0 % |
| core, spillover (1-4, 9-12) | M3 | 8000 | -0.563 to -0.030 | 0.0 % to 91.6 % | 100.0 % to 100.0 % | 100.0 % to 100.0 % |
| core, no spillover (5-8, 13-16) | M0 | 8000 | 0.132 to 0.469 | 0.6 % to 38.8 % | 9.1 % to 46.9 % | 92.9 % to 100.0 % |
| core, no spillover (5-8, 13-16) | M0-G | 8000 | 0.039 to 0.202 | 52.9 % to 71.0 % | 43.8 % to 86.6 % | 100.0 % to 100.0 % |
| core, no spillover (5-8, 13-16) | M0-GI (post hoc) | 8000 | 0.000 to 0.052 | 87.5 % to 95.2 % | 59.2 % to 100.0 % | 100.0 % to 100.0 % |
| core, no spillover (5-8, 13-16) | M0′ | 8000 | -0.000 to 0.002 | 94.7 % to 95.5 % | 93.2 % to 100.0 % | 100.0 % to 100.0 % |
| core, no spillover (5-8, 13-16) | M1 | 8000 | -0.048 to -0.008 | 51.3 % to 91.1 % | 90.8 % to 100.0 % | 100.0 % to 100.0 % |
| core, no spillover (5-8, 13-16) | M1-BC | 8000 | -0.063 to 0.006 | 30.2 % to 82.0 % | 88.6 % to 100.0 % | 96.4 % to 100.0 % |
| core, no spillover (5-8, 13-16) | M2 | 8000 | -0.119 to 0.048 | 0.0 % to 93.5 % | 99.3 % to 100.0 % | 100.0 % to 100.0 % |
| core, no spillover (5-8, 13-16) | M3 | 8000 | -0.530 to 0.150 | 0.2 % to 87.9 % | 23.5 % to 100.0 % | 100.0 % to 100.0 % |
| null (17-20) | M0 | 8000 | -0.001 to 0.001 | 94.8 % to 95.3 % | 4.7 % to 5.2 % | 100.0 % to 100.0 % |
| null (17-20) | M0-G | 8000 | -0.001 to 0.001 | 94.8 % to 95.3 % | 4.7 % to 5.2 % | 100.0 % to 100.0 % |
| null (17-20) | M0-GI (post hoc) | 8000 | -0.001 to 0.001 | 94.8 % to 95.3 % | 4.7 % to 5.2 % | 100.0 % to 100.0 % |
| null (17-20) | M0′ | 8000 | -0.065 to -0.031 | 37.9 % to 77.0 % | 23.1 % to 62.1 % | 100.0 % to 100.0 % |
| null (17-20) | M1 | 8000 | -0.007 to 0.002 | 90.3 % to 93.7 % | 6.3 % to 9.7 % | 100.0 % to 100.0 % |
| null (17-20) | M1-BC | 8000 | -0.010 to 0.002 | 71.4 % to 85.3 % | 14.7 % to 28.6 % | 99.7 % to 100.0 % |
| null (17-20) | M2 | 8000 | -0.116 to -0.052 | 0.0 % to 23.9 % | 76.1 % to 100.0 % | 100.0 % to 100.0 % |
| null (17-20) | M3 | 8000 | -0.122 to -0.057 | 38.1 % to 88.1 % | 11.9 % to 61.9 % | 100.0 % to 100.0 % |

**Scenarios with spillover.** (φ0 = 0.4; Figure 1, `tables/table_1_core_tau.tex`.) M0 is close to unbiased for τ (bias −0.002 to 0.001) with coverage 94.5 % to 95.4 %, under both random and feedback adoption. The empirical SE (0.037 to 0.082) matches the mean model SE (0.037 to 0.081), and power against H0: τ = 0 is 77.1 % to 100 %. κ and φ are recovered as well (`tables/table_2_core_spill.tex`, Figures 2 and 3). M0-G and M0-GI pick the same root as M0 in every one of these replications.

Under the strict pre-specified criterion (|bias| ≤ 2 MCSE and coverage in [0.93, 0.97]), M0 passes for τ in seven of the eight spillover scenarios. It fails narrowly in scenario 1 (bias −0.0017, bias/MCSE −2.18). For φ, the bias/MCSE ratio exceeds 2 in scenarios 1, 2, 3, 9 and 10 and in all four null scenarios, with a maximum of 5.43 (scenario 9). These biases are at most 0.0035, below 1 % of φ0 = 0.4. They are detectable only because R = 8000 makes the MCSE about 0.0006, and they shrink with N. Coverage stays within [0.93, 0.97] everywhere, and the relative error of the model SE is between −3.5 % and +0.9 % for τ, κ and φ in these scenarios. We report this as a small finite-sample bias of the φ estimate, not as a pass.

The comparators behave as follows:

- M1 (Poisson FE with the same pressure offset) has bias −0.027 to −0.017 and coverage 78.1 % to 88.2 %.
- M1-BC returns an estimate in only 14.7 % (scenario 2) and 23.9 % (scenario 4) of replications. The split-panel jackknife needs within-unit policy variation in both half panels, and feedback adoption with T = 26 often puts every switch in one half. Where it converges, its coverage is 25.7 % to 82.7 %.
- M2 (two-step log R DiD) has bias −0.110 to 0.013 and coverage from 0 % to 94.0 %.
- M3 (TWFE on log(Y+1)) has bias −0.563 to −0.030.
- M4 (pooled renewal Poisson) has bias 0.008 to 0.139.
- M0′, which ignores spillover, has bias −0.248 to −0.114.
- The infeasible Oracle has coverage 91.3 % to 94.8 %.

**Null block.** (Figure 4, `tables/table_3_null.tex`.) The rejection rate of the 5 % t(G−1) test of τ = 0 is 4.7 % to 5.2 % for M0, 5.7 % to 6.4 % for M0-OI and 6.3 % to 9.7 % for M1. It is 14.7 % to 28.6 % for M1-BC, 76.1 % to 100 % for M2, and 11.9 % to 61.9 % for M3.

**Scenarios without spillover.** (φ0 = 0; Figures 1, 2 and 5, `tables/table_6_roots.tex`.) Here M0 fails. Its τ bias is 0.132 to 0.469 and its coverage 0.6 % to 38.8 %. The φ estimate is biased upward by 0.232 to 0.888 and the κ estimate downward by 0.150 to 0.457.

The cause is the estimating equations, not the code. In about half of the replications (47.2 % to 51.7 %) the profiled φ-moment has more than one root in [0, 3], and the local solver started at φ̄ = 0.5 converges to a spurious interior root. The M0-OI Hansen J test detects this in only 6.4 % to 27.2 % of replications.

The pre-specified global search M0-G reduces the τ bias to 0.039 to 0.202 but leaves coverage at 52.9 % to 71.0 %. It selects the boundary φ = 0 in 28.6 % to 44.8 % of replications.

The post-hoc selection M0-GI gives τ bias 0.000 to 0.052 and coverage 87.5 % to 95.2 %, and selects the boundary in 63.7 % to 80.9 % of replications. M0′, which imposes κ = φ = 0 and is therefore correctly specified in these scenarios, has bias −0.000 to 0.002 and coverage 94.7 % to 95.5 %. The Oracle also under-covers in some of these scenarios (80.2 % to 94.5 %).

**Misspecification.** (Figure 6, `tables/table_4_misspec.tex`.)

- X1 (wrong generation interval): M0 τ bias 0.010, coverage 94.6 %; J rejects in 100 % of replications.
- X2 (reporting drop at adoption): τ bias −0.057, coverage 59.4 %. M1 is worse, with bias −0.099 and coverage 0.1 %. J rejects in 100 %.
- X3 (heterogeneous τ(i)): bias 0.005 and coverage 94.0 % against both the switcher mean and the first-order pseudo-true value.
- X4 (susceptible depletion): bias −0.000, coverage 95.2 %.
- X5 (dynamic effects): the constant-τ M0 is not designed for this case and has bias 0.219 and coverage 0 %. The event-time version M0-ET recovers every bin, with bias −0.000 and coverage 95.3 % to 95.8 %. J has no power against X5 (rejection 4.7 %).

**Cases averted ψ.** (`tables/table_5_psi.tex`.) In the spillover scenarios, the median relative bias of the ψ estimate for the treated units is −0.004 to 0.006. Coverage of the realised-path target ψ† is 85.7 % to 95.4 %, and it is below 93 % in scenarios 1, 2, 4 and 10 (92.4 %, 85.7 %, 91.9 %, 90.1 %). Coverage of the population ψ is 72.1 % to 92.5 %.

Without spillover, the wrong root carries over to ψ, with median relative bias up to 2.474 (scenario 14). Under misspecification, ψ is fragile:

- X1: relative bias 3.386.
- X2: relative bias 1004.524, coverage 50.8 %.
- X3: relative bias 0.779.
- X4: population coverage 53.5 %.
- X5: relative bias −0.994, coverage 0 %.

**Application.** (Figure 7, `tables/table_7_application.tex`.) The panel has 3091 counties in G = 49 states over T = 26 weeks (2020-06-01 to 2020-11-23). The median weekly count is 27, 7.4 % of county-weeks are zero, 56.8 % of county-weeks are under a mask order, and 1768 counties switch during the window.

Table 2. Application to US county COVID-19 cases and mask orders, main specification: point estimates with 95 % t(G−1) confidence intervals (G = 49 states) for the direct effect τ, the neighbour-policy effect κ and the coupling φ; NA = parameter not in the model.

| analysis | method | τ (95 % CI) | κ (95 % CI) | φ (95 % CI) |
|---|---|---|---|---|
| main | M0 | 0.040 (−0.217, 0.296) | −0.506 (−1.018, 0.006) | 0.772 (0.590, 0.954) |
| main | M0-OI | −0.025 (−0.249, 0.199) | −0.517 (−0.983, −0.050) | 0.590 (0.475, 0.705) |
| main | M1 | −0.024 (−0.109, 0.061) | −0.069 (−0.163, 0.025) | 0.567 (0.355, 0.779) |
| main | M2 | −0.084 (−0.205, 0.036) | NA | NA |
| main | M3 | 0.087 (−0.188, 0.361) | NA | NA |
| main | M4 | 0.004 (−0.030, 0.038) | 0.004 (−0.034, 0.043) | 0.027 (0.020, 0.035) |

M0-G and M0-GI select the same root as M0 in the main specification.

- **Over-identification.** The M0-OI Hansen J statistic is 17.873 (df 2, p = 0.0001), so the over-identifying restrictions are rejected.
- **Sensitivity of τ under M0.** Across the pre-specified sensitivity analyses the τ estimate runs from −0.177 (policy lagged one week) to 0.120 (w1 = 1.0), and every 95 % interval contains 0.
- **M0-G selections.** In five sensitivity analyses M0-G selects a different root from M0:
  - the boundary φ = 0 under the one-week lag (the τ estimate −2.178) and ω × 4 (−1.274);
  - an interior root with the φ estimate ≤ 0.010 for w1 = 0.6, 0.7 and 0.9.

  M0-GI agrees with M0 in every sensitivity analysis.
- **Event-time estimates.** The M0-ET bins e0 to e6 are all negative (−0.232 to −0.577), and the joint Wald test of all bins equal to zero gives F = 2.642 (df 7, 42), p = 0.0234. Their average over e0–e5, −0.412 (95 % CI −1.154, 0.330), is descriptive and was not pre-specified.
- **Lead placebo.** The coefficient of the lead indicator is −0.332 (se 0.259).
- **Cases averted.** The ψ estimate summed over all counties is 3.41e5 (se 3.97e5) at 4 weeks and 7.17e10 (se 5.16e11) at 26 weeks. The 26-week value exceeds the 7.67e6 cases observed in treated counties by four orders of magnitude.

![Figure 1: τ bias and coverage in the 16 core scenarios](figures/fig_01_core_tau.png)

Figure 1. Bias (mean ± 2 MCSE) and coverage of 95 % intervals for τ, by method, in the 16 core scenarios (R = 8000 for every method, including the post-hoc M0-GI).

![Figure 2: φ bias and coverage](figures/fig_02_core_phi.png)

Figure 2. Bias and coverage for the transmission coupling φ, including the boundary φ0 = 0.

![Figure 3: κ bias and coverage](figures/fig_03_core_kappa.png)

Figure 3. Bias and coverage for the neighbour-policy effect κ.

![Figure 4: type I error](figures/fig_04_null_typeI.png)

Figure 4. Rejection rate of the 5 % test of τ = 0 under the null (± 2 MCSE).

![Figure 5: roots](figures/fig_05_roots.png)

Figure 5. Share of replications with multiple roots and boundary selections.

![Figure 6: misspecification](figures/fig_06_misspec.png)

Figure 6. (a) τ bias under X1–X5; (b) M0-ET event-time path under X5.

![Figure 7: application](figures/fig_07_application.png)

Figure 7. (a) τ with 95 % t(G−1) intervals by method and sensitivity analysis on the US county data; (b) profiled φ-moment and the M0-G criterion Q in the main specification.

## Analysis

**Aim 1 (M0 under feedback and spillover).** With φ0 away from the boundary, M0 is close to unbiased with nominal coverage under random and feedback adoption, negative-binomial overdispersion and small counts (N = 500, T = 12). Its size under the null is 4.7 % to 5.2 %. Feedback did not degrade M0 relative to random adoption. This matches the sequential-exogeneity argument: the instruments are functions of the past only.

The strict criterion of 2 MCSE is not met for φ in 9 of the 12 scenarios with φ0 > 0, and for τ in one of them. The size of these biases (≤ 0.0035 for φ, 0.0017 for τ) is negligible next to the empirical SE (0.037 or more for τ). It is the expected O(1/T) finite-sample bias of a non-linear GMM estimator. We do not claim exact unbiasedness.

**The boundary φ0 = 0 is a limitation of the method.** When φ0 = 0 the neighbour term carries no signal. The just-identified M0 equations then have spurious interior roots in about half of the replications, and a local solver lands on them. Three observations indicate that this is not an implementation error:

- M0′, the same code path with φ fixed at 0, is unbiased with nominal coverage in these scenarios;
- the grid search finds more than one root in about half of the replications, and M0-GI, which picks among the same candidate roots, removes most of the bias;
- no replication with φ0 = 0.4 has more than one root.

This is a global-identification failure of the estimating equations. Two consequences follow.

- **Reporting.** A reported M0 estimate should come with the φ-profile.
- **Selection.** The pre-specified selection rule (M0-G, min over roots of the M0-OI criterion) only halves the bias, because that criterion is itself flat near φ = 0.

The post-hoc rule M0-GI (87.5 % to 95.2 % coverage) shows that the information to choose the right root is in the data. Because it was designed after seeing interim results, it must be validated on a fresh design before it can replace M0-G.

**Aim 2 (ψ).** Where the θ estimate is consistent, the ψ estimate is close to unbiased for ψ†. However, the delta-method interval under-covers in the T = 26, N = 500 scenarios (85.7 % and 90.1 %), even for the realised-path target, for which the theory predicts nominal coverage. A plausible cause is that the delta method is inaccurate for the long counterfactual products: ψ depends on exp(−τ) compounded over up to 26 weeks of renewal. The population target is under-covered further, as Remark 3 predicts. ψ amplifies any error in the θ estimate (X1, X2, X3, and the no-spillover roots).

**Aim 3 (comparators).**

- *M1.* Poisson FE has a design-dependent negative bias (−0.048 to −0.008 in the core block) and has coverage below 93 % in all 16 core scenarios (51.3 % to 91.1 %). The bias correction M1-BC does not repair this: it is undefined whenever a half panel has no switchers, and otherwise inflates the variance without matching SEs.
- *M2 and M3.* The two-step log R DiD and TWFE on log(Y+1) fail in direction and size, including type I error up to 100 % (M2) and 61.9 % (M3) under the null.

The contribution would have been weakened if M1 and M1-BC were unbiased with correct coverage. They are not, so the comparison supports M0, conditional on φ0 being away from 0. We did not check the predicted decomposition of the M2 bias into its two theoretical components or the predicted sign of the M3 bias scenario by scenario. No claim is made about them.

**Aim 4 (J test and robustness).** J has full power against X1 and X2, low power against the wrong-root problem (6.4 % to 27.2 %), and no power against dynamic effects (X5). M0 is robust to heterogeneous effects (X3) and to susceptible depletion at the attack rates of this design (X4). It is biased under differential reporting (X2); no comparator has better coverage there (M1 0.1 %, M2 0.0 %, M1-BC 17.6 %). Under dynamic effects (X5) the constant-τ summary averages very different event-time effects, so M0-ET should be the default report when effects may build up.

**Application.** On US county data the transmission-scale estimate of the direct mask-order effect is imprecise, and its interval covers both no effect and a reduction in the reproduction number of about 20 %: the τ estimate = 0.040 (−0.217, 0.296). Three findings argue against reading it as "no effect":

- J rejects the model (p = 0.0001).
- The event-time bins are all negative and jointly different from zero, which suggests that effects build up over time. This is the situation in which, in X5, the constant-τ estimate was strongly biased toward zero.
- The ψ estimate grows to 7.17e10 at 26 weeks. This is consistent with a counterfactual renewal without orders that grows without bound over six months, so ψ beyond a few weeks is not informative.

The φ estimate of 0.772 is well away from the boundary, and M0-GI agrees with M0. However, M0-G jumps to near-boundary roots in five sensitivity analyses, so root selection matters in practice as well. The lead placebo cannot be used as a pre-trend test here, because under feedback a future order is correlated with current transmission (Proposition 2(c)).

## Limitations

- **Boundary.** M0 is not reliable at or near φ0 = 0. The pre-specified fix M0-G is only partial, and the rule that works (M0-GI) is post hoc. Its sampling properties were assessed on the same designs that motivated it (R = 8000 in all 25 scenarios). It still has a small positive bias in the N = 500 scenarios (0.015 in scenario 5, 0.052 in scenario 13), where coverage is 91.2 % and 87.5 %.
- **Finite-sample bias.** The φ estimate has a small finite-sample bias that the strict 2-MCSE criterion detects at R = 8000.
- **ψ inference.** The ψ intervals under-cover for T = 26 and N = 500 even for ψ†. ψ is very sensitive to misspecification, and in the application it is not usable beyond short horizons.
- **Design range.** One DGP family was used: a two-lag generation interval, a 6-nearest-neighbour graph, 50 clusters and adoption rules of one type. Results for few clusters (G < 30), long generation intervals or other network shapes are not covered.
- **Comparators not attempted.** The M1-BC comparator follows Dhaene–Jochmans, which assumes stationarity. Its frequent non-existence under staggered adoption is a property of that estimator in this design, not a bug, but other bias corrections were not tried.
- **Theory checks not run.** The M2 bias decomposition and the predicted sign of the M3 bias were not verified.
- **Application data.** The application uses state/territorial mask-order data mapped to counties, reported cases, and a fixed generation interval. The J test rejects the model, so the point estimates should be read as conditional on a misspecified model.
- **Process.** The `run_all.sh` driver stopped with a shell syntax error after the simulation stage, because the file was edited while the running instance was reading it. No rows were lost, and the later stages were rerun from the saved rows.

## Reproduction

All commands are run from `experiment/`.

1. Install the dependencies:
   ```bash
   pip install -r requirements.txt
   ```
2. Fetch the data. This step verifies the SHA-256 checksums in `data/checksums.tsv`:
   ```bash
   bash fetch_data.sh
   bash fetch_fonts.sh
   ```
3. Run everything (deps → data → fonts → prep → tests → sim → supp → app → agg → tables → figures):
   ```bash
   QDG_WORKERS=16 bash run_all.sh
   ```
   The simulation is resumable: rerun the same command after an interruption. `MAX_MINUTES` caps the start of new blocks.
4. Rebuild the outputs from existing rows:
   ```bash
   STAGES="agg tables figures" bash run_all.sh
   ```
   This regenerates `results.json`, `metrics.csv`, `metrics_key_numbers.md`, `tables/` (with `tables/manifest.json`) and `figures/` (with `figures/manifest.json`).
5. Run the tests:
   ```bash
   python3 -m pytest
   ```
   Expected result: 40 tests pass.

Expected wall-clock on 16 cores: about 3.4 h for the main Monte Carlo stage (201.4 min measured), roughly 1.5–2 h for the supplementary M0-GI stage on 16 workers (an extrapolation: only its last 9320 of 20000 tasks were timed on 16 workers, 47.2 min), under 1 min for the application (23 s measured), and about 1 min for aggregation, tables and figures (49 s + 11 s measured). Run sequentially by `run_all.sh`, the total is therefore about 5–5.5 h; in this run the supplementary stage overlapped the main stage at low priority.
