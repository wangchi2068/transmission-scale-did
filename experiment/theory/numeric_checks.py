"""Symbolic and numerical checks of the derivations in Web Appendix B.

Every check is deterministic (fixed seed) and printed as PASS/FAIL with the quantity compared.
Run from the repository root:
  python -I experiment/theory/numeric_checks.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy import optimize

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "experiment" / "code"))

import dgp  # noqa: E402
import qdgmm  # noqa: E402

RESULTS = []


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, bool(ok), detail))
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}")


# ---------------------------------------------------------------- S1 Lemma 2 martingale decomposition
def s1_quasi_difference() -> None:
    """Verify Lemma 2: g_g(theta_0) = sum_{s=1}^T sum_{i in g} a_{is} (Y_{is} - mu_{is}) exactly."""
    rng = np.random.default_rng(20270101)
    T, N, G = 7, 40, 8
    cl = np.repeat(np.arange(G), N // G)
    alpha = rng.normal(0.0, 0.3, size=N)
    gamma = np.r_[0.0, rng.normal(0.0, 0.2, size=T - 1)]
    D = (rng.random((T, N)) > 0.6).astype(float)
    Db = rng.uniform(0.0, 1.0, size=(T, N))
    lam = rng.uniform(5.0, 40.0, size=(T, N))
    m = np.exp(gamma[:, None] - 0.2 * D - 0.1 * Db) * lam
    mu = np.exp(alpha[None, :]) * m
    Y = rng.poisson(mu).astype(float)
    z = rng.normal(size=(3, T - 1, N))

    r = m[:-1] / m[1:]
    u = Y[1:] * r - Y[:-1]
    eps = Y - mu

    g_direct = np.zeros((3, G))
    g_mds = np.zeros((3, G))
    for g in range(G):
        sel = cl == g
        g_direct[:, g] = np.einsum("qtn,tn->q", z[:, :, sel], u[:, sel])
        for s in range(T):
            a = np.zeros((3, sel.sum()))
            if s >= 1:
                a += z[:, s - 1, sel] * r[s - 1, sel][None, :]
            if s <= T - 2:
                a -= z[:, s, sel]
            g_mds[:, g] += a @ eps[s, sel]

    max_err = float(np.max(np.abs(g_direct - g_mds)))
    check("S1 Lemma 2 martingale sum identity", np.allclose(g_direct, g_mds, rtol=1e-12, atol=1e-10),
          f"max |diff|={max_err:.2e}")


# ---------------------------------------------------------------- S3 closed-form profiling (eq:profile)
def s3_profiling() -> None:
    """Verify eq. (4): profiled delta_t(beta) solves all T-1 week-effect moment equations exactly."""
    rng = np.random.default_rng(20270102)
    T, N = 8, 60
    geo = dgp.make_geography(N, G=10, cols=2)
    p = dgp.DGPParams(N=N, T=T, adopt="random")
    sim = dgp.simulate(p, geo, [20270102, 0], [20270102, 1])
    panel = dgp.build_panel(sim["Yall"], sim["Dall"], geo.W, sim["omega"], geo.cl)
    mod = qdgmm.QDModel(panel, qdgmm.spec_m0(panel))
    beta = np.array([-0.18, -0.07, 0.35])
    m = mod.moments(beta)
    week_moments = m["gi"][:T - 1].sum(1)
    max_res = float(np.max(np.abs(week_moments)))
    # Also check eq. (4) against direct log-ratio formula
    c = mod._core(beta)
    q = np.exp(beta[0] * (panel["D"][:-1] - panel["D"][1:]) +
               beta[1] * (panel["Db"][:-1] - panel["Db"][1:])) * c["lam"][:-1] / c["lam"][1:]
    ed_formula = (mod.s * panel["Y"][1:] * q).sum(1) / (mod.s * panel["Y"][:-1]).sum(1)
    check("S3 eq. (4) profiled delta solves week moments",
          max_res < 1e-10 and np.allclose(np.exp(m["delta"]), ed_formula, rtol=1e-12),
          f"max |week moment|={max_res:.2e}")


# ---------------------------------------------------------------- S4 expansion of b_it
def s4_b_expansion() -> None:
    from scipy.special import gammaln

    def nb_log_mean_offset(mu, k, c):
        q = 1.0 / (1.0 + mu / k)
        ks = np.arange(0, 400)
        pmf = np.exp(gammaln(ks + k) - gammaln(k) - gammaln(ks + 1) + k * np.log(q) + ks * np.log1p(-q))
        pmf = pmf / pmf.sum()
        return float(np.sum(pmf * np.log(ks + c)) - np.log(mu))

    k, c = 10.0, 0.5
    exact2 = nb_log_mean_offset(2.0, k, c)
    exact5 = nb_log_mean_offset(5.0, k, c)
    approx5 = (c - 0.5 + c / k) / 5.0 - 1.0 / (2.0 * k)
    ok = abs(exact2 - (-0.00556)) < 2e-4 and abs(exact5 - approx5) < 0.012
    check("S4 expansion accurate for mu>=5 and exact value at mu=2",
          ok, "mu=2 exact %.5f (paper -0.00556) | mu=5 exact %.5f approx %.5f"
          % (exact2, exact5, approx5))


# ---------------------------------------------------------------- S5 working variance identity
def s5_working_variance() -> None:
    """Verify Remark 1: r^2 mu_{t+1} = r mu_t and Poisson Var(u_it | F_{t-1}) = mu_it (1 + E[r_it | F_{t-1}])."""
    rng = np.random.default_rng(20270103)
    mu1 = 18.0
    m1 = 9.0
    alpha = np.log(mu1 / m1)
    # Two possible future policy states at t+1 (e.g., D=0 or D=1) with equal conditional probability
    m2_vals = np.array([12.0, 8.0])
    r_vals = m1 / m2_vals
    mu2_vals = np.exp(alpha) * m2_vals
    check("S5 algebraic identity r^2 mu_{t+1} = r mu_t",
          np.allclose(r_vals ** 2 * mu2_vals, r_vals * mu1, rtol=1e-12))

    # Monte Carlo check of conditional variance Var(Y_{t+1} r_t - Y_t | F_{t-1}) under Poisson law
    n = 400_000
    y1 = rng.poisson(mu1, size=n).astype(float)
    idx = rng.integers(0, 2, size=n)
    r = r_vals[idx]
    mu2 = mu2_vals[idx]
    y2 = rng.poisson(mu2).astype(float)
    u = y2 * r - y1
    emp_var = float(np.var(u, ddof=1))
    theo_var = float(mu1 * (1.0 + np.mean(r_vals)))
    rel_err = abs(emp_var - theo_var) / theo_var
    check("S5 Remark 1 Poisson working variance Var(u|F_{t-1}) = mu(1 + E[r|F_{t-1}])",
          rel_err < 0.01, f"emp={emp_var:.3f}, theo={theo_var:.3f}, rel_err={rel_err:.4f}")


# ---------------------------------------------------------------- S7 Proposition 3 heterogeneity weights
def s7_weights() -> None:
    """Verify Proposition 3: constant-effect QD-GMM root lies in [min tau_i, max tau_i] and matches 1st-order avg."""
    rng = np.random.default_rng(20270105)
    N, T = 4000, 6
    tau_i = rng.uniform(-0.45, -0.05, size=N)
    alpha_i = rng.normal(0.0, 0.2, size=N)
    gamma_t = np.linspace(0.0, 0.25, T)
    # Absorbing staggered adoption
    adopt_t = rng.integers(1, T + 1, size=N)
    D = (np.arange(T)[:, None] >= adopt_t[None, :]).astype(float)
    lam = rng.uniform(10.0, 40.0, size=(T, N))
    s = 1.0 / (lam[:-1] + 1.0)
    mu = np.exp(alpha_i[None, :] + gamma_t[:, None] + tau_i[None, :] * D) * lam
    sw = D[1:] - D[:-1]

    M_t = (s * mu[:-1]).mean(1)
    A_t = (s * D[:-1] * mu[:-1]).mean(1)
    K_t = (s * mu[:-1] * sw).mean(1)
    H_t = (s * mu[:-1] * sw * tau_i[None, :]).mean(1)

    def h_prof(tau):
        C_t = (s * mu[:-1] * sw * (np.exp(tau_i[None, :] - tau) - 1.0)).mean(1)
        return float(-np.sum(A_t * C_t / (M_t + C_t)))

    tau_tilde = float(optimize.brentq(h_prof, tau_i.min(), tau_i.max()))
    w = np.where(M_t > 0, A_t / M_t, 0.0)
    tau_first_order = float(np.sum(w * H_t) / np.sum(w * K_t))
    check("S7 Proposition 3 root in support & close to non-negative weighted average",
          tau_i.min() <= tau_tilde <= tau_i.max() and abs(tau_tilde - tau_first_order) < 0.01,
          f"support=[{tau_i.min():.3f},{tau_i.max():.3f}] root={tau_tilde:.4f} first_order={tau_first_order:.4f}")


# ---------------------------------------------------------------- S8 T=2 Poisson FE expected score (eq:t2fe)
def s8_t2fe() -> None:
    """Verify eq. (6) for T=2 under BOTH randomised adoption and threshold feedback adoption."""
    rng = np.random.default_rng(20270104)
    n = 200_000
    lam1 = rng.uniform(5, 50, n)
    omega = rng.uniform(0.5, 5, n)
    y1 = rng.poisson(lam1).astype(float)
    lam2_obs = y1 + omega
    tau0 = float(np.log(0.8))
    h0 = 0.0

    def lhs(h, d_val, D_vec):
        # eq. (6): E[1(D_{i2}=d) * Lambda_{i2} * (Lambda_{i1} exp(h0 + tau0*d) - exp(h) Y_{i1}) / (Lambda_{i1} + exp(h) Lambda_{i2})]
        ind = (D_vec == d_val)
        return float(np.mean(ind * lam2_obs * (lam1 * np.exp(h0 + tau0 * d_val) - np.exp(h) * y1)
                             / (lam1 + np.exp(h) * lam2_obs)))

    # (a) Randomised adoption: D_{i2} independent of Y_{i1}
    D_rand = rng.binomial(1, 0.5, size=n).astype(float)
    lhs_r0 = lhs(h0, 0.0, D_rand)
    lhs_r1 = lhs(h0 + tau0, 1.0, D_rand)
    h0_rand = float(optimize.brentq(lambda h: lhs(h, 0.0, D_rand), -1.0, 1.0))
    h1_rand = float(optimize.brentq(lambda h: lhs(h, 1.0, D_rand), -1.0, 1.0))
    check("S8(a) randomised T=2 eq. (6) negative at true h for BOTH d=0 and d=1",
          lhs_r0 < 0.0 and lhs_r1 < 0.0 and abs((h1_rand - h0_rand) - tau0) > 1e-3,
          f"LHS(d=0)={lhs_r0:+.3e}, LHS(d=1)={lhs_r1:+.3e}, h1-h0={h1_rand - h0_rand:+.4f} vs tau0={tau0:+.4f}")

    # (b) Feedback adoption: D_{i2} = 1(Y_{i1} > 30)
    D_fb = (y1 > 30.0).astype(float)
    lhs_f0 = lhs(h0, 0.0, D_fb)
    lhs_f1 = lhs(h0 + tau0, 1.0, D_fb)
    h0_fb = float(optimize.brentq(lambda h: lhs(h, 0.0, D_fb), -1.0, 1.0))
    h1_fb = float(optimize.brentq(lambda h: lhs(h, 1.0, D_fb), -1.0, 1.0))
    check("S8(b) feedback T=2 eq. (6): d=0 positive, d=1 negative, downward bias",
          lhs_f0 > 0.0 and lhs_f1 < 0.0 and (h1_fb - h0_fb) < tau0 - 0.05,
          f"LHS(d=0)={lhs_f0:+.3e}, LHS(d=1)={lhs_f1:+.3e}, h1-h0={h1_fb - h0_fb:+.4f} vs tau0={tau0:+.4f}")


def main() -> int:
    s1_quasi_difference()
    s3_profiling()
    s4_b_expansion()
    s5_working_variance()
    s7_weights()
    s8_t2fe()
    ok = all(r[1] for r in RESULTS)
    print(f"\n{sum(r[1] for r in RESULTS)}/{len(RESULTS)} checks passed")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
