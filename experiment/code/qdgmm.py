"""Quasi-differenced sequential-moment GMM (QD-GMM) on the multiplicative renewal panel (theory/methods.tex).

Model:   E(Y_it | F_{t-1}, alpha_i) = exp(alpha_i + gamma_t + b'X_it) * Lambda_it(phi),
         Lambda_it(phi) = L_it + phi * N_it + omega_i.
Moment:  u_it = Y_{i,t+1} m_it / m_{i,t+1} - Y_it,   E(z_it u_it) = 0 for F_{t-1}-measurable z_it.

The linear index b'X is general: (tau D, kappa Dbar) for M0, event-time bins for M0-ET, D only for M0'.
Week effects delta_t = gamma_{t+1} - gamma_t are profiled in closed form (eq:profile); the remaining
K (+1 if phi is free) equations are solved with Powell's hybrid method using the exact profiled Jacobian.
Variance: V = G/(G-1) * Gamma^{-1} Omega Gamma^{-T} / G with the full (delta, beta) Jacobian; t_{G-1} intervals.
Over-identified version (M0-OI): two-step GMM over (delta, beta) with weight Omega_tilde^{-1} at the
just-identified estimate; Hansen J ~ chi^2_{q-p} (Corollary 1).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np
from scipy import optimize, stats


@dataclass
class QDSpec:
    lin: Dict[str, np.ndarray]                       # name -> T x N regressor in the linear index
    phi_free: bool = True
    phi_fixed: float = 0.0
    phi_bar: float = 0.5                             # reference value for the scale weight / phi instrument
    lin_inst: Optional[Dict[str, np.ndarray]] = None # name -> T x N instrument (default: the regressor)
    extra_inst: List[np.ndarray] = field(default_factory=list)  # (T-1) x N raw extra instruments (times s)
    extra_names: List[str] = field(default_factory=list)
    G: Optional[int] = None                          # clusters of the panel the spec was built on

    @property
    def names(self) -> List[str]:
        return list(self.lin) + (["phi"] if self.phi_free else [])


class QDModel:
    """Pre-computes the static arrays of one panel/specification."""

    def __init__(self, panel: dict, spec: QDSpec):
        self.p = panel
        self.spec = spec
        Y = panel["Y"]
        self.Y = Y
        self.T, self.N = Y.shape
        self.L, self.Nb, self.om = panel["L"], panel["Nb"], panel["omega"]
        self.cl, self.G = panel["cl"], int(panel["G"])
        self.X = [np.asarray(v, float) for v in spec.lin.values()]
        self.K = len(self.X)
        lam_bar = self.L + spec.phi_bar * self.Nb + self.om
        self.s = 1.0 / (lam_bar[:-1] + 1.0)
        inst = spec.lin_inst or {}
        Z = [self.s * np.asarray(inst.get(nm, x), float)[:-1] for nm, x in spec.lin.items()]
        if spec.phi_free:
            Z.append(self.s * spec.phi_bar * self.Nb[:-1] / lam_bar[:-1])
        self.Zpar = Z                                  # just-identified parameter instruments
        self.Zext = [self.s * np.asarray(z, float) for z in spec.extra_inst]
        self.dX = [x[:-1] - x[1:] for x in self.X]
        self.np_ = self.K + int(spec.phi_free)
        # scale for convergence tolerance
        self.fscale = max(float(np.sum(self.s * Y[:-1])) / self.N, 1e-12)

    # ------------------------------------------------------------------------------------------
    def _split(self, beta):
        beta = np.asarray(beta, float)
        b = beta[:self.K]
        phi = beta[self.K] if self.spec.phi_free else self.spec.phi_fixed
        return b, phi

    def _core(self, beta, delta=None):
        """Return dict with v, Yr, dx list, Ed (exp delta)."""
        b, phi = self._split(beta)
        lam = self.L + phi * self.Nb + self.om
        if np.any(lam <= 0) or not np.all(np.isfinite(beta)):
            return None
        eta = np.zeros_like(self.Y)
        for bk, x in zip(b, self.X):
            if bk != 0.0:
                eta = eta + bk * x
        lr = np.clip(eta[:-1] - eta[1:], -50, 50)
        q = np.exp(lr) * lam[:-1] / lam[1:]
        Yq = self.Y[1:] * q
        if delta is None:
            den = (self.s * self.Y[:-1]).sum(1)
            num = (self.s * Yq).sum(1)
            if np.any(den <= 0) or np.any(num <= 0):
                return None
            Ed = num / den
        else:
            Ed = np.exp(np.asarray(delta, float))
        Yr = Yq / Ed[:, None]
        v = Yr - self.Y[:-1]
        dx = list(self.dX)
        if self.spec.phi_free:
            r = self.Nb / lam
            dx.append(r[:-1] - r[1:])
        return dict(v=v, Yr=Yr, dx=dx, Ed=Ed, lam=lam)

    def _blocks(self, c, Zlist):
        """Sums for the Jacobian: a_t, b_tj, c_mt, H_mj (see module doc)."""
        sYr = self.s * c["Yr"]
        a = sYr.sum(1)                                             # (T-1,)
        b = np.stack([(sYr * d).sum(1) for d in c["dx"]], 1)       # (T-1, p)
        cm = np.stack([(z * c["Yr"]).sum(1) for z in Zlist], 0)    # (m, T-1)
        H = np.array([[(z * c["Yr"] * d).sum() for d in c["dx"]] for z in Zlist])  # (m, p)
        return a, b, cm, H

    def profiled_eq(self, beta):
        c = self._core(beta)
        if c is None:
            big = 1e6 * self.fscale
            return np.full(self.np_, big), np.eye(self.np_)
        F = np.array([(z * c["v"]).sum() for z in self.Zpar]) / self.N
        a, b, cm, H = self._blocks(c, self.Zpar)
        Jp = (H - cm @ (b / a[:, None])) / self.N
        return F, Jp

    # ------------------------------------------------------------------------------------------
    def solve(self, x0, xtol=1e-10, maxfev=400):
        sol = optimize.root(self.profiled_eq, np.asarray(x0, float), jac=True, method="hybr",
                            options={"xtol": xtol, "maxfev": maxfev})
        F, _ = self.profiled_eq(sol.x)
        ok = bool(np.all(np.isfinite(sol.x)) and np.max(np.abs(F)) < 1e-7 * self.fscale)
        return sol.x, ok, float(np.max(np.abs(F)) / self.fscale)

    def moments(self, beta, delta=None, extra=False):
        """Unit contributions g_i (q x N), full Jacobian sum J (q x p_full) and delta."""
        c = self._core(beta, delta)
        if c is None:
            return None
        Tm1 = self.T - 1
        Zl = self.Zpar + (self.Zext if extra else [])
        gi = np.vstack([self.s * c["v"], np.stack([(z * c["v"]).sum(0) for z in Zl], 0)])
        a, b, cm, H = self._blocks(c, Zl)
        P = Tm1 + self.np_
        Q = Tm1 + len(Zl)
        J = np.zeros((Q, P))
        J[np.arange(Tm1), np.arange(Tm1)] = -a
        J[:Tm1, Tm1:] = b
        J[Tm1:, :Tm1] = -cm
        J[Tm1:, Tm1:] = H
        return dict(gi=gi, J=J, delta=np.log(c["Ed"]), c=c)

    def cluster_sum(self, gi):
        return np.vstack([np.bincount(self.cl, r, self.G) for r in gi])   # (q, G)

    def inference(self, beta):
        m = self.moments(beta)
        if m is None:
            return None
        G = self.G
        gg = self.cluster_sum(m["gi"])
        Gam = m["J"] / G
        Om = gg @ gg.T / G
        Gi = np.linalg.inv(Gam)
        V = Gi @ Om @ Gi.T / G * G / (G - 1)
        Tm1 = self.T - 1
        Vb = V[Tm1:, Tm1:]
        return dict(V=V, Vb=Vb, delta=m["delta"], Gam=Gam, Om=Om, gg=gg, cond=float(np.linalg.cond(Gam)))


def tcrit(G, level=0.95):
    return float(stats.t.ppf(0.5 + level / 2, G - 1))


def fit(panel: dict, spec: QDSpec, starts: Optional[Sequence[Sequence[float]]] = None) -> dict:
    """Just-identified QD-GMM. Returns est/se dicts, ok flag, V_beta, starts tried."""
    mod = QDModel(panel, spec)
    names = spec.names
    if starts is None:
        base = [0.0] * mod.K + ([spec.phi_bar] if spec.phi_free else [])
        alt1 = [-0.2] * mod.K + ([0.2] if spec.phi_free else [])
        alt2 = [0.0] * mod.K + ([1.0] if spec.phi_free else [])
        starts = [base, alt1, alt2]
    best = None
    for k, x0 in enumerate(starts):
        x, ok, res = mod.solve(x0)
        if best is None or (ok and not best[1]) or (ok == best[1] and res < best[2]):
            best = (x, ok, res, k)
        if ok:
            break
    x, ok, res, k = best
    out = dict(names=names, ok=ok, resid=res, start_used=k, beta=x)
    inf = mod.inference(x) if np.all(np.isfinite(x)) else None
    if inf is None or not np.all(np.isfinite(inf["Vb"])) or np.any(np.diag(inf["Vb"]) < 0):
        out["ok"] = False
        out.update(est={n: float(v) for n, v in zip(names, x)}, se={n: float("nan") for n in names}, Vb=None)
        return out
    se = np.sqrt(np.diag(inf["Vb"]))
    out.update(est={n: float(v) for n, v in zip(names, x)}, se={n: float(v) for n, v in zip(names, se)},
               Vb=inf["Vb"], delta=inf["delta"], cond=inf["cond"], model=mod)
    return out


def fit_two_step(panel: dict, spec: QDSpec, first: dict) -> dict:
    """Two-step GMM over (delta, beta) with extra instruments; Hansen J test (Corollary 1)."""
    mod = QDModel(panel, spec)
    names = spec.names
    n_extra = len(mod.Zext)
    out = dict(names=names, ok=False, J=float("nan"), J_df=n_extra, J_p=float("nan"))
    if not first.get("ok") or n_extra == 0:
        return out
    beta1 = np.array([first["est"][n] for n in names])
    m1 = mod.moments(beta1, extra=True)
    if m1 is None:
        return out
    G = mod.G
    gg1 = mod.cluster_sum(m1["gi"])
    Om1 = gg1 @ gg1.T / G
    try:
        Wm = np.linalg.inv(Om1)
        C = np.linalg.cholesky(Wm)          # W = C C^T
    except np.linalg.LinAlgError:
        return out
    Tm1 = mod.T - 1
    theta1 = np.concatenate([m1["delta"], beta1])

    def resid(th):
        m = mod.moments(th[Tm1:], delta=th[:Tm1], extra=True)
        if m is None:
            return np.full(C.shape[0], 1e6)
        gbar = mod.cluster_sum(m["gi"]).sum(1) / G
        return np.sqrt(G) * (C.T @ gbar)

    def jac(th):
        m = mod.moments(th[Tm1:], delta=th[:Tm1], extra=True)
        if m is None:
            return np.zeros((C.shape[0], th.size))
        return np.sqrt(G) * (C.T @ (m["J"] / G))

    try:
        sol = optimize.least_squares(resid, theta1, jac=jac, method="lm", xtol=1e-12, ftol=1e-12, max_nfev=200)
    except Exception:  # noqa: BLE001 - numerical failure is recorded, not raised
        return out
    th = sol.x
    m = mod.moments(th[Tm1:], delta=th[:Tm1], extra=True)
    if m is None:
        return out
    gbar = mod.cluster_sum(m["gi"]).sum(1) / G
    Jstat = float(G * gbar @ Wm @ gbar)
    Gam = m["J"] / G
    try:
        V = np.linalg.inv(Gam.T @ Wm @ Gam) / G
    except np.linalg.LinAlgError:
        return out
    Vb = V[Tm1:, Tm1:]
    se = np.sqrt(np.clip(np.diag(Vb), 0, None))
    out.update(ok=bool(sol.success and np.isfinite(Jstat)), J=Jstat, J_p=float(stats.chi2.sf(Jstat, n_extra)),
               est={n: float(v) for n, v in zip(names, th[Tm1:])}, se={n: float(v) for n, v in zip(names, se)},
               cond_Omega=float(np.linalg.cond(Om1)))
    return out


# ------------------------------------------------------------------------------------------------
# Specification builders
# ------------------------------------------------------------------------------------------------
def spec_m0(panel: dict, phi_bar: float = 0.5) -> QDSpec:
    return QDSpec(lin={"tau": panel["D"], "kappa": panel["Db"]}, phi_free=True, phi_bar=phi_bar)


def spec_m0prime(panel: dict, phi_bar: float = 0.0) -> QDSpec:
    """kappa = phi = 0 imposed (Theorem 3(d)); scale weight uses the local pressure only."""
    return QDSpec(lin={"tau": panel["D"]}, phi_free=False, phi_fixed=0.0, phi_bar=phi_bar)


def event_bins(panel: dict, E: int = 6, min_cells: int = 1):
    """Event-time indicators 1{t - A_i = e} for e = 0..E-1 and 1{t - A_i >= E} (absorbing policy).
    Bins with fewer than min_cells treated cells in weeks 0..T-2 are dropped (not identified)."""
    D, ev = panel["D"], panel["ev"]
    T = D.shape[0]
    et = np.where(D == 1, np.arange(T)[:, None] - ev[None, :], -1)
    bins = {}
    for e in range(E + 1):
        B = ((np.minimum(et, E) == e) & (D == 1)).astype(float)
        if B[:-1].sum() >= min_cells and B.sum() > 0:
            bins["tau_e%d" % e] = B
    return bins


def spec_m0_et(panel: dict, E: int = 6, phi_bar: float = 0.5) -> QDSpec:
    lin = event_bins(panel, E)
    lin["kappa"] = panel["Db"]
    return QDSpec(lin=lin, phi_free=True, phi_bar=phi_bar, G=int(panel["G"]))


def pre_exposed(panel: dict) -> np.ndarray:
    """Units already under an absorbing policy in the pre-sample weeks (any Dall[0..S-1] = 1, given as panel["Dpre"]).
    Their ev = 0 only because panel["ev"] is computed inside the window, so without a flag they are
    pooled with units that switch on in the first window week."""
    Dpre = np.asarray(panel.get("Dpre", np.zeros((0, panel["D"].shape[1]))), float)
    return Dpre.any(0)


def attach_presample(panel: dict, Dpre: np.ndarray) -> dict:
    """Return panel with the pre-sample policy path Dpre (S x N) attached as panel["Dpre"]."""
    out = dict(panel)
    out["Dpre"] = np.asarray(Dpre, float)
    return out


def spec_m0_et2(panel: dict, E: int = 6, phi_bar: float = 0.5) -> QDSpec:
    """M0-ET excluding the pre-exposed units (those already under the policy in the pre-sample weeks).

    Requires panel["Dpre"] (S x N pre-sample path, see attach_presample). Raises KeyError if any
    pre-exposed unit is in the panel: their ev = 0 only because ev is computed inside the window and their
    regressors do not move while treated, so inside this panel they are collinear with an intercept-like
    block and no event-time coefficient is identified; the sample is cut instead of pooled or split.

    Units with any Dall[:S] = 1 leave the sample (application.make_panel(..., drop_presample=True)); the
    neighbour share Db = W D is formed in the full sample first, so the removed units still enter their
    neighbours' Db. Bins as in spec_m0_et: tau_e0 = switches on in window week 0, tau_e1..E-1, open-ended
    tau_eE, plus kappa. The spec records its own cluster count (G) for the Wald test."""
    if "Dpre" not in panel:
        raise KeyError("spec_m0_et2 needs panel['Dpre'] (pre-sample policy path); use attach_presample(panel, Dpre)")
    if pre_exposed(panel).any():
        raise ValueError("spec_m0_et2: panel still contains pre-exposed units; build the panel without them "
                         "(application.make_panel(..., drop_presample=True))")
    lin = event_bins(panel, E)
    lin["kappa"] = panel["Db"]
    return QDSpec(lin=lin, phi_free=True, phi_bar=phi_bar, G=int(panel["G"]))


def spec_m0_oi(panel: dict, phi_bar: float = 0.5) -> QDSpec:
    """M0 plus two predetermined instruments: lagged import share N_{t-1}/Lambda_{t-1}(phi_bar)
    (zero at the first quasi-difference) and lagged growth log((Y_{t-1}+0.5)/(Y_{t-2}+0.5))."""
    sp = spec_m0(panel, phi_bar)
    L, Nb, om = panel["L"], panel["Nb"], panel["omega"]
    lam_bar = L + phi_bar * Nb + om
    share = Nb / lam_bar
    T, N = panel["Y"].shape
    z1 = np.zeros((T - 1, N))
    z1[1:] = share[:-2]
    Yfull = np.vstack([panel["Ypre"], panel["Y"]])          # rows: -S..T-1
    S = panel["Ypre"].shape[0]
    # for quasi-difference row t (0..T-2): Y_{t-1} = Yfull[S+t-1], Y_{t-2} = Yfull[S+t-2]
    idx = np.arange(T - 1) + S
    z2 = np.log((Yfull[idx - 1] + 0.5) / (Yfull[idx - 2] + 0.5))
    sp.extra_inst = [z1, z2]
    sp.extra_names = ["lag_import_share", "lag_growth"]
    return sp


# ------------------------------------------------------------------------------------------------
# M0-G: global root search over the parameter space phi in [0, PHI_MAX] (A5/A6) with root selection
# ------------------------------------------------------------------------------------------------
PHI_GRID = (0.0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0)
ROOT_DEDUP_ATOL = 1e-8
ROOT_DEDUP_RTOL = 1e-7


def _append_unique_root(candidates, beta) -> bool:
    """Append an interior solution unless its full parameter vector matches an existing root."""
    beta = np.asarray(beta, float)
    for _, existing in candidates:
        if np.allclose(existing, beta, rtol=ROOT_DEDUP_RTOL, atol=ROOT_DEDUP_ATOL):
            return False
    candidates.append(("interior", beta))
    return True


def _solve_given_phi(mod: QDModel, phi: float, x0):
    """Solve the K linear-index equations with phi fixed; returns (b, ok, F_phi/fscale)."""
    K = mod.K

    def f(b):
        F, J = mod.profiled_eq(np.r_[b, phi])
        return F[:K], J[:K, :K]
    sol = optimize.root(f, np.asarray(x0, float), jac=True, method="hybr", options={"xtol": 1e-10, "maxfev": 400})
    F, _ = mod.profiled_eq(np.r_[sol.x, phi])
    ok = bool(np.all(np.isfinite(sol.x)) and np.max(np.abs(F[:K])) < 1e-7 * mod.fscale)
    return sol.x, ok, float(F[K] / mod.fscale)


def selection_criterion(mod_oi: QDModel, beta) -> float:
    """Q(beta) = G gbar' Omega(beta)^{-1} gbar over all moments (week, parameter and the two predetermined
    extra instruments of M0-OI), week effects profiled; used only to choose among candidate solutions."""
    m = mod_oi.moments(beta, extra=True)
    if m is None:
        return float("inf")
    G = mod_oi.G
    gg = mod_oi.cluster_sum(m["gi"])
    gbar = gg.sum(1) / G
    Om = gg @ gg.T / G
    try:
        return float(G * gbar @ np.linalg.solve(Om, gbar))
    except np.linalg.LinAlgError:
        return float("inf")


SEL_PHI_REFS = (0.05, 2.0)


def spec_selection_inv(panel: dict, phi_bar: float = 0.5, refs=SEL_PHI_REFS) -> QDSpec:
    """Selection spec of M0-GI (post-hoc addition, run-sim-application; experiment_design.md deviation 7):
    M0 moments plus the phi-instrument evaluated at other reference values, N_t/Lambda_t(phi_r) for phi_r in
    refs. Under A1-A4 these are valid at the true root whatever phi_r, while the spurious root of the M0
    equations at phi0 = 0 sits near phi_bar and moves with it, so they discriminate between the candidates."""
    sp = spec_m0(panel, phi_bar)
    L, Nb, om = panel["L"], panel["Nb"], panel["omega"]
    sp.extra_inst = [(Nb / (L + pr * Nb + om))[:-1] for pr in refs]
    sp.extra_names = ["phiinst_%g" % pr for pr in refs]
    return sp


def fit_global(panel: dict, phi_bar: float = 0.5, grid=PHI_GRID, selection: str = "oi") -> dict:
    """M0-G. (1) Profile the phi-equation over a grid of phi in [0, max(grid)], solving (tau, kappa) at each
    grid value; (2) refine every sign change with the full three-equation solver; (3) add the boundary
    candidate phi = 0 (constrained estimate); (4) choose the candidate with the smallest selection criterion
    Q (all moments incl. the M0-OI extra instruments; selection="inv" -> M0-GI, which uses the
    phi-instruments at other reference values instead, see spec_selection_inv). Interior choice: usual M0 inference. Boundary choice:
    (tau, kappa) inference with phi fixed at 0 (conditional on the boundary; flagged at_boundary)."""
    sp = spec_m0(panel, phi_bar)
    mod = QDModel(panel, sp)
    if selection == "oi":        # M0-G (pre-specified): extras = M0-OI's lagged import share and lagged growth
        mod_oi = QDModel(panel, spec_m0_oi(panel, phi_bar))
    elif selection == "inv":     # M0-GI (post-hoc): extras = phi-instruments at other reference values
        mod_oi = QDModel(panel, spec_selection_inv(panel, phi_bar))
    else:
        raise ValueError("selection must be 'oi' or 'inv'")
    names = sp.names
    prof = []
    x = np.zeros(mod.K)
    for ph in grid:
        b, ok, Fp = _solve_given_phi(mod, ph, x)
        if ok:
            x = b
        prof.append((ph, b, ok, Fp))
    cands = []
    for (p0, b0, ok0, F0), (p1, b1, ok1, F1) in zip(prof[:-1], prof[1:]):
        if ok0 and ok1 and np.sign(F0) != np.sign(F1):
            w = abs(F0) / (abs(F0) + abs(F1))
            xs = np.r_[(1 - w) * b0 + w * b1, (1 - w) * p0 + w * p1]
            xr, okr, res = mod.solve(xs)
            if okr and xr[-1] >= 0:
                _append_unique_root(cands, xr)
    b0, ok0, F0 = prof[0][1], prof[0][2], prof[0][3]
    if ok0:
        cands.append(("boundary", np.r_[b0, 0.0]))
    out = dict(names=names, ok=False, n_roots=sum(c[0] == "interior" for c in cands),
               roots=[float(c[1][-1]) for c in cands if c[0] == "interior"],
               profile=[(float(p), float(F)) for p, _, okp, F in prof if okp])
    if not cands:
        out.update(est={n: float("nan") for n in names}, se={n: float("nan") for n in names}, Vb=None,
                   reason="no candidate")
        return out
    Q = [selection_criterion(mod_oi, c[1]) for c in cands]
    k = int(np.argmin(Q))
    kind, beta = cands[k]
    out.update(Q=[float(q) for q in Q], chosen=kind, at_boundary=kind == "boundary", beta=beta)
    if kind == "interior":
        inf = mod.inference(beta)
        Vb = None if inf is None else inf["Vb"]
        se = np.sqrt(np.diag(Vb)) if Vb is not None and np.all(np.diag(Vb) >= 0) else np.full(len(names), np.nan)
    else:
        mod_fix = QDModel(panel, QDSpec(lin=dict(sp.lin), phi_free=False, phi_fixed=0.0, phi_bar=phi_bar))
        inf = mod_fix.inference(beta[:-1])
        Vs = None if inf is None else inf["Vb"]
        Vb = None
        if Vs is not None:
            Vb = np.zeros((len(names), len(names)))
            Vb[:-1, :-1] = Vs
        se = np.r_[np.sqrt(np.diag(Vs)), np.nan] if Vs is not None else np.full(len(names), np.nan)
    out.update(ok=bool(np.all(np.isfinite(beta)) and Vb is not None),
               est={n: float(v) for n, v in zip(names, beta)}, se={n: float(v) for n, v in zip(names, se)}, Vb=Vb)
    return out
