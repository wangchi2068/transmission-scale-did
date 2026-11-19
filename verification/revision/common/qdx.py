"""Grouped extension of the QD-GMM estimator used by the revision analyses (review Comment 1, 3, 4, 5).

QDX generalises qdgmm.QDModel in two directions:
  * week effects by group:   delta_{t,r}, r = wgroup[i] (e.g. Census division x week); profiled in closed form
                              within each group, Ed_{t,r} = sum_{i in r} s Yq / sum_{i in r} s Y_t;
  * group-specific phi:       Lambda_it = L_it + phi_{pg(i)} N_it + omega_i, one phi-instrument per group
                              (s * phi_bar * N / Lambda_bar restricted to the group).
With wgroup = phigroup = None it reproduces qdgmm.QDModel(spec_m0) exactly (checked in validate.py).

Over-identification is measured in score (LM / C(alpha)) form because the grouped week moments outnumber the
49 clusters, so the full cluster covariance Omega is singular and the Hansen two-step J cannot be formed:
    psi_g = g_{e,g} - Gamma_e Gamma_p^{-1} g_{p,g}      (extra moments net of their projection on the
                                                          just-identified moments, evaluated at the JI root)
    J_LM  = G psibar' S_psi^{-1} psibar,  S_psi = sum_g psi_g psi_g' / G   (asymptotically chi2(#extras))
At the just-identified root gbar_p = 0, so psibar = gbar_e.  "Q*" below is the same statistic used as a
root-selection criterion (evaluated at each candidate root).
"""
from __future__ import annotations

import numpy as np
from scipy import optimize, stats


class QDX:
    def __init__(self, panel, lin, wgroup=None, phigroup=None, phi_bar=0.5, extras=(), extra_names=(),
                 phi_free=True, phi_fixed=None):
        self.p = panel
        Y = panel["Y"]
        self.Y = Y
        self.T, self.N = Y.shape
        self.L, self.Nb, self.om = panel["L"], panel["Nb"], panel["omega"]
        self.cl, self.G = np.asarray(panel["cl"]), int(panel["G"])
        self.names_lin = list(lin)
        self.X = [np.asarray(v, float) for v in lin.values()]
        self.K = len(self.X)
        self.wg = np.zeros(self.N, int) if wgroup is None else np.asarray(wgroup, int)
        self.R = int(self.wg.max()) + 1
        self.pg = np.zeros(self.N, int) if phigroup is None else np.asarray(phigroup, int)
        self.P = int(self.pg.max()) + 1
        self.phi_free = phi_free
        self.phi_fixed = np.zeros(self.P) if phi_fixed is None else np.broadcast_to(np.asarray(phi_fixed, float), (self.P,)).copy()
        lam_bar = self.L + phi_bar * self.Nb + self.om
        self.s = 1.0 / (lam_bar[:-1] + 1.0)
        Z = [self.s * x[:-1] for x in self.X]
        if phi_free:
            base = self.s * phi_bar * self.Nb[:-1] / lam_bar[:-1]
            for p in range(self.P):
                Z.append(base * (self.pg == p)[None, :])
        self.Zpar = Z
        self.Zext = [self.s * np.asarray(z, float) for z in extras]
        self.extra_names = list(extra_names)
        self.dX = [x[:-1] - x[1:] for x in self.X]
        self.np_ = self.K + (self.P if phi_free else 0)
        self.names = self.names_lin + ((["phi"] if self.P == 1 else ["phi_%d" % p for p in range(self.P)]) if phi_free else [])
        self.fscale = max(float(np.sum(self.s * Y[:-1])) / self.N, 1e-12)
        self.Tm1 = self.T - 1
        self.nW = self.Tm1 * self.R
        # one-hot of week groups for fast group sums: (N, R)
        self.Oh = np.zeros((self.N, self.R))
        self.Oh[np.arange(self.N), self.wg] = 1.0

    # --------------------------------------------------------------------------------------------
    def _split(self, beta):
        beta = np.asarray(beta, float)
        b = beta[:self.K]
        phi = beta[self.K:] if self.phi_free else self.phi_fixed
        return b, phi

    def gsum(self, A):
        """(T-1, N) -> (T-1, R) sums within week groups."""
        return A @ self.Oh

    def _core(self, beta, delta=None):
        b, phi = self._split(beta)
        if not np.all(np.isfinite(beta)):
            return None
        lam = self.L + phi[self.pg][None, :] * self.Nb + self.om
        if np.any(lam <= 0):
            return None
        eta = np.zeros_like(self.Y)
        for bk, x in zip(b, self.X):
            if bk != 0.0:
                eta = eta + bk * x
        q = np.exp(np.clip(eta[:-1] - eta[1:], -50, 50)) * lam[:-1] / lam[1:]
        Yq = self.Y[1:] * q
        if delta is None:
            den = self.gsum(self.s * self.Y[:-1])
            num = self.gsum(self.s * Yq)
            if np.any(den <= 0) or np.any(num <= 0):
                return None
            Ed = num / den                                 # (T-1, R)
        else:
            Ed = np.exp(np.asarray(delta, float).reshape(self.Tm1, self.R))
        Yr = Yq / Ed[:, self.wg]
        v = Yr - self.Y[:-1]
        dx = list(self.dX)
        if self.phi_free:
            r = self.Nb / lam
            rr = r[:-1] - r[1:]
            for p in range(self.P):
                dx.append(rr * (self.pg == p)[None, :])
        return dict(v=v, Yr=Yr, dx=dx, Ed=Ed, lam=lam)

    def _blocks(self, c, Zlist):
        sYr = self.s * c["Yr"]
        a = self.gsum(sYr)                                                      # (T-1, R)
        b = np.stack([self.gsum(sYr * d) for d in c["dx"]], -1)                 # (T-1, R, p)
        cm = np.stack([self.gsum(z * c["Yr"]) for z in Zlist], 0)               # (m, T-1, R)
        H = np.array([[(z * c["Yr"] * d).sum() for d in c["dx"]] for z in Zlist])  # (m, p)
        return a, b, cm, H

    def profiled_eq(self, beta):
        c = self._core(beta)
        if c is None:
            return np.full(self.np_, 1e6 * self.fscale), np.eye(self.np_)
        F = np.array([(z * c["v"]).sum() for z in self.Zpar]) / self.N
        a, b, cm, H = self._blocks(c, self.Zpar)
        Jp = (H - np.einsum("mtr,trj->mj", cm, b / a[:, :, None])) / self.N
        return F, Jp

    def solve(self, x0, xtol=1e-10, maxfev=600):
        sol = optimize.root(self.profiled_eq, np.asarray(x0, float), jac=True, method="hybr",
                            options={"xtol": xtol, "maxfev": maxfev})
        F, _ = self.profiled_eq(sol.x)
        ok = bool(np.all(np.isfinite(sol.x)) and np.max(np.abs(F)) < 1e-7 * self.fscale)
        return sol.x, ok, float(np.max(np.abs(F)) / self.fscale)

    def fit(self, starts=None):
        if starts is None:
            ph = [0.5] * (self.np_ - self.K)
            starts = [[0.0] * self.K + ph, [-0.2] * self.K + [0.2] * len(ph), [0.0] * self.K + [1.0] * len(ph)]
        best = None
        for x0 in starts:
            x, ok, res = self.solve(x0)
            if best is None or (ok and not best[1]) or (ok == best[1] and res < best[2]):
                best = (x, ok, res)
            if ok:
                break
        return best

    # --------------------------------------------------------------------------------------------
    def moments(self, beta, extra=True):
        """Cluster sums gg (q x G) and Jacobian sums J (q x P_full); rows = [week-group, params, extras]."""
        c = self._core(beta)
        if c is None:
            return None
        Zl = self.Zpar + (self.Zext if extra else [])
        sv = self.s * c["v"]
        # week-group unit contributions summed to clusters: row index (t, r) -> t*R + r
        G = self.G
        gw = np.zeros((self.Tm1, self.R, G))
        for t in range(self.Tm1):
            np.add.at(gw[t], (self.wg, self.cl), sv[t])
        gg_w = gw.reshape(self.nW, G)
        gg_z = np.vstack([np.bincount(self.cl, (z * c["v"]).sum(0), G) for z in Zl])
        gg = np.vstack([gg_w, gg_z])
        a, b, cm, H = self._blocks(c, Zl)
        P = self.nW + self.np_
        Q = self.nW + len(Zl)
        J = np.zeros((Q, P))
        J[np.arange(self.nW), np.arange(self.nW)] = -a.reshape(-1)
        J[:self.nW, self.nW:] = b.reshape(self.nW, -1)
        J[self.nW:, :self.nW] = -cm.reshape(len(Zl), -1)
        J[self.nW:, self.nW:] = H
        return dict(gg=gg, J=J, c=c)

    def inference(self, beta, m=None):
        m = self.moments(beta, extra=False) if m is None else m
        if m is None:
            return None
        G = self.G
        P = self.nW + self.np_
        Gam = m["J"][:P] / G
        gg = m["gg"][:P]
        try:
            Gi = np.linalg.inv(Gam)
        except np.linalg.LinAlgError:
            return None
        infl = -Gi @ gg                                   # (P, G) cluster influence terms (sum = 0 at root)
        V = infl @ infl.T / G / G * G / (G - 1)
        return dict(V=V, Vb=V[self.nW:, self.nW:], infl_b=infl[self.nW:], Gam=Gam, Gi=Gi, gg=gg)

    def j_lm(self, beta, centered=False, sel=None, m=None):
        """Score-form over-identification statistic for the extra instruments (see module doc).
        sel: indices of the extra instruments to use (default all)."""
        m = self.moments(beta, extra=True) if m is None else m
        if m is None or not self.Zext:
            return dict(J_LM=float("nan"), df=len(self.Zext), p=float("nan"))
        G = self.G
        P = self.nW + self.np_
        Gam_p = m["J"][:P] / G
        Gam_e = m["J"][P:] / G
        gp, ge = m["gg"][:P], m["gg"][P:]
        if sel is not None:
            ge, Gam_e = ge[list(sel)], Gam_e[list(sel)]
        try:
            psi = ge - Gam_e @ np.linalg.solve(Gam_p, gp)
        except np.linalg.LinAlgError:
            return dict(J_LM=float("nan"), df=len(self.Zext), p=float("nan"))
        pbar = psi.sum(1) / G
        R = psi - pbar[:, None] if centered else psi
        S = R @ R.T / G
        try:
            J = float(G * pbar @ np.linalg.solve(S, pbar))
        except np.linalg.LinAlgError:
            J = float("nan")
        df = ge.shape[0]
        return dict(J_LM=J, df=df, p=float(stats.chi2.sf(J, df)) if np.isfinite(J) else float("nan"),
                    gbar_e=(ge.sum(1) / G).tolist())

    def q_full(self, beta):
        """Paper's Q = G gbar' Omega^{-1} gbar over all moments (only defined when Omega is non-singular)."""
        m = self.moments(beta, extra=True)
        if m is None:
            return float("nan")
        G = self.G
        gbar = m["gg"].sum(1) / G
        Om = m["gg"] @ m["gg"].T / G
        if Om.shape[0] >= G:
            return float("nan")
        try:
            return float(G * gbar @ np.linalg.solve(Om, gbar))
        except np.linalg.LinAlgError:
            return float("nan")

    # --------------------------------------------------------------------------------------------
    def solve_given_phi(self, phi_vec, x0):
        K = self.K

        def f(b):
            F, J = self.profiled_eq(np.r_[b, phi_vec])
            return F[:K], J[:K, :K]
        sol = optimize.root(f, np.asarray(x0, float), jac=True, method="hybr", options={"xtol": 1e-10, "maxfev": 400})
        F, _ = self.profiled_eq(np.r_[sol.x, phi_vec])
        ok = bool(np.all(np.isfinite(sol.x)) and np.max(np.abs(F[:K])) < 1e-7 * self.fscale)
        return sol.x, ok, F[K:] / self.fscale

    def global_roots(self, grid):
        """Scalar-phi only (P == 1): profile F_phi over grid, refine sign changes, add the phi=0 boundary."""
        assert self.P == 1 and self.phi_free
        prof = []
        x = np.zeros(self.K)
        for ph in grid:
            b, ok, Fp = self.solve_given_phi(np.array([ph]), x)
            if ok:
                x = b
            prof.append((ph, b, ok, float(Fp[0])))
        cands = []
        for (p0, b0, ok0, F0), (p1, b1, ok1, F1) in zip(prof[:-1], prof[1:]):
            if ok0 and ok1 and np.sign(F0) != np.sign(F1):
                w = abs(F0) / (abs(F0) + abs(F1))
                xs = np.r_[(1 - w) * b0 + w * b1, (1 - w) * p0 + w * p1]
                xr, okr, _ = self.solve(xs)
                if okr and xr[-1] >= 0 and not any(np.allclose(xr, c[1], rtol=1e-7, atol=1e-8) for c in cands):
                    cands.append(("interior", xr))
        if prof[0][2]:
            cands.append(("boundary", np.r_[prof[0][1], 0.0]))
        return cands, prof


def tsum(V, i, j):
    return float(np.sqrt(max(V[i, i] + V[j, j] + 2 * V[i, j], 0.0)))


# ------------------------------------------------------------------------------------------------
# Census geography (state FIPS -> division / region)
# ------------------------------------------------------------------------------------------------
DIVISIONS = {
    1: (9, 23, 25, 33, 44, 50), 2: (34, 36, 42), 3: (17, 18, 26, 39, 55), 4: (19, 20, 27, 29, 31, 38, 46),
    5: (10, 11, 12, 13, 24, 37, 45, 51, 54), 6: (1, 21, 28, 47), 7: (5, 22, 40, 48),
    8: (4, 8, 16, 30, 32, 35, 49, 56), 9: (6, 41, 53),
}
REGION_OF_DIV = {1: 1, 2: 1, 3: 2, 4: 2, 5: 3, 6: 3, 7: 3, 8: 4, 9: 4}   # 1 NE, 2 MW, 3 S, 4 W
REGION_NAMES = {1: "Northeast", 2: "Midwest", 3: "South", 4: "West"}
DIV_NAMES = {1: "New England", 2: "Middle Atlantic", 3: "East North Central", 4: "West North Central",
             5: "South Atlantic", 6: "East South Central", 7: "West South Central", 8: "Mountain", 9: "Pacific"}


def division_of_fips(fips):
    m = {f: d for d, fs in DIVISIONS.items() for f in fs}
    return np.array([m[int(f)] for f in np.asarray(fips)])


def groups_for_units(state_fips):
    """Return 0-based (division index, region index) arrays for unit-level state FIPS codes."""
    div = division_of_fips(state_fips)
    reg = np.array([REGION_OF_DIV[d] for d in div])
    return div - 1, reg - 1


def oi_extras(panel, phi_bar=0.5):
    """M0-OI extra instruments (raw, before multiplication by s): lagged import share and lagged growth."""
    L, Nb, om = panel["L"], panel["Nb"], panel["omega"]
    share = Nb / (L + phi_bar * Nb + om)
    T, N = panel["Y"].shape
    z1 = np.zeros((T - 1, N))
    z1[1:] = share[:-2]
    Yfull = np.vstack([panel["Ypre"], panel["Y"]])
    S = panel["Ypre"].shape[0]
    idx = np.arange(T - 1) + S
    z2 = np.log((Yfull[idx - 1] + 0.5) / (Yfull[idx - 2] + 0.5))
    return [z1, z2], ["lag_import_share", "lag_growth"]


def inv_extras(panel, refs=(0.05, 2.0)):
    """M0-GI extra instruments: phi-instrument N/Lambda(phi_r) at other reference values."""
    L, Nb, om = panel["L"], panel["Nb"], panel["omega"]
    return [(Nb / (L + pr * Nb + om))[:-1] for pr in refs], ["phiinst_%g" % pr for pr in refs]
