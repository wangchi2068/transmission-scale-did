"""Scenario table (theory/simulation_spec.md: 16 core + 4 null + 5 misspecification = 25 scenarios).

Scenario ids are fixed (they enter the seed: seed = BASE_SEED + id * 1e5 + rep); RUN_ORDER gives the
priority in which the runner schedules them (core feedback with spillover first, then the rest of the core
block, then the null block, then the misspecification block).

Pilot-calibration deviation (documented in experiment_design.md): the unit-effect components of the spec
(a_g, b_i ~ N(0, 0.2^2), 0.3 pi_g) are multiplied by ALPHA_SCALE = 0.25, because with the spec values and
no susceptible depletion, persistent log-growth heterogeneity produced county-week counts of up to 1e12
per unit population scale at T = 26 (pilot table in logs/pilot_alpha_tails.txt).
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, replace
from typing import Dict, List

import numpy as np

from dgp import DGPParams

BASE_SEED = 20261112
PSI_SEED_OFFSET = 900_000_000          # independent stream for the population psi truth
ALPHA_SCALE = 0.25
TAU0 = float(np.log(0.8))
HERE = os.path.dirname(os.path.abspath(__file__))
CALIB_PATH = os.path.join(os.path.dirname(HERE), "calibration.json")


@dataclass(frozen=True)
class Scenario:
    sid: int
    name: str
    block: str              # core | null | misspec
    params: DGPParams
    w_est: tuple = (0.85, 0.15)
    psi_pairs: int = 5000

    def to_dict(self):
        return dict(sid=self.sid, name=self.name, block=self.block, params=self.params.to_dict(),
                    w_est=list(self.w_est))


def load_calibration(path: str = CALIB_PATH) -> dict:
    if not os.path.exists(path):
        raise FileNotFoundError("calibration.json missing: run `python src/calibrate.py` first")
    with open(path) as fh:
        return json.load(fh)


def _base(N, T, adopt, phi, kappa, tau, k, gbar) -> DGPParams:
    return DGPParams(N=N, T=T, adopt=adopt, tau=tau, kappa=kappa, phi=phi, k=k, gbar=gbar,
                     alpha_state_sd=0.2 * ALPHA_SCALE, alpha_unit_sd=0.2 * ALPHA_SCALE,
                     alpha_pi=0.3 * ALPHA_SCALE)


def build_scenarios(calib: dict = None) -> Dict[int, Scenario]:
    calib = load_calibration() if calib is None else calib
    gb = calib["gbar"]                                   # keys "T{T}_{adopt}"
    out: Dict[int, Scenario] = {}
    sid = 1
    # core block: ids 1..16
    for adopt in ("feedback", "random"):
        for spill in (True, False):
            for N in (500, 2000):
                for T in (12, 26):
                    phi, kap = (0.4, -0.1) if spill else (0.0, 0.0)
                    nm = "core_N%d_T%d_%s_%s" % (N, T, "fb" if adopt == "feedback" else "rd",
                                                 "spill" if spill else "nospill")
                    out[sid] = Scenario(sid, nm, "core",
                                        _base(N, T, adopt, phi, kap, TAU0, 10.0, gb["T%d_%s" % (T, adopt)]))
                    sid += 1
    # null block: ids 17..20
    for adopt in ("feedback", "random"):
        for k in (10.0, np.inf):
            nm = "null_N500_T12_%s_k%s" % ("fb" if adopt == "feedback" else "rd", "10" if k == 10 else "inf")
            out[sid] = Scenario(sid, nm, "null", _base(500, 12, adopt, 0.4, -0.1, 0.0, k, gb["T12_%s" % adopt]))
            sid += 1
    # misspecification block: ids 21..25 (base N = 2000, T = 26, feedback, spillover)
    base = _base(2000, 26, "feedback", 0.4, -0.1, TAU0, 10.0, gb["T26_feedback"])
    out[21] = Scenario(21, "X1_wrong_generation_interval", "misspec", replace(base, w_true=(0.7, 0.3)))
    out[22] = Scenario(22, "X2_reporting_drops_after_adoption", "misspec", replace(base, report_p=0.5, report_drop=0.85))
    out[23] = Scenario(23, "X3_heterogeneous_effects", "misspec", replace(base, tau_sd=0.1))
    out[24] = Scenario(24, "X4_susceptible_depletion", "misspec", replace(base, depletion_pop=calib["depletion_pop"]))
    out[25] = Scenario(25, "X5_dynamic_effects", "misspec",
                       replace(base, dyn=(-0.1, -0.15, -0.2, -0.25, -0.3, -0.3)))
    return out


def run_order(scen: Dict[int, Scenario]) -> List[int]:
    """Priority: core feedback+spillover, core feedback no-spill, core random+spill, core random no-spill,
    null, misspecification. Within a group, the cheaper N = 500 scenarios first."""
    def key(s: Scenario):
        if s.block == "core":
            g = (0 if s.params.adopt == "feedback" else 2) + (0 if s.params.phi > 0 else 1)
        elif s.block == "null":
            g = 4
        else:
            g = 5
        return (g, s.params.N, s.params.T, s.sid)
    return [s.sid for s in sorted(scen.values(), key=key)]


def seed_for(sid: int, rep: int, base: int = BASE_SEED) -> int:
    return base + sid * 100_000 + rep


def default_workers() -> int:
    """Worker count: env QDG_WORKERS, else the usable cores (CPU affinity, capped by a cgroup v2 CPU quota)."""
    env = os.environ.get("QDG_WORKERS")
    if env:
        return max(1, int(env))
    try:
        n = len(os.sched_getaffinity(0))
    except AttributeError:
        n = os.cpu_count() or 1
    try:
        with open("/sys/fs/cgroup/cpu.max") as fh:
            q, per = fh.read().split()[:2]
        if q != "max":
            n = min(n, max(1, int(float(q) / float(per))))
    except (OSError, ValueError):
        pass
    return max(1, n)
