import sys, numpy as np
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', '..', 'src'))
import simulate as S, dgp, qdgmm, scenarios as sc
from scipy import optimize
sid = int(sys.argv[1]); scn = S.scenario_table()[sid]; p = scn.params; geo = S.geography(p.N)
for r in range(4):
    seed = sc.seed_for(sid, r); sim = dgp.simulate(p, geo, [seed,0],[seed,1])
    pan = dgp.build_panel(sim['Yall'], sim['Dall'], geo.W, sim['omega'], geo.cl)
    out = []
    for pb in (0.1, 0.25, 0.5, 1.0, 2.0):
        mod = qdgmm.QDModel(pan, qdgmm.spec_m0(pan, phi_bar=pb))
        x, ok, res = mod.solve([0, 0, pb])
        x2, ok2, _ = mod.solve([p.tau, p.kappa, max(p.phi, 0.01)])
        out.append((pb, np.round(x, 3), ok, np.round(x2, 3), ok2))
    print('rep', r)
    for o in out: print('   ', o)
    # profile of phi moment at phi_bar=0.5
    mod = qdgmm.QDModel(pan, qdgmm.spec_m0(pan, 0.5))
    prof = []
    for ph in (0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5):
        f = lambda b: mod.profiled_eq(np.r_[b, ph])[0][:2]
        sol = optimize.root(f, [p.tau, p.kappa], method='hybr')
        F = mod.profiled_eq(np.r_[sol.x, ph])[0]
        prof.append((ph, round(F[2]/mod.fscale, 5), np.round(sol.x, 2), bool(np.max(np.abs(F[:2]))<1e-7*mod.fscale)))
    print('   profile', prof)
