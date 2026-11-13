import sys, os, numpy as np, time
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', '..', 'src'))
from multiprocessing import Pool
import simulate as S, dgp, qdgmm, scenarios as sc
def one(a):
    sid, r = a
    scn = S.scenario_table()[sid]; p = scn.params; geo = S.geography(p.N)
    seed = sc.seed_for(sid, r); sim = dgp.simulate(p, geo, [seed,0],[seed,1])
    pan = dgp.build_panel(sim['Yall'], sim['Dall'], geo.W, sim['omega'], geo.cl)
    f = qdgmm.fit(pan, qdgmm.spec_m0(pan))
    t = time.time(); g = qdgmm.fit_global(pan); tg = time.time() - t
    k = ('tau','kappa','phi')
    return sid, f['ok'], [f['est'][x] for x in k], [f['se'][x] for x in k], g['ok'], [g['est'][x] for x in k], [g['se'][x] for x in k], g['n_roots'], g.get('at_boundary'), tg
if __name__ == '__main__':
    sids = [int(x) for x in sys.argv[1].split(',')]; R = int(sys.argv[2])
    with Pool(4) as P:
        out = P.map(one, [(s, r) for s in sids for r in range(R)])
    t = qdgmm.tcrit(50)
    for s in sids:
        o = [x for x in out if x[0]==s]; p = S.scenario_table()[s].params; tr = np.array([p.tau, p.kappa, p.phi])
        for lab, io, ie, isd in (('M0 ', 1, 2, 3), ('M0G', 4, 5, 6)):
            ok = np.array([x[io] for x in o]); E = np.array([x[ie] for x in o])[ok]; SE = np.array([x[isd] for x in o])[ok]
            cov = np.nanmean(np.where(np.isfinite(SE), np.abs(E - tr) <= t*SE, np.nan), 0)
            print(s, S.scenario_table()[s].name, lab, 'ok %.2f' % ok.mean(), 'bias', np.round(E.mean(0)-tr,3), 'empSD', np.round(E.std(0),3), 'cov', np.round(cov,2))
        print('    nroots dist', np.bincount([x[7] for x in o]), 'boundary share', np.mean([bool(x[8]) for x in o]), 'time M0G %.2f' % np.mean([x[9] for x in o]))
