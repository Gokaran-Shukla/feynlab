"""Two-particle Green function of a 2-level discrete bath: segment CT-HYB (removal estimator) vs Lehmann ED.
Usage: python ed_vs_cthyb_g2.py N_MEAS   (ED about 5 min, CT-HYB 16 s for N_MEAS = 400000)"""
import numpy as np, time, sys
from feynlab.core.matsubara import MatsubaraGrid
from feynlab.qmc.cthyb import CTHYBSolver
from feynlab.vertex import ImpurityED, vertices
beta,U,mu=5.0,2.0,1.0
eb=np.array([-1.0,1.0]); V=np.array([0.6,0.6])
nn,no=4,3
t0=time.time()
ed=ImpurityED(-mu,U,eb,V); gs_ed=ed.generalized_susceptibilities(beta,n_nu=nn,n_om=no); print("ED",time.time()-t0,flush=True)
g=MatsubaraGrid(beta,1024)
D=np.sum(V[None,:]**2/(g.iw[:,None]-eb[None,:]),axis=1)
t0=time.time()
r=CTHYBSolver(g,U,n_iw_meas=64).solve(D,mu,n_meas=int(sys.argv[1]),meas_every=10,n_chains=4,g2=(nn,no))
print("QMC",time.time()-t0, "docc",r['docc'],ed.double_occupancy(beta),flush=True)
q=r['gs']
for key in ("chi_uu","chi_ud"):
    d=q[key]-gs_ed[key]; e=np.abs(q[key+"_err"])+1e-12
    print(key,"max|diff|",np.abs(d).max(),"max|chi|",np.abs(gs_ed[key]).max(),"median |diff|/err",np.median(np.abs(d)/e),"95%",np.percentile(np.abs(d)/e,95))
cm=lambda gs: np.real((gs['chi_uu']-gs['chi_ud']).sum(axis=(1,2)))/beta**2
print("chi_m(w) ED ",np.round(cm(gs_ed),5)); print("chi_m(w) QMC",np.round(cm(q),5))
