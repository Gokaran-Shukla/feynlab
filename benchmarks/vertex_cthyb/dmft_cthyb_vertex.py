"""Square-lattice DMFT with segment CT-HYB (continuous bath), then G2 / vertex on a box of 2*n_nu fermionic and
2*n_om+1 bosonic frequencies; checks the Schwinger-Dyson equation. Writes /tmp/g2/ctqmc_b*_U*_n*.pkl.
Usage: python dmft_cthyb_vertex.py BETA U N_NU N_OM N_DMFT_ITER N_MEAS_G2   e.g. 10 1.0 16 15 12 2000000"""
import numpy as np, time, sys, pickle
from feynlab.core.matsubara import MatsubaraGrid
from feynlab.qmc.cthyb import CTHYBSolver
from feynlab.core.lattice import square_lattice_eps
from feynlab.vertex import vertices, eom_sigma
beta,U=float(sys.argv[1]),float(sys.argv[2]); nn,no=int(sys.argv[3]),int(sys.argv[4])
eps=square_lattice_eps(16,0.25).ravel(); mu=U/2
g=MatsubaraGrid(beta,1024); iw=g.iw
nmeas_iw=64
D=0.25*0.25*4/iw  # initial guess ~ t^2 z / iw
S=CTHYBSolver(g,U,n_iw_meas=nmeas_iw)
sig_full=np.full(len(iw),U/2+0j)
t0=time.time()
for it in range(int(sys.argv[5])):
    r=S.solve(D,mu,n_meas=300000,meas_every=10,n_chains=2,seed=11+it)
    sig=r['Sigma']
    # high-frequency tail: U/2 + U^2 n(1-n/2)... half filling: U/2 + U^2/4 / iw
    sig_full=U/2+U**2/4/iw; sig_full[:nmeas_iw-8]=sig[:nmeas_iw-8]
    Gl=np.mean(1/(iw[:,None]+mu-eps[None,:]-sig_full[:,None]),axis=1)
    Dn=iw+mu-sig_full-1/Gl
    diff=np.abs(Dn-D)[:20].max(); D=0.5*Dn+0.5*D
    print(it,"diff",diff,"docc",r['docc'],"ImSig0",sig[0].imag,time.time()-t0,flush=True)
r=S.solve(D,mu,n_meas=int(sys.argv[6]),meas_every=5,n_chains=2,seed=999,g2=(nn,no))
print("G2 run",time.time()-t0,flush=True)
gs=r['gs']; v=vertices(gs,beta)
Seom=eom_sigma(v['F_ud'],gs,beta,U,1.0)
nus=gs['nu']; i0=nn
print("SDE check: Sigma_EOM vs Sigma_QMC (first 4 positive nu)")
for j in range(4): print(nus[i0+j], Seom[i0+j], r['Sigma'][j])
gs.pop("G_of"); r["gs"]=gs; __import__("os").makedirs("/tmp/g2", exist_ok=True)
pickle.dump(dict(r=r,D=D,sig_full=sig_full,v={k:v[k] for k in v},beta=beta,U=U),open(f"/tmp/g2/ctqmc_b{beta}_U{U}_n{nn}.pkl","wb"))
