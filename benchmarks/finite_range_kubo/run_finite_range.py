"""Finite-range (Gaussian) disorder: Kubo-Streda NCA (FiniteRangeHall) vs Boltzmann (analytic sigma_xx and
SemiclassicalHall side jump + anomalous distribution + Gaussian skew) for the massive Dirac model.
Usage: python run_finite_range.py N_I_V0SQ N K d1 d2 ...   e.g. 0.01 1024 2.0 0.0 0.5 1.0 2.0 3.0 (about 1 min per d)"""
import numpy as np, time, sys
from feynlab.transport.models import MassiveDirac
from feynlab.transport.boltzmann import SemiclassicalHall, dirac_finite_range_sigma_xx
from feynlab.transport.finite_range import FiniteRangeHall
dm=MassiveDirac(1.0,0.5); mu=1.0; nV2=float(sys.argv[1]); N=int(sys.argv[2]); K=float(sys.argv[3])
for rd in [float(x) for x in sys.argv[4:]]:
    t0=time.time()
    fr=FiniteRangeHall(dm,mu,n_i=1.0,V0=np.sqrt(nV2),d=rd,K=K,N=N)
    r=fr.decompose(intrinsic_value=dm.ref_intrinsic(mu))
    sxx=fr.sigma_xx()
    sc=SemiclassicalHall(dm,mu,1,n_i=1.0,V0=np.sqrt(nV2),range_d=rd,n_theta=256).compute()
    ext_sc=sc['side_jump']+sc['anomalous_distribution']+sc['skew_gaussian']
    print(f"d={rd}: sxx Kubo {sxx:.4f} Boltz {dirac_finite_range_sigma_xx(1,0.5,mu,1.0,np.sqrt(nV2),rd):.4f} | ext Kubo {r['extrinsic_gaussian']:.5f} semicl {ext_sc:.5f} (sj {sc['side_jump']:.4f} ad {sc['anomalous_distribution']:.4f} gs {sc['skew_gaussian']:.4f}) Gam {r['Gamma_broadening']:.4f} t {time.time()-t0:.0f}s gmres {fr.gmres_info}", flush=True)
