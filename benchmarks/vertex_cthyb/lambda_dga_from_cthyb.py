"""lambda-DGA (16x16 square lattice) from a CT-HYB vertex pickle written by dmft_cthyb_vertex.py.
Usage: python lambda_dga_from_cthyb.py /tmp/g2/ctqmc_b10.0_U1.0_n16.pkl"""
import numpy as np, pickle, sys
from feynlab.vertex import lambda_dga
from feynlab.core.lattice import square_lattice_eps
d=pickle.load(open(sys.argv[1],"rb")); beta,U=d['beta'],d['U']; N=16
eps=square_lattice_eps(N,0.25); v=d['v']; gs=d['r']['gs']
sf=d['sig_full']; nfull=len(sf); allp=(2*np.arange(nfull)+1)*np.pi/beta
def sd(x):
    x=np.asarray(x); n=np.round((np.abs(x)*beta/np.pi-1)/2).astype(int)
    out=sf[np.minimum(n,nfull-1)].copy(); out[x<0]=np.conj(out[x<0]); return out
Gl=lambda x: np.mean(1/(1j*x[None,None,:]+U/2-eps[...,None]-sd(x)[None,None,:]),axis=(0,1))
gs2=dict(nu=gs['nu'],om=gs['om'],G_of=Gl)
w0=len(gs['om'])//2; nn=len(gs['nu'])//2
for kw in (dict(lambda_mode=None,n_tail=256),dict(n_tail=256),dict(n_tail=256,docc=d['r']['docc'])):
    r=lambda_dga(eps,U/2,sd,v,gs2,beta,U,1.0,**kw)
    iq=r['qs'].index((N//2,N//2))
    print({k:val for k,val in kw.items() if k!='n_tail'},"lam %.4f target %.4f chiQ %.3f ImS(pi,0) %.4f ImS(pi/2,pi/2) %.4f"%(r['lam_m'],r['sum_target'],r['chi_m_lambda'][iq,w0].real,r['Sigma'][N//2,0,nn].imag,r['Sigma'][N//4,N//4,nn].imag))
