import numpy as np
from scipy import ndimage
from common import *
f=np.load(f'{OUT}/factors.npz'); cf=np.load(f'{OUT}/cover_frac.npz')
fr={k:cf[f'f{k}'] for k in [10,20,30,40,50,60,80,90,95]}
LDD={10:.020,20:.032,30:.015,40:.150,60:.800,95:.003}
def Wmap(sealed):
    w=sum(fr[k]*v for k,v in LDD.items())+fr[50]*sealed; vf=sum(fr[k] for k in LDD)+fr[50]
    return np.clip(w/np.maximum(vf,1e-6),0.001,1)
dem,_,_=rd('dem'); z=np.where(np.isfinite(dem),dem,0)
res=z-ndimage.uniform_filter(z,3); mu=ndimage.uniform_filter(res,5); RI=np.sqrt(np.maximum(ndimage.uniform_filter(res**2,5)-mu**2,0))
Ws={'IC_C_sealed1':Wmap(1.0),'IC_C_sealed0.15':Wmap(0.15),'IC_rough':np.clip(1-RI/np.nanmax(RI[land]),0.001,1),'IC_topo':np.ones(land.shape)}
S=np.clip(np.tan(np.radians(f['slope'])),0.005,1); sink=land&((fr[80]+fr[90])>=0.5)
di,dj,dist=downstream(); mask=land&~sink; order,rec=topo_order(di,dj,mask)
n=land.size; Sf=S.ravel(); dl=dist.ravel(); out={}
for name,W in Ws.items():
    Wf=W.ravel(); sW=Wf.copy(); sS=Sf.copy(); cnt=np.ones(n)
    for i in order:
        j=rec[i]
        if j>=0: sW[j]+=sW[i]; sS[j]+=sS[i]; cnt[j]+=cnt[i]
    Dup=(sW/cnt)*(sS/cnt)*np.sqrt(cnt*cell*cell); Ddn=np.zeros(n)
    for i in order[::-1]:
        j=rec[i]; Ddn[i]=dl[i]/(Wf[i]*Sf[i])+(Ddn[j] if j>=0 else 0)
    with np.errstate(divide='ignore',invalid='ignore'): IC=np.log10(Dup/Ddn).reshape(land.shape)
    IC[~mask]=np.nan; out[name]=IC
np.savez(f'{OUT}/ic.npz',**out,sink=sink); print('ok')
