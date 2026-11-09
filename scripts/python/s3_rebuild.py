import numpy as np, pandas as pd, rasterio
from rasterio.warp import reproject, Resampling
from scipy.spatial import cKDTree
from scipy import ndimage
from pyproj import Transformer
from common import *
Uu=D+'/'
fa,_,_=rd('flowacc'); sl,_,_=rd('slope'); dem,_,_=rd('dem')
# ---- WorldCover class fractions on analysis grid
with rasterio.open(Uu+'Sichang_C_WorldCover_30m.tif') as r: t30=r.transform; crs=r.crs; sh=r.shape
with rasterio.open(Uu+'Sichang_WorldCover2021_10m.tif') as r: cl=r.read(1); t10=r.transform
def to_grid(a10):
    g=np.full(sh,np.nan,'float32'); reproject(a10.astype('float32'),g,src_transform=t10,src_crs=crs,dst_transform=t30,dst_crs=crs,resampling=Resampling.average,src_nodata=np.nan,dst_nodata=np.nan)
    o=np.full(land.shape,np.nan,'float32'); reproject(g,o,src_transform=t30,src_crs=crs,dst_transform=tr,dst_crs=crs,resampling=Resampling.nearest,src_nodata=np.nan,dst_nodata=np.nan); return np.nan_to_num(o.astype(float))
fr={k:to_grid(cl==k) for k in [10,20,30,40,50,60,80,90,95]}
LDD={10:.020,20:.032,30:.015,40:.150,60:.800,95:0.0}
C=sum(fr[k]*v for k,v in LDD.items())
np.savez(f'{OUT}/cover_frac.npz',built=fr[50],water=fr[80]+fr[90],tree=fr[10],**{f'f{k}':fr[k] for k in fr})
# ---- K nearest soil point
S=pd.read_csv(Uu+'soil_samples.csv')
for c in ['K_si_OM0.5','K_si_OM1.0','K_si_OM2.0']: S.loc[S[c]<=0,c]=0.0057
tf=Transformer.from_crs(4326,32647,always_xy=True); ex,ny=tf.transform(S.E.values,S.N.values)
rr,cc=np.indices(A.shape); X=tr.c+(cc+.5)*tr.a; Y=tr.f+(rr+.5)*tr.e
_,idx=cKDTree(np.c_[ex,ny]).query(np.c_[X.ravel(),Y.ravel()]); idx=idx.reshape(A.shape)
K=S['K_si_OM1.0'].values[idx]; K05=S['K_si_OM0.5'].values[idx]; K2=S['K_si_OM2.0'].values[idx]
R=6240.0; th=np.radians(sl)
LS_MB=((fa+1)*cell/22.13)**0.4*(np.sin(th)/0.0896)**1.3
z=np.where(np.isfinite(dem),dem,np.nanmin(dem)); p=np.pad(z,1,mode='edge')
dzdx=((p[:-2,2:]+2*p[1:-1,2:]+p[2:,2:])-(p[:-2,:-2]+2*p[1:-1,:-2]+p[2:,:-2]))/(8*cell)
dzdy=((p[2:,:-2]+2*p[2:,1:-1]+p[2:,2:])-(p[:-2,:-2]+2*p[:-2,1:-1]+p[:-2,2:]))/(8*cell)
asp=np.arctan2(dzdy,-dzdx); x=np.abs(np.sin(asp))+np.abs(np.cos(asp))
def LS_DG(cap=300.0):
    Ain=fa*cell*cell
    if cap is not None: Ain=np.minimum(Ain,cap*cell*x)
    s=np.sin(th); beta=(s/0.0896)/(3*np.power(s,0.8)+0.56); m=beta/(1+beta)
    L=((Ain+cell**2)**(m+1)-Ain**(m+1))/(cell**(m+2)*x**m*22.13**m)
    S_=np.where(np.tan(th)<0.09,10.8*s+0.03,16.8*s-0.5); return L*S_
np.savez(f'{OUT}/factors.npz',K=K,K_OM05=K05,K_OM2=K2,C=C,R=R,LS_MB=LS_MB,LS_DG=LS_DG(None),LS_DGcap=LS_DG(300.0),slope=sl,fa=fa)
# ---- IC (W = C, sealed = 1)
vegfrac=sum(fr[k] for k in LDD)+fr[50]
W=np.clip((C+fr[50]*1.0)/np.maximum(vegfrac,1e-6),0.001,1)
Sg=np.clip(np.tan(th),0.005,1); water=fr[80]+fr[90]; sink=land&(water>=0.5)
di,dj,dist=downstream(); mask=land&~sink; order,rec=topo_order(di,dj,mask)
n=land.size; Wf=W.ravel(); Sf=Sg.ravel(); dl=dist.ravel(); sW=Wf.copy(); sS=Sf.copy(); cnt=np.ones(n)
for i in order:
    j=rec[i]
    if j>=0: sW[j]+=sW[i]; sS[j]+=sS[i]; cnt[j]+=cnt[i]
Dup=(sW/cnt)*(sS/cnt)*np.sqrt(cnt*cell*cell); Ddn=np.zeros(n)
for i in order[::-1]:
    j=rec[i]; Ddn[i]=dl[i]/(Wf[i]*Sf[i])+(Ddn[j] if j>=0 else 0)
with np.errstate(divide='ignore',invalid='ignore'): IC=np.log10(Dup/Ddn).reshape(land.shape)
IC[~mask]=np.nan; np.savez(f'{OUT}/ic.npz',IC_C_sealed1=IC,sink=sink)
# ---- verify per-unit vs saved table
unit=np.load(f'{OUT}/unit_5ha.npy'); T=pd.read_csv(f'{OUT}/units_5ha_mc.csv').set_index('uid')
erod=np.clip(1-fr[50]-water,0,1); Ax=R*K*LS_DG(300.0)*C
rows=[]
for u in T.index:
    m=(unit==u)&land; rows.append((u,np.nansum(Ax[m])/erod[m].sum(),np.nanmean(IC[m])))
V=pd.DataFrame(rows,columns=['uid','Ae','IC']).set_index('uid').join(T[['Ae_DGcap','IC_C_sealed1']])
from scipy.stats import spearmanr
print('supply rel diff median %.3f max %.3f; rho %.4f'%((V.Ae/V.Ae_DGcap-1).abs().median(),(V.Ae/V.Ae_DGcap-1).abs().max(),spearmanr(V.Ae,V.Ae_DGcap)[0]))
print('IC abs diff max %.3f; rho %.4f'%((V.IC-V.IC_C_sealed1).abs().max(),spearmanr(V.IC,V.IC_C_sealed1)[0]))
