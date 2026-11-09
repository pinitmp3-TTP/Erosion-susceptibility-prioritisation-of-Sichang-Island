import numpy as np, rasterio, os
# Package layout: <BASE>/inputs, <BASE>/outputs, <BASE>/interpretation, <BASE>/scripts/python
BASE=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D=os.path.join(BASE,'inputs'); OUT=os.path.join(BASE,'outputs'); INTERP=os.path.join(BASE,'interpretation')
os.makedirs(OUT,exist_ok=True)
BARE_FILE=os.path.join(D,'S2_bare_frequency_2019_2024_30m.tif'); BARE_BAND=2   # band 2 = freq20 (NDVI < 0.20 and BSI > 0)
def rd(n, d=D):
    with rasterio.open(f'{d}/{n}.tif') as r: a=r.read(1).astype(float); t=r.transform; p=r.profile
    a[(a<-1e30)|(a==-9999)]=np.nan; return a,t,p
A,tr,prof=rd('grid_land_mask'); land=np.isfinite(A)   # analysis grid (30.66 m, UTM 47N) and land mask
fdir=rasterio.open(f'{D}/flowdir.tif').read(1).astype(int)
cell=tr.a
DIRS={1:(0,1),2:(1,1),4:(1,0),8:(1,-1),16:(0,-1),32:(-1,-1),64:(-1,0),128:(-1,1)}
def downstream():
    nr,nc=land.shape; di=np.full(land.shape,-1); dj=np.full(land.shape,-1); dist=np.zeros(land.shape)
    for code,(a,b) in DIRS.items():
        m=fdir==code; r,c=np.nonzero(m); rr,cc=r+a,c+b; ok=(rr>=0)&(rr<nr)&(cc>=0)&(cc<nc)
        di[r[ok],c[ok]]=rr[ok]; dj[r[ok],c[ok]]=cc[ok]; dist[r,c]=cell*(np.sqrt(2) if a and b else 1)
    return di,dj,dist
def topo_order(di,dj,mask):
    nr,nc=mask.shape; idx=np.arange(nr*nc).reshape(mask.shape)
    rec=np.full(nr*nc,-1); r,c=np.nonzero(mask)
    ok=di[r,c]>=0; dn=np.where(ok, idx[np.clip(di[r,c],0,None),np.clip(dj[r,c],0,None)], -1)
    dn[~mask.ravel()[np.where(dn>=0,dn,0)]]=-1
    rec[idx[r,c]]=dn
    indeg=np.zeros(nr*nc,int); np.add.at(indeg, rec[rec>=0],1)
    cur=[i for i in idx[mask] if indeg[i]==0]; order=[]
    while cur:
        order+=cur; nxt=[]
        for i in cur:
            j=rec[i]
            if j>=0:
                indeg[j]-=1
                if indeg[j]==0: nxt.append(j)
        cur=nxt
    return np.array(order), rec
