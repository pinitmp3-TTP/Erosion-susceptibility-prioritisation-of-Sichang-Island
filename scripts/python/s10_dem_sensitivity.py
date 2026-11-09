# DEM sensitivity: repeat flow routing, slope, LS, IC and the ranking with FABDEM (forests and buildings removed),
# keeping the 41 units fixed. A control run with GLO-30 from the same export checks that the routing itself
# reproduces the main ranking. Uses scripts s3_rebuild.py, s2_icvariants.py and s4_table.py unchanged.
import numpy as np, pandas as pd, rasterio, os, shutil, subprocess, sys, heapq, json, glob
from scipy.stats import spearmanr
HERE=os.path.dirname(os.path.abspath(__file__)); BASE=os.path.dirname(os.path.dirname(HERE))
IN=os.path.join(BASE,'inputs'); OUT=os.path.join(BASE,'outputs')
ref=rasterio.open(os.path.join(IN,'dem.tif')); prof=ref.profile; nr,nc=ref.shape; cell=ref.transform.a
land=rasterio.open(os.path.join(IN,'grid_land_mask.tif')).read(1)==1
src=rasterio.open(os.path.join(IN,'Sichang_FABDEM_grid.tif'))
D8=[(1,0,1),(2,1,1),(4,1,0),(8,1,-1),(16,0,-1),(32,-1,-1),(64,-1,0),(128,-1,1)]   # ESRI code, drow, dcol
def hydro(z):
    z=np.where(np.isfinite(z),z,0.0).astype(float)
    # priority-flood fill of land cells, draining to the sea (non-land cells)
    f=np.full(z.shape,np.inf); done=~land.copy(); pq=[]
    for r,c in zip(*np.nonzero(land)):
        if r in (0,nr-1) or c in (0,nc-1) or any(not land[min(max(r+a,0),nr-1),min(max(c+b,0),nc-1)] for _,a,b in D8):
            f[r,c]=z[r,c]; done[r,c]=True; heapq.heappush(pq,(z[r,c],r,c))
    while pq:
        h,r,c=heapq.heappop(pq)
        for _,a,b in D8:
            rr,cc=r+a,c+b
            if 0<=rr<nr and 0<=cc<nc and not done[rr,cc]:
                f[rr,cc]=max(z[rr,cc],h+1e-3); done[rr,cc]=True; heapq.heappush(pq,(f[rr,cc],rr,cc))
    zf=np.where(land,f,np.nanmin(z[land])-1.0)
    fd=np.zeros(z.shape,'int32')
    for r,c in zip(*np.nonzero(land)):
        best=0; code=0
        for k,a,b in D8:
            rr,cc=r+a,c+b
            if 0<=rr<nr and 0<=cc<nc:
                drop=(zf[r,c]-zf[rr,cc])/(np.hypot(a,b))
                if drop>best: best,code=drop,k
        fd[r,c]=code
    # flow accumulation (upstream land cells, excluding self)
    fa=np.zeros(z.shape); dn={}
    for k,a,b in D8:
        for r,c in zip(*np.nonzero(land&(fd==k))):
            rr,cc=r+a,c+b
            if 0<=rr<nr and 0<=cc<nc and land[rr,cc]: dn[(r,c)]=(rr,cc)
    indeg={}
    for v in dn.values(): indeg[v]=indeg.get(v,0)+1
    stack=[(r,c) for r,c in zip(*np.nonzero(land)) if indeg.get((r,c),0)==0]; acc=np.ones(z.shape)
    while stack:
        u=stack.pop(); v=dn.get(u)
        if v is None: continue
        acc[v]+=acc[u]; indeg[v]-=1
        if indeg[v]==0: stack.append(v)
    fa=np.where(land,acc-1,0)
    # Horn slope (degrees) on the unfilled DEM
    p=np.pad(z,1,mode='edge')
    dzdx=((p[:-2,2:]+2*p[1:-1,2:]+p[2:,2:])-(p[:-2,:-2]+2*p[1:-1,:-2]+p[2:,:-2]))/(8*cell)
    dzdy=((p[2:,:-2]+2*p[2:,1:-1]+p[2:,2:])-(p[:-2,:-2]+2*p[:-2,1:-1]+p[:-2,2:]))/(8*cell)
    sl=np.degrees(np.arctan(np.hypot(dzdx,dzdy)))
    return z,fd,fa,sl
def run(name,band):
    z=src.read(band)[:nr,:nc].astype(float); z[z<-1000]=np.nan
    z,fd,fa,sl=hydro(z)
    di=os.path.join(OUT,'dem_sens',name,'inputs'); do=os.path.join(OUT,'dem_sens',name,'outputs'); os.makedirs(di,exist_ok=True); os.makedirs(do,exist_ok=True)
    for f in glob.glob(os.path.join(IN,'*')): 
        if not os.path.basename(f).startswith(('dem.','slope.','flowdir.','flowacc.')): shutil.copy(f,di)
    p=prof.copy(); p.update(dtype='float32',nodata=-9999,count=1)
    for nm,a in [('dem',z),('slope',np.where(land,sl,-9999)),('flowacc',fa)]:
        with rasterio.open(os.path.join(di,nm+'.tif'),'w',**p) as o: o.write(a.astype('float32'),1)
    p2=prof.copy(); p2.update(dtype='int32',nodata=0,count=1)
    with rasterio.open(os.path.join(di,'flowdir.tif'),'w',**p2) as o: o.write(fd,1)
    for f in glob.glob(os.path.join(OUT,'unit_5ha.npy'))+glob.glob(os.path.join(OUT,'units_5ha.*'))+[os.path.join(OUT,'units_5ha_mc.csv')]: shutil.copy(f,do)
    env=dict(os.environ,SICHANG_INPUTS=di,SICHANG_OUT=do)
    for s in [['s3_rebuild.py'],['s2_icvariants.py'],['s4_table.py','5']]:
        subprocess.run([sys.executable]+s,cwd=HERE,env=env,check=True,capture_output=True)
    return pd.read_csv(os.path.join(do,'units_5ha_table.csv')).set_index('uid')
def score(t):
    z=lambda v:(v-v.mean())/v.std(ddof=0); cs=0.5*z(t.Ae_DGcap)+0.5*z(t.IC_C_sealed1); return cs
M=pd.read_csv(os.path.join(OUT,'units_5ha_scored.csv')).set_index('uid'); top0=set(M.label[M.cls==5])
res={}; diff=src.read(1)[:nr,:nc]-src.read(2)[:nr,:nc]
res['fab_minus_glo_land']=dict(mean=round(float(diff[land].mean()),2),p5=round(float(np.percentile(diff[land],5)),2),min=round(float(diff[land].min()),2))
for name,band in [('glo30_control',2),('fabdem',1)]:
    t=run(name,band); cs=score(t); top=set(M.label[cs.rank(ascending=False)<=8])
    res[name]=dict(rho_CS=round(spearmanr(cs,M.CS)[0],3),rho_supply=round(spearmanr(t.Ae_DGcap,M.Ae_DGcap)[0],3),rho_IC=round(spearmanr(t.IC_C_sealed1,M.IC_C_sealed1)[0],3),
                   rho_slope=round(spearmanr(t.slope,M.slope)[0],3),top8_shared=len(top&top0),top8=sorted(top),
                   mean_slope_top0=round(t.slope[M.cls==5].mean(),1),mean_slope_rest=round(t.slope[M.cls<5].mean(),1),
                   rho_supply_IC=round(spearmanr(t.Ae_DGcap,t.IC_C_sealed1)[0],2),island_mean_slope=round(float(t.slope.mean()),1))
R=json.load(open(os.path.join(OUT,'v4_results.json'))); R['dem_sensitivity']=res; json.dump(R,open(os.path.join(OUT,'v4_results.json'),'w'),indent=1,ensure_ascii=False)
print(json.dumps(res,indent=1))
