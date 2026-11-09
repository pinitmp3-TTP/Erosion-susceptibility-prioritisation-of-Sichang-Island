import numpy as np, pandas as pd, geopandas as gpd, rasterio, sys
from rasterio.warp import reproject, Resampling
from shapely.ops import unary_union
from shapely.geometry import shape
from rasterio import features
from common import *
THR=sys.argv[1] if len(sys.argv)>1 else '5'
Uu=D+'/'
unit=np.load(f'{OUT}/unit_{THR}ha.npy'); U=gpd.read_file(f'{OUT}/units_{THR}ha.shp').set_index('uid')
f=np.load(f'{OUT}/factors.npz'); ic=np.load(f'{OUT}/ic.npz'); cf=np.load(f'{OUT}/cover_frac.npz')
dem,_,_=rd('dem'); fa=f['fa']; sl=f['slope']
def regrid(path,res=Resampling.average,band=1):
    with rasterio.open(path) as r: a=r.read(band).astype('float32'); t=r.transform; c=r.crs; nd=r.nodata
    if nd is not None: a[a==nd]=np.nan
    o=np.full(land.shape,np.nan,'float32'); reproject(a,o,src_transform=t,src_crs=c,dst_transform=tr,dst_crs=prof['crs'],resampling=res,src_nodata=np.nan,dst_nodata=np.nan); return o.astype(float)
Cn=regrid(Uu+'Sichang_C_NDVI_2021_30m.tif',Resampling.nearest)
Cn=np.where(np.isfinite(Cn),Cn*np.clip(1-cf['built']-cf['water'],0,1),np.nan)   # sealed surfaces C=0
bare=regrid(BARE_FILE,band=BARE_BAND); R,K,C=f['R'],f['K'],f['C']
erod=np.clip(1-cf['built']-cf['water'],0,1)
LSs={'MB':f['LS_MB'],'DG':f['LS_DG'],'DGcap':f['LS_DGcap']}
diag=np.isin(fdir,[2,8,32,128]); step=np.where(diag,cell*np.sqrt(2),cell)
island=unary_union([shape(g) for g,v in features.shapes(land.astype('uint8'),mask=land,transform=tr) if v==1])
rows=[]
for uid,u in U.iterrows():
    m=(unit==uid)&land; r=dict(uid=uid,label=u.label,kind=u.kind,area_ha=m.sum()*cell*cell/1e4); e=erod[m]; r['sealed_frac']=1-e.mean()
    for k,LS in LSs.items():
        Ax=R*K*LS*C; r[f'A_{k}']=np.nanmean(Ax[m]); r[f'Ae_{k}']=np.nansum(Ax[m])/max(e.sum(),1e-6)
    r['Ae_DGcap_K05']=np.nansum((R*f['K_OM05']*LSs['DGcap']*C)[m])/e.sum(); r['Ae_DGcap_K2']=np.nansum((R*f['K_OM2']*LSs['DGcap']*C)[m])/e.sum()
    r['Ae_DGcap_ndvi']=np.nansum((R*K*LSs['DGcap']*Cn)[m])/e.sum()
    for k in ['IC_C_sealed1','IC_C_sealed0.15','IC_rough','IC_topo']: r[k]=np.nanmean(ic[k][m])
    bm=m&(cf['built']<0.1); r['bare_nobuilt']=np.nanmean(bare[bm]) if bm.sum() else np.nan; r['bare_all']=np.nanmean(bare[m]); r['slope']=np.nanmean(sl[m])
    if u.kind=='basin':
        poly=u.geometry.intersection(island); A_km2=m.sum()*cell*cell/1e6; P=poly.length/1e3
        hull=np.array(poly.convex_hull.exterior.coords); Lb=max(np.hypot(*(hull[j]-hull[k])) for j in range(len(hull)) for k in range(j+1,len(hull)))/1e3
        z=dem[m]; H=np.nanmax(z)-np.nanmin(z); L=step[m&(fa>20)].sum()/1e3
        r.update(Dd=L/A_km2,H=H,Rh=H/(Lb*1000),HI=(np.nanmean(z)-np.nanmin(z))/H,Rn=H/1000*L/A_km2,Ff=A_km2/Lb**2,Rc=4*np.pi*A_km2/P**2)
    rows.append(r)
pd.DataFrame(rows).set_index('uid').to_csv(f'{OUT}/units_{THR}ha_table.csv'); print('ok',THR)
