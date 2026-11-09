import numpy as np, pandas as pd, geopandas as gpd, rasterio, json
import sys, os
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.join(HERE,'..','..','python'))
from common import *
A0=A
F=np.load(f'{OUT}/factors.npz'); I=np.load(f'{OUT}/ic.npz'); cf=np.load(f'{OUT}/cover_frac.npz')
A=F['R']*F['K']*F['LS_DGcap']*F['C']; A[~land]=np.nan
IC=I['IC_C_sealed1'].copy(); IC[~land]=np.nan
from rasterio.warp import reproject, Resampling
def regrid(path,res=Resampling.average,band=1):
    with rasterio.open(path) as r: a=r.read(band).astype('float32'); t=r.transform; c=r.crs; nd=r.nodata
    if nd is not None: a[a==nd]=np.nan
    o=np.full(land.shape,np.nan,'float32'); reproject(a,o,src_transform=t,src_crs=c,dst_transform=tr,dst_crs=prof['crs'],resampling=res,src_nodata=np.nan,dst_nodata=np.nan); return o.astype(float)
bare=regrid(BARE_FILE,band=BARE_BAND).astype(float); bare[~land]=np.nan; bare[cf['built']>=0.10]=np.nan
p=prof.copy(); p.update(count=3,dtype='float32',nodata=-9999,compress='deflate')
with rasterio.open(HERE+'/sichang_rasters_v4.tif','w',**p) as o:
    for i,a in enumerate([A,IC,bare],1):
        o.write(np.where(np.isfinite(a),a,-9999).astype('float32'),i)
print('A mean',np.nanmean(A),'IC',np.nanmin(IC),np.nanmax(IC),'bare',np.nanmax(bare))
U=gpd.read_file(f'{OUT}/units_5ha.shp'); S=pd.read_csv(f'{OUT}/units_5ha_scored.csv')
tb=pd.read_csv(f'{OUT}/turbidity_unit_response.csv')[['label','resp']]
S=S.merge(tb,on='label',how='left'); S['rank']=S.CS.rank(ascending=False).astype(int)
names={1:'Very low',2:'Low',3:'Moderate',4:'High',5:'Very high'}
S['prio_lab']=S.cls.map(names); S['driver']=S.driver.map({'supply+transfer':'Supply + transfer','supply':'Supply','transfer':'Transfer'})
print(U.columns.tolist())
key='label' if 'label' in U.columns else 'uid'
keep=dict(label='label',kind='kind',area_ha='area_ha',slope='slope_deg',Ae_DGcap='supply',IC_C_sealed1='IC',zS='z_supply',zT='z_trans',CS='CS',cls='priority',prio_lab='prio_lab',rank='rank',P_top='P_top',driver='driver',bare_nobuilt='bare_freq',resp='turb_resp',Dd='Dd',H='H',Rh='Rh',HI='HI',Ff='Ff',Rc='Rc')
X=S[list(keep)].rename(columns=keep)
G=U[[key,'geometry']].merge(X,left_on=key,right_on='label' if key=='label' else 'label',how='inner') if key=='label' else U[['uid','geometry']].merge(S[['uid']+list(keep)].rename(columns=keep),on='uid')
G=G.to_crs(4326); print(len(G), G.columns.tolist())
G.to_file(HERE+'/sichang_units_v4.shp')
