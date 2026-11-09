# Export the main rasters (GeoTIFF, analysis grid, UTM 47N) and unit tables for users of the data set
import numpy as np, pandas as pd, geopandas as gpd, rasterio
from common import *
f=np.load(f'{OUT}/factors.npz'); ic=np.load(f'{OUT}/ic.npz'); cf=np.load(f'{OUT}/cover_frac.npz')
R=float(f['R']); os.makedirs(f'{OUT}/rasters',exist_ok=True)
p=prof.copy(); p.update(count=1,dtype='float32',nodata=-9999,compress='deflate')
lay={'K_OM1':f['K'],'C_LDD':f['C'],'LS_DGcap':f['LS_DGcap'],'LS_DG':f['LS_DG'],'LS_MB':f['LS_MB'],
     'A_DGcap':R*f['K']*f['LS_DGcap']*f['C'],'A_DG':R*f['K']*f['LS_DG']*f['C'],'A_MB':R*f['K']*f['LS_MB']*f['C'],
     'built_fraction':cf['built'],'water_fraction':cf['water']}
lay.update({k:ic[k] for k in ['IC_C_sealed1','IC_C_sealed0.15','IC_rough','IC_topo']})
for k,a in lay.items():
    a=np.where(land&np.isfinite(a),a,-9999).astype('float32')
    with rasterio.open(f'{OUT}/rasters/{k}.tif','w',**p) as o: o.write(a,1)
T=pd.read_csv(f'{OUT}/units_5ha_scored.csv'); key=pd.read_csv(INTERP+'/random_points_key.csv')[['uid','ucode']].drop_duplicates()
U=gpd.read_file(f'{OUT}/units_5ha.shp')[['uid','geometry']].merge(T,on='uid').merge(key,on='uid',how='left')
U.to_file(f'{OUT}/units_5ha_scored.gpkg')
for t in ['2','10']: gpd.read_file(f'{OUT}/units_{t}ha.shp').to_file(f'{OUT}/units_{t}ha.gpkg')
print('exported',len(lay),'rasters;',len(U),'units')
