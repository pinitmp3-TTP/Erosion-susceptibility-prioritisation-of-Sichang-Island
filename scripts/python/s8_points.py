import numpy as np, pandas as pd, geopandas as gpd
from shapely.geometry import Point
from common import *
rng=np.random.default_rng(20261114)
unit=np.load(f'{OUT}/unit_5ha.npy'); U=gpd.read_file(f'{OUT}/units_5ha.shp').set_index('uid')
T=pd.read_csv(f'{OUT}/units_5ha_mc.csv').set_index('uid')
codes={u:f'U{k:02d}' for k,u in enumerate(rng.permutation(U.index.values),1)}
rows=[]
for uid,g in U.iterrows():
    n=10 if g.area_ha>=3 else 6
    poly=g.geometry.buffer(-15)
    if poly.is_empty or poly.area<1: poly=g.geometry
    minx,miny,maxx,maxy=poly.bounds; pts=[]; tries=0
    while len(pts)<n and tries<20000:
        tries+=1; p=Point(rng.uniform(minx,maxx),rng.uniform(miny,maxy))
        if poly.contains(p) and all(p.distance(q)>=40 for q in pts): pts.append(p)
    for p in pts: rows.append(dict(uid=uid,label=g.label,ucode=codes[uid],x=p.x,y=p.y))
P=gpd.GeoDataFrame(rows,geometry=[Point(r['x'],r['y']) for r in rows],crs=U.crs)
P=P.sample(frac=1,random_state=7).reset_index(drop=True); P['pid']=[f'P{k:03d}' for k in range(1,len(P)+1)]
W=P.to_crs(4326); P['lon']=W.geometry.x.round(6); P['lat']=W.geometry.y.round(6)
key=P[['pid','uid','label','ucode','lon','lat','x','y']].merge(T[['cls','CS','P_top']],left_on='uid',right_index=True)
key.to_csv(f'{OUT}/random_points_key_regenerated.csv',index=False)
P[['pid','lon','lat','geometry']].to_file(f'{OUT}/inv_points.gpkg')
print(len(P))
