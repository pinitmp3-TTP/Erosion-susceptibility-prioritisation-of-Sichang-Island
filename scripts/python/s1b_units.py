import numpy as np, pandas as pd, geopandas as gpd, sys
from rasterio import features
from scipy import ndimage
from shapely.geometry import shape
from common import *
THR=float(sys.argv[1]) if len(sys.argv)>1 else 5.0
outlet=np.load(f'{OUT}/outlet.npy')
u,cnt=np.unique(outlet[land],return_counts=True); ha=cnt*cell*cell/1e4
big=u[ha>=THR]; unit=np.zeros(land.shape,int)
cen=[(np.mean(np.nonzero(outlet==o)[0]),o) for o in big]
for k,(_,o) in enumerate(sorted(cen),1): unit[outlet==o]=k
face=land&(unit==0); lab,n=ndimage.label(face,structure=np.ones((3,3)))
fs=ndimage.sum(np.ones_like(lab),lab,range(1,n+1))*cell*cell/1e4
cenf=sorted([(ndimage.center_of_mass(lab==i)[0],i) for i in range(1,n+1) if fs[i-1]>=1.0])
for k,(_,i) in enumerate(cenf,1): unit[lab==i]=1000+k
for i in range(1,n+1):
    if fs[i-1]>=1.0: continue
    m=lab==i; ring=ndimage.binary_dilation(m,np.ones((3,3)))&~m&land&(unit>0)
    v,c=np.unique(unit[ring],return_counts=True); unit[m]=v[np.argmax(c)]
assert (unit[land]>0).all()
geoms=[(int(v),shape(g)) for g,v in features.shapes(unit.astype('int32'),mask=land,transform=tr)]
gdf=gpd.GeoDataFrame({'uid':[v for v,_ in geoms]},geometry=[g for _,g in geoms],crs=prof['crs']).dissolve('uid').reset_index()
gdf['kind']=np.where(gdf.uid<1000,'basin','face'); gdf['label']=[f'B{v}' if v<1000 else f'F{v-1000}' for v in gdf.uid]
gdf['area_ha']=[(unit==v).sum()*cell*cell/1e4 for v in gdf.uid]
gdf.to_file(f'{OUT}/units_{THR:g}ha.shp'); np.save(f'{OUT}/unit_{THR:g}ha.npy',unit)
print(f'{len(big)} basins, {len(cenf)} faces')
