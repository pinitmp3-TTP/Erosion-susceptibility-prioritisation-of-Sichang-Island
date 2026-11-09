# Agreement between the two blind erosion-feature inventories (distance, unit and tile scale)
import numpy as np, pandas as pd, geopandas as gpd, json, re
from shapely.geometry import Polygon
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
from sklearn.metrics import cohen_kappa_score as kap
from common import *
FA=gpd.read_file(INTERP+'/erosion_features_interpreter1.gpkg'); FB=gpd.read_file(INTERP+'/erosion_features_interpreter2.gpkg')
xa=np.c_[FA.geometry.x,FA.geometry.y]; xb=np.c_[FB.geometry.x,FB.geometry.y]
dA,_=cKDTree(xb).query(xa); dB,_=cKDTree(xa).query(xb)          # nearest feature of the other interpreter
NEW='ภาพเก่ายังไม่มี'                                              # "not present on the older image"
inv=dict(A_n=len(FA),B_n=len(FB),A_classes=FA['class'].value_counts().to_dict(),B_classes=FB['class'].value_counts().to_dict(),
         A_new=int((FA.seen==NEW).sum()),B_new=int((FB.seen==NEW).sum()),B_conf_ge2=int(FB.conf.astype(str).str[0].astype(int).ge(2).sum()),
         match30=int((dA<=30).sum()),match50=int((dA<=50).sum()),match100_A=int((dA<=100).sum()),match100_B=int((dB<=100).sum()))
# unit scale (41 units)
unit=np.load(f'{OUT}/unit_5ha.npy'); U=gpd.read_file(f'{OUT}/units_5ha.shp').set_index('uid')
def cnt(F):
    c=((F.geometry.x-tr.c)/tr.a).astype(int).values; r=((F.geometry.y-tr.f)/tr.e).astype(int).values
    u=pd.Series(unit[r,c]); return u[u>=0].value_counts().reindex(U.index,fill_value=0)
ca,cb=cnt(FA),cnt(FB); rho,p=spearmanr(ca,cb)
inv.update(unit_kappa=round(kap(ca>0,cb>0),2),unit_rho=round(rho,2),unit_p=round(p,2),A_units=int((ca>0).sum()),B_units=int((cb>0).sum()),both_units=int(((ca>0)&(cb>0)).sum()))
# tile scale (108 sweep tiles)
k=open(INTERP+'/sweep_tiles_300m.kml',encoding='utf-8').read(); tiles=[]
for name,coords in re.findall(r'<name>(T\d+)</name>.*?<LinearRing><coordinates>(.*?)</coordinates>',k,flags=re.S):
    tiles.append(dict(tile=name,geometry=Polygon([tuple(map(float,c.split(',')[:2])) for c in coords.split()])))
Tl=gpd.GeoDataFrame(tiles,crs=4326).to_crs(FA.crs)
ta=Tl.geometry.apply(lambda g:FA.within(g).any()); tb=Tl.geometry.apply(lambda g:FB.within(g).any())
inv.update(tile_kappa=round(kap(ta,tb),2),tile_both=int((ta&tb).sum()),tiles=len(Tl),A_tiles=int(ta.sum()),B_tiles=int(tb.sum()))
roc=pd.read_csv(f'{OUT}/roc_features_results.csv'); inv['roc']=[{k:(round(v,2) if isinstance(v,float) else v) for k,v in r.items()} for r in roc.to_dict('records')]
R=json.load(open(f'{OUT}/v4_results.json')); R['inventory']=inv; json.dump(R,open(f'{OUT}/v4_results.json','w'),indent=1,ensure_ascii=False)
print(json.dumps({k:v for k,v in inv.items() if k!='roc'},ensure_ascii=False))
