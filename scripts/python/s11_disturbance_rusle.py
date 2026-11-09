# Disturbance-aware supply (method fixed before results were examined):
# on the erodible (non-sealed) fraction of each cell, C' = (1 - f) * C_veg + f * 0.80,
# where f is the Sentinel-2 bare-soil frequency (NDVI < 0.20, BSI > 0) and 0.80 the LDD value for bare land.
# IC is kept unchanged. The new ranking is evaluated against nearshore turbidity, interpreted disturbance
# and the two feature inventories, but not against bare-soil frequency (which it now contains).
import numpy as np, pandas as pd, geopandas as gpd, rasterio, json
from rasterio.warp import reproject, Resampling
from scipy.stats import spearmanr, mannwhitneyu
from sklearn.metrics import roc_auc_score
from common import *
f=np.load(f'{OUT}/factors.npz'); cf=np.load(f'{OUT}/cover_frac.npz'); unit=np.load(f'{OUT}/unit_5ha.npy')
T=pd.read_csv(f'{OUT}/units_5ha_scored.csv').set_index('uid')
s=rasterio.open(BARE_FILE); a=s.read(BARE_BAND).astype('float32'); fb=np.full(land.shape,np.nan,'float32')
reproject(a,fb,src_transform=s.transform,src_crs=s.crs,dst_transform=tr,dst_crs=prof['crs'],resampling=Resampling.average,src_nodata=np.nan,dst_nodata=np.nan)
fb=np.nan_to_num(fb.astype(float))
erod=np.clip(1-cf['built']-cf['water'],0,1); Cveg=np.where(erod>0,f['C']/np.maximum(erod,1e-9),0)
Cd=erod*((1-fb)*Cveg+fb*0.80)
A=f['R']*f['K']*f['LS_DGcap']*Cd
sup=pd.Series({u:np.nansum(A[(unit==u)&land])/max(erod[(unit==u)&land].sum(),1e-6) for u in T.index})
z=lambda v:(v-v.mean())/v.std(ddof=0)
CSd=0.5*z(sup)+0.5*z(T.IC_C_sealed1); rk=CSd.rank(ascending=False)
topD=set(T.label[rk<=8]); top0=set(T.label[T.cls==5]); T['CSd']=CSd; T['supd']=sup
X=dict(C_mean_island=round(float(np.nanmean(f['C'][land])),4),Cd_mean_island=round(float(np.nanmean(Cd[land])),4),
       rho_supply_with_main=round(spearmanr(sup,T.Ae_DGcap)[0],2),rho_CS_with_main=round(spearmanr(CSd,T.CS)[0],2),
       top8_shared=len(topD&top0),top8=sorted(topD),new_in_top=sorted(topD-top0),
       rho_supply_slope=round(spearmanr(sup,T.slope)[0],2))
# evaluation 1: turbidity response
tb=pd.read_csv(f'{OUT}/turbidity_unit_response.csv').set_index('label'); TT=T.set_index('label').join(tb[['resp']])
X['turbidity']=dict(rho_main=round(spearmanr(TT.resp,TT.CS)[0],2),rho_dist=round(spearmanr(TT.resp,TT.CSd)[0],2),p_dist=round(spearmanr(TT.resp,TT.CSd)[1],3),
    top_med=round(TT.resp[TT.index.isin(topD)].median(),2),
    rest_med=round(TT.resp[~TT.index.isin(topD)].median(),2))
# evaluation 2: interpreted disturbance (unit share of disturbed points)
P=pd.read_csv(f'{OUT}/randompts_AB.csv'); g=P.groupby('label').D.mean().to_frame().join(TT[['CS','CSd']])
X['disturbed_share']=dict(rho_main=round(spearmanr(g.D,g.CS)[0],2),rho_dist=round(spearmanr(g.D,g.CSd)[0],2),p_dist=round(spearmanr(g.D,g.CSd)[1],3))
# evaluation 3: feature inventories (unit compound score at features vs at 189 negatives)
key=pd.read_csv(INTERP+'/random_points_key.csv').set_index('pid')
hasE=set(P.pid[P.E==1]); neg=key.loc[[p for p in P.pid.unique() if p not in hasE]]
def unit_at(x,y):
    rr=((np.asarray(y)-tr.f)/tr.e).astype(int); cc=((np.asarray(x)-tr.c)/tr.a).astype(int); return unit[rr,cc]
mapU=lambda u,col:pd.Series(u).map(T[col]).values
rng=np.random.default_rng(3)
def auc(pos,negv):
    y=np.r_[np.ones(len(pos)),np.zeros(len(negv))]; sc=np.r_[pos,negv]; ok=np.isfinite(sc); a=roc_auc_score(y[ok],sc[ok]); bs=[]
    for _ in range(2000):
        i=rng.integers(0,len(pos),len(pos)); j=rng.integers(0,len(negv),len(negv)); yy=np.r_[np.ones(len(i)),np.zeros(len(j))]; ss=np.r_[pos[i],negv[j]]; k=np.isfinite(ss)
        if len(set(yy[k]))==2: bs.append(roc_auc_score(yy[k],ss[k]))
    return [round(a,2),round(np.percentile(bs,2.5),2),round(np.percentile(bs,97.5),2)]
un=unit_at(neg.x,neg.y)
X['inventory_auc']={}
for nm,fn in [('I1','erosion_features_interpreter1.gpkg'),('I2','erosion_features_interpreter2.gpkg')]:
    F=gpd.read_file(INTERP+'/'+fn); uf=unit_at(F.geometry.x,F.geometry.y)
    X['inventory_auc'][nm]=dict(main=auc(mapU(uf,'CS'),mapU(un,'CS')),disturbance_aware=auc(mapU(uf,'CSd'),mapU(un,'CSd')))
R=json.load(open(f'{OUT}/v4_results.json')); R['disturbance_aware']=X; json.dump(R,open(f'{OUT}/v4_results.json','w'),indent=1,ensure_ascii=False)
print(json.dumps(X,indent=1))
