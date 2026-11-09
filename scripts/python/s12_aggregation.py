# Aggregation test: unit supply as mean (main), 90th percentile, and exceedance area
# (share of erodible area with A above the island-wide 90th percentile of the main A),
# for the main and the disturbance-aware C. IC unchanged. Evaluated as in s11.
import numpy as np, pandas as pd, geopandas as gpd, rasterio, json
from rasterio.warp import reproject, Resampling
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
from common import *
f=np.load(f'{OUT}/factors.npz'); cf=np.load(f'{OUT}/cover_frac.npz'); unit=np.load(f'{OUT}/unit_5ha.npy')
T=pd.read_csv(f'{OUT}/units_5ha_scored.csv').set_index('uid')
s=rasterio.open(BARE_FILE); a=s.read(BARE_BAND).astype('float32'); fb=np.full(land.shape,np.nan,'float32')
reproject(a,fb,src_transform=s.transform,src_crs=s.crs,dst_transform=tr,dst_crs=prof['crs'],resampling=Resampling.average,src_nodata=np.nan,dst_nodata=np.nan)
fb=np.nan_to_num(fb.astype(float))
erod=np.clip(1-cf['built']-cf['water'],0,1); Cveg=np.where(erod>0,f['C']/np.maximum(erod,1e-9),0)
RKLS=f['R']*f['K']*f['LS_DGcap']
Ae={'main':np.where(erod>0,RKLS*Cveg,np.nan),'dist':np.where(erod>0,RKLS*((1-fb)*Cveg+fb*0.80),np.nan)}   # per erodible area
thr=np.nanpercentile(Ae['main'][land&(erod>0)],90)
def unitscore(A,how):
    out={}
    for u in T.index:
        m=(unit==u)&land&(erod>0); v=A[m]; w=erod[m]
        if how=='mean': out[u]=np.nansum(v*w)/w.sum()
        elif how=='p90': out[u]=np.nanpercentile(v,90)
        else: out[u]=w[v>thr].sum()/w.sum()
    return pd.Series(out)
z=lambda v:(v-v.mean())/v.std(ddof=0)
tb=pd.read_csv(f'{OUT}/turbidity_unit_response.csv').set_index('label')['resp']
P=pd.read_csv(f'{OUT}/randompts_AB.csv'); dsh=P.groupby('label').D.mean()
key=pd.read_csv(INTERP+'/random_points_key.csv').set_index('pid'); hasE=set(P.pid[P.E==1]); neg=key.loc[[p for p in P.pid.unique() if p not in hasE]]
def uat(x,y): return unit[((np.asarray(y)-tr.f)/tr.e).astype(int),((np.asarray(x)-tr.c)/tr.a).astype(int)]
un=uat(neg.x,neg.y); FA=gpd.read_file(INTERP+'/erosion_features_interpreter1.gpkg'); FB=gpd.read_file(INTERP+'/erosion_features_interpreter2.gpkg')
ua,ub=uat(FA.geometry.x,FA.geometry.y),uat(FB.geometry.x,FB.geometry.y)
def auc(cs,uf): 
    m=pd.Series(cs); pos=m.reindex(uf).values; ng=m.reindex(un).values
    return round(roc_auc_score(np.r_[np.ones(len(pos)),np.zeros(len(ng))],np.r_[pos,ng]),2)
top0=set(T.index[T.cls==5]); res={}
for c in ['main','dist']:
    for how in ['mean','p90','exceed']:
        sup=unitscore(Ae[c],how); cs=0.5*z(sup)+0.5*z(T.IC_C_sealed1); top=set(cs.sort_values(ascending=False).index[:8])
        lab=T.label; csl=pd.Series(cs.values,index=lab.values)
        res[f'{c}_{how}']=dict(top8_shared=len(top&top0),new=sorted(lab[list(top-top0)]),rho_CS_main=round(spearmanr(cs,T.CS)[0],2),
            rho_sup_slope=round(spearmanr(sup,T.slope)[0],2),turb=round(spearmanr(csl.reindex(tb.index),tb)[0],2),
            disturbed=round(spearmanr(csl.reindex(dsh.index),dsh)[0],2),auc_I1=auc(cs,ua),auc_I2=auc(cs,ub))
R=json.load(open(f'{OUT}/v4_results.json')); R['aggregation_test']=dict(threshold_t_ha_yr=round(float(thr),1),**res); json.dump(R,open(f'{OUT}/v4_results.json','w'),indent=1,ensure_ascii=False)
print('thr',round(thr,1)); print(pd.DataFrame(res).T.to_string())
