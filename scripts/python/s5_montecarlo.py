# Monte Carlo analysis of ranking stability (1000 draws)
import numpy as np, pandas as pd, geopandas as gpd
from common import *
N=1000; rng=np.random.default_rng(20261101)
unit=np.load(f'{OUT}/unit_5ha.npy'); U=gpd.read_file(f'{OUT}/units_5ha.shp').set_index('uid')
f=np.load(f'{OUT}/factors.npz'); cf=np.load(f'{OUT}/cover_frac.npz'); T=pd.read_csv(f'{OUT}/units_5ha_table.csv').set_index('uid')
uids=U.index.values; labels=U.label.values; R=float(f['R'])
erod=np.clip(1-cf['built']-cf['water'],0,1)
CLS=[10,20,30,40,60,95]
CR={10:(0.003,0.043),20:(0.032,0.05),30:(0.015,0.05),40:(0.15,0.25),60:(0.45,0.8),95:(0.0,0.003)}
Ks={0.5:f['K_OM05'],1.0:f['K'],2.0:f['K_OM2']}; LSs={'MB':f['LS_MB'],'DG':f['LS_DG'],'DGcap':f['LS_DGcap']}
# precompute S[ls][om][class] = per-unit sum(R*K*LS*frac_class)/erodible area  (supply is linear in each class C)
E=np.array([erod[(unit==u)&land].sum() for u in uids])
S={ls:{om:np.array([[np.nansum((R*K*LS*cf[f'f{c}'])[(unit==u)&land]) for u in uids] for c in CLS])/E for om,K in Ks.items()} for ls,LS in LSs.items()}
ICV=['IC_C_sealed1','IC_C_sealed0.15','IC_rough','IC_topo']; ICm=T.loc[uids,ICV].values
z=lambda x:(x-x.mean())/x.std(ddof=0)
top=np.zeros(len(uids)); ranks=np.zeros((N,len(uids)))
for i in range(N):
    om=rng.uniform(0.5,2.0); ls=rng.choice(list(LSs)); icv=rng.integers(4); w=rng.uniform(0.3,0.6); lg=rng.random()<0.5
    Cv=np.array([rng.uniform(*CR[c]) if CR[c][0]==0 else np.exp(rng.uniform(*np.log(CR[c]))) for c in CLS])
    # K interpolated linearly in OM between the 0.5/1.0/2.0 % rasters
    lo,hi=(0.5,1.0) if om<=1.0 else (1.0,2.0); a=(om-lo)/(hi-lo)
    sup=Cv@((1-a)*S[ls][lo]+a*S[ls][hi])
    if lg: sup=np.log(sup)
    cs=w*z(sup)+(1-w)*z(ICm[:,icv])
    r=pd.Series(cs).rank(ascending=False,method='first').values; ranks[i]=r
    top+=pd.qcut(pd.Series(cs).rank(method='first'),5,labels=False).values==4
out=pd.DataFrame(dict(uid=uids,label=labels,P_top=top/N,rank_med=np.median(ranks,0),rank_p5=np.percentile(ranks,5,0),rank_p95=np.percentile(ranks,95,0)))
out.to_csv(f'{OUT}/units_5ha_mc_rerun.csv',index=False); print(out.sort_values('P_top',ascending=False).head(12).to_string())
