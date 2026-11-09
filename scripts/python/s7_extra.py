# Additional diagnostics: unit-level tests, IC vs geometry, partial correlation, potential supply, first-round kappa
import numpy as np, pandas as pd, json
from scipy import ndimage
from scipy.stats import spearmanr, mannwhitneyu, t as tdist
from sklearn.metrics import cohen_kappa_score as kap
from common import *
R=json.load(open(f'{OUT}/v4_results.json')); X={}
T=pd.read_csv(f'{OUT}/units_5ha_scored.csv').set_index('uid'); unit=np.load(f'{OUT}/unit_5ha.npy')
f=np.load(f'{OUT}/factors.npz'); cf=np.load(f'{OUT}/cover_frac.npz')
rs=lambda a,b:(round(spearmanr(a,b)[0],2),round(spearmanr(a,b)[1],3))
# 1. unit-level tests of the random-point indicators (8 top-class units vs 33 others)
P=pd.read_csv(f'{OUT}/randompts_AB.csv'); g=P.groupby('label')[['closed','D','bare10','sealed']].mean()
g=g.join(T.set_index('label')[['cls']]); top=g.cls==5
X['unit_level']={c:dict(top_med=round(g[c][top].median(),2),rest_med=round(g[c][~top].median(),2),
    p=round(mannwhitneyu(g[c][top],g[c][~top],alternative='two-sided').pvalue,4)) for c in ['closed','D','bare10','sealed']}
# 2. IC against unit area and mean distance to the coast
dist=ndimage.distance_transform_edt(land)*cell
dc=pd.Series({u:dist[(unit==u)&land].mean() for u in T.index})
X['IC_vs_area']=rs(T.IC_C_sealed1,T.area_ha); X['IC_vs_distcoast']=rs(T.IC_C_sealed1,dc)
X['supply_vs_distcoast']=rs(T.Ae_DGcap,dc)
# 3. partial Spearman correlation of supply and IC controlling for mean slope
def pcor(x,y,z):
    rk=lambda v:pd.Series(v).rank().values; x,y,z=rk(x),rk(y),rk(z)
    rx=x-np.polyval(np.polyfit(z,x,1),z); ry=y-np.polyval(np.polyfit(z,y,1),z)
    r=np.corrcoef(rx,ry)[0,1]; n=len(x); tt=r*np.sqrt((n-3)/(1-r**2)); return round(r,2),round(2*tdist.sf(abs(tt),n-3),3)
X['partial_supply_IC_given_slope']={k:pcor(T.Ae_DGcap,T[k],T.slope) for k in ['IC_C_sealed1','IC_C_sealed0.15','IC_rough','IC_topo']}
# 4. potential supply (C = 1 on erodible land) and the ranking it gives
erod=np.clip(1-cf['built']-cf['water'],0,1); Apot=f['R']*f['K']*f['LS_DGcap']*erod
pot=pd.Series({u:np.nansum(Apot[(unit==u)&land])/erod[(unit==u)&land].sum() for u in T.index})
z=lambda v:(v-v.mean())/v.std(ddof=0)
CSp=0.5*z(pot)+0.5*z(T.IC_C_sealed1); topP=set(T.label[CSp.rank(ascending=False)<=8]); topA=set(T.label[T.cls==5])
X['potential_supply']=dict(rho_actual=rs(pot,T.Ae_DGcap),rho_CS=rs(CSp,T.CS),top8_shared=len(topP&topA),top8_potential=sorted(topP),
    forest_share_C=round(float(np.nansum((cf['f10']*0.020)[land])/np.nansum(f['C'][land])),2))
# 5. first-round agreement (Interpreter 1 round 1 used February 2024 imagery; redone on 2021 imagery)
def ld(p): return pd.read_csv(p,dtype=str).set_index('pid')
a1=ld(INTERP+'/random_points_interpreter1_round1.csv'); b=ld(INTERP+'/random_points_interpreter2.csv').loc[a1.index]; af=ld(INTERP+'/random_points_interpreter1.csv').loc[a1.index]
grp=lambda c:c.map({'1':'veg','2':'veg','3':'veg','4':'bare','5':'bare','6':'sealed','7':'rock','8':'water'})
def kk(a,b):
    ca,cb=a.cover.str[0],b.cover.str[0]
    return dict(cover8=round(kap(ca,cb),2),grouped=round(kap(grp(ca),grp(cb)),2),closed=round(kap(ca=='1',cb=='1'),2),sealed=round(kap(ca=='6',cb=='6'),2))
X['kappa_round1']=kk(a1,b); X['round1_img']=a1.img.value_counts().to_dict(); X['round1_changed_cover']=int((a1.cover!=af.cover).sum())
X['split_remaining']=dict(I1=int(len(ld(INTERP+'/random_points_interpreter1.csv'))-40),I2=int(len(ld(INTERP+'/random_points_interpreter2.csv'))-40))
R['extra']=X; json.dump(R,open(f'{OUT}/v4_results.json','w'),indent=1,ensure_ascii=False)
print(json.dumps(X,indent=1,ensure_ascii=False))
