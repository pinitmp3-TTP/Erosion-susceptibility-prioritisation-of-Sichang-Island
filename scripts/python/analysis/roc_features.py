import numpy as np, pandas as pd, geopandas as gpd, rasterio, sys, os
from scipy import ndimage
from sklearn.metrics import roc_auc_score
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')); from common import *
f=np.load(f'{OUT}/factors.npz'); ic=np.load(f'{OUT}/ic.npz'); cf=np.load(f'{OUT}/cover_frac.npz')
unit=np.load(f'{OUT}/unit_5ha.npy'); T=pd.read_csv(f'{OUT}/units_5ha_mc.csv').set_index('uid')
cs=np.full(land.shape,np.nan)
for u,v in T.CS.items(): cs[unit==u]=v
ptop=np.full(land.shape,np.nan)
for u,v in T.P_top.items(): ptop[unit==u]=v
built=(cf['built']>=0.25)&land
dbuilt=ndimage.distance_transform_edt(~built)*cell
maps={'compound score (unit)':cs,'supply A (pixel)':f['R']*f['K']*f['LS_DGcap']*f['C'],'IC (pixel)':ic['IC_C_sealed1'],
      'slope only':f['slope'],'near built-up (−distance)':-dbuilt}
with rasterio.open(BARE_FILE) as r: bare_src=r
def sample(gx,gy):
    rr,cc=rasterio.transform.rowcol(tr,gx,gy); rr=np.clip(np.array(rr),0,land.shape[0]-1); cc=np.clip(np.array(cc),0,land.shape[1]-1)
    d={k:m[rr,cc] for k,m in maps.items()}
    with rasterio.open(BARE_FILE) as r: d['S2 bare frequency']=np.array([v[BARE_BAND-1] for v in r.sample(zip(gx,gy))])
    d['on_land']=land[rr,cc]; d['label']=[None]*len(gx)
    return pd.DataFrame(d)
# negatives: interpreted random points with no erosion, not water
cols=['pid','lon','lat','img','cover','bare','dist','eros','where','old_img','seen','conf','notes']
def pts(p):   # A = Interpreter 1, B = Interpreter 2
    d=pd.read_csv(p,dtype=str); d[['lon','lat']]=d[['lon','lat']].astype(float); return d[cols]
PA=pts(INTERP+'/random_points_interpreter1.csv'); PB=pts(INTERP+'/random_points_interpreter2.csv')
P=pd.concat([PA.assign(who='A'),PB.assign(who='B')]); P=P[P.cover.notna()]
neg=P[(P.eros=='ไม่มี')&~P.cover.astype(str).str.startswith('8')].drop_duplicates('pid')
key=pd.read_csv(INTERP+'/random_points_key.csv').set_index('pid')
neg=neg.join(key[['x','y']],on='pid')
N=sample(neg.x.values,neg.y.values); N=N[N.on_land]
FA=gpd.read_file(INTERP+'/erosion_features_interpreter1.gpkg').assign(who='A'); FB=gpd.read_file(INTERP+'/erosion_features_interpreter2.gpkg').assign(who='B')
F=pd.concat([FA,FB]); F=gpd.GeoDataFrame(F,geometry='geometry',crs=FA.crs)
Fs=sample(F.geometry.x.values,F.geometry.y.values); Fs['who']=F.who.values; Fs['new']=(F.seen.values=='ภาพเก่ายังไม่มี'); Fs['conf']=F.conf.values
print(f'negatives (random points, no erosion, on land): {len(N)}; features A {len(FA)}, B {len(FB)}')
rng=np.random.default_rng(3)
def auc_ci(pos,neg,nb=2000):
    y=np.r_[np.ones(len(pos)),np.zeros(len(neg))]; s=np.r_[pos,neg]; ok=np.isfinite(s); a=roc_auc_score(y[ok],s[ok])
    bs=[]
    for _ in range(nb):
        p=rng.choice(pos,len(pos)); q=rng.choice(neg,len(neg)); yy=np.r_[np.ones(len(p)),np.zeros(len(q))]; ss=np.r_[p,q]; o=np.isfinite(ss)
        bs.append(roc_auc_score(yy[o],ss[o]))
    return a,*np.percentile(bs,[2.5,97.5])
sets={'A (17)':Fs[Fs.who=='A'],'B (12)':Fs[Fs.who=='B'],'A+B (29)':Fs,'A+B new since 2013/2017':Fs[Fs.new]}
keys=[k for k in maps]+['S2 bare frequency']
rows=[]
for sname,S in sets.items():
    for k in keys:
        a,lo,hi=auc_ci(S[k].values,N[k].values); rows.append((sname,len(S),k,a,lo,hi))
R=pd.DataFrame(rows,columns=['features','n','map','AUC','lo','hi'])
pd.set_option('display.width',200)
print(R.pivot_table(index='map',columns='features',values='AUC').round(2).to_string())
print(R.assign(ci=R.apply(lambda r:f'{r.AUC:.2f} [{r.lo:.2f}, {r.hi:.2f}]',axis=1)).pivot(index='map',columns='features',values='ci').to_string())
R.to_csv(f'{OUT}/roc_features_results.csv',index=False)
# where are the features: priority class
Fs['cls']=[T.cls.get(u) for u in unit[rasterio.transform.rowcol(tr,F.geometry.x.values,F.geometry.y.values)]]
print('\nfeatures by priority class (5=top):'); print(pd.crosstab(Fs.who,Fs.cls))
print('features within 60 m of built-up: A %d/17, B %d/12'%((Fs[Fs.who=="A"]['near built-up (−distance)']>=-60).sum(),(Fs[Fs.who=="B"]['near built-up (−distance)']>=-60).sum()))
print('negatives within 60 m of built-up: %.0f%%'%(100*(N['near built-up (−distance)']>=-60).mean()))
