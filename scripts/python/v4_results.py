import numpy as np, pandas as pd, json, openpyxl
from scipy.stats import spearmanr, pearsonr, mannwhitneyu
from sklearn.metrics import cohen_kappa_score as kap
from common import *
rng=np.random.default_rng(1); RES={}
z=lambda x:(x-x.mean())/x.std(ddof=1)
qcls=lambda s: pd.qcut(s.rank(method='first'),5,labels=[1,2,3,4,5]).astype(int)
def rci(x,y,nperm=5000):
    x=np.asarray(x,float); y=np.asarray(y,float); ok=np.isfinite(x)&np.isfinite(y); x,y=x[ok],y[ok]
    r=spearmanr(x,y)[0]; n=len(x); zf=np.arctanh(r); se=1.06/np.sqrt(n-3)
    p=(np.abs([spearmanr(x,rng.permutation(y))[0] for _ in range(nperm)])>=abs(r)).mean()
    return dict(rho=round(r,2),lo=round(np.tanh(zf-1.96*se),2),hi=round(np.tanh(zf+1.96*se),2),p=round(p,3),n=n)
f=np.load(f'{OUT}/factors.npz'); sl=f['slope']
# ---- island-level soil loss by LS form
for k in ['MB','DG','DGcap']:
    Ax=f['R']*f['K']*f[f'LS_{k}']*f['C']; Ax[~land]=np.nan
    RES[f'A_island_{k}']=dict(mean=round(np.nanmean(Ax),1),median=round(np.nanmedian(Ax),1),tol_pct=round(100*np.mean(Ax[land]>12.5),1),
        share_from_gt35=round(100*np.nansum(Ax[land&(sl>np.degrees(np.arctan(.35)))])/np.nansum(Ax),1))
RES['steep_area_pct']=round(100*np.mean(sl[land]>np.degrees(np.arctan(.35))),1)
# ---- units
T=pd.read_csv(f'{OUT}/units_5ha_table.csv').set_index('label'); MC=pd.read_csv(f'{OUT}/units_5ha_mc.csv').set_index('label')
T=T.join(MC[['P_top','rank_med','rank_p5','rank_p95']])
RES['units']=dict(n=len(T),basins=int((T.kind=='basin').sum()),faces=int((T.kind=='face').sum()),
  basin_ha=round(T[T.kind=='basin'].area_ha.sum(),1),face_ha=round(T[T.kind=='face'].area_ha.sum(),1),
  basin_ha_range=[round(T[T.kind=='basin'].area_ha.min(),1),round(T[T.kind=='basin'].area_ha.max(),1)],face_ha_range=[round(T[T.kind=='face'].area_ha.min(),1),round(T[T.kind=='face'].area_ha.max(),1)])
SUP,TR='Ae_DGcap','IC_C_sealed1'
RES['supply_unit_range']=[round(T[SUP].min(),1),round(T[SUP].max(),1)]
RES['supply_vs_IC']={k:rci(T[SUP],T[k]) for k in ['IC_C_sealed1','IC_C_sealed0.15','IC_rough','IC_topo']}
RES['supply_variants_vs_IC']={k:rci(T[k],T[TR],2000) for k in ['A_MB','Ae_MB','Ae_DG','Ae_DGcap_ndvi']}
B=T[T.kind=='basin'].copy(); MORPH=['Dd','H','Rh','HI','Rn','Ff','Rc']
Zm=B[MORPH].apply(z); Zm[['Ff','Rc']]*=-1; B['M']=Zm.mean(1)
RES['morph']=dict(M_vs_IC=rci(B.M,B[TR],2000),M_vs_supply=rci(B.M,B[SUP],2000),r_H_Rh=round(pearsonr(B.H,B.Rh)[0],2),r_Ff_Rc=round(pearsonr(B.Ff,B.Rc)[0],2),n=len(B))
T['zS']=z(T[SUP]); T['zT']=z(T[TR]); T['CS']=0.5*T.zS+0.5*T.zT; T['cls']=qcls(T.CS)
T['driver']=np.where((T.zS>=.5)&(T.zT>=.5),'supply+transfer',np.where(T.zS>T.zT,'supply','transfer'))
top=T[T.cls==5].sort_values('CS',ascending=False)
RES['top']=dict(units=top.index.tolist(),ha=round(top.area_ha.sum(),1),pct=round(100*top.area_ha.sum()/T.area_ha.sum(),1),
  drivers=top.driver.value_counts().to_dict(),P_top_min=round(top.P_top.min(),2),next_P_top=round(T[T.cls<5].P_top.max(),2),
  high=T[T.cls==4].sort_values('CS',ascending=False).index.tolist())
RES['top_eq_slope_top8']=len(set(T.slope.nlargest(len(top)).index)&set(top.index))
RES['CS_vs_slope']=round(spearmanr(T.CS,T.slope)[0],2); RES['supply_vs_slope']=round(spearmanr(T[SUP],T.slope)[0],2); RES['IC_vs_slope']=round(spearmanr(T[TR],T.slope)[0],2)
k=len(top); ts=set(top.index)
obsS=len(ts&set(T[qcls(T[SUP])==5].index)); obsT=len(ts&set(T[qcls(T[TR])==5].index))
rp=2*np.sin(np.pi*spearmanr(T[SUP],T[TR])[0]/6); sims=[]
for _ in range(20000):
    x=rng.multivariate_normal([0,0],[[1,rp],[rp,1]],len(T)); c=0.5*z(pd.Series(x[:,0]))+0.5*z(pd.Series(x[:,1]))
    tc=set(np.argsort(-c.values)[:k]); sims.append((len(tc&set(np.argsort(-x[:,0])[:k])),len(tc&set(np.argsort(-x[:,1])[:k]))))
sims=np.array(sims); RES['single_block']=dict(k=k,supply_only=obsS,transfer_only=obsT,null_mean=sims.mean(0).round(1).tolist(),null_95=np.percentile(sims,[2.5,97.5],axis=0).T.tolist())
RES['faces_in_top']=int((top.kind=='face').sum())
RES['bare']={c:rci(T[b],T[c],2000) for b in ['bare_nobuilt'] for c in ['CS',SUP,TR]}
RES['bare_all_CS']=rci(T.bare_all,T.CS,2000)
# thresholds
th={}
for t in ['2','10']:
    X=pd.read_csv(f'{OUT}/units_{t}ha_table.csv'); th[t]=dict(n=len(X),**{k:rci(X[SUP],X[k],1000)['rho'] for k in ['IC_C_sealed1','IC_topo']})
RES['thresholds']=th
# ---- turbidity
d=pd.read_csv(D+'/Sichang_turbidity_outlets.csv'); d=d[d.T_median>0]
key=pd.read_csv(INTERP+'/random_points_key.csv')[['ucode','label']].drop_duplicates()
w=d.pivot_table(index=['date','ucode'],columns='kind',values='T_median').dropna().reset_index()
w=w.merge(d.groupby('date')[['P1','P3','P7']].first(),left_on='date',right_index=True).merge(key,on='ucode'); w['dT']=np.log(w.outlet)-np.log(w.control)
tur={'images':int(d.date.nunique()),'paired_dates':int(w.date.nunique())}
for name,wet in [('P3>=10',w.P3>=10),('P3>=5',w.P3>=5),('P3>=20',w.P3>=20),('P1>=5',w.P1>=5)]:
    g=w.assign(st=np.where(wet,'wet',np.where(w.P7<1,'dry','mid'))); g=g[g.st!='mid']
    Rr=g.groupby(['label','st']).dT.median().unstack(); n=g[g.st=='wet'].groupby('label').size(); Rr=Rr[n.reindex(Rr.index).fillna(0)>=3]
    Rr['resp']=Rr.wet-Rr.dry; Rr=Rr.join(T[['CS','cls']])
    tur[name]=dict(units=len(Rr),wet_dates=int(g[g.st=='wet'].date.nunique()),dry_dates=int(g[g.st=='dry'].date.nunique()),rho_CS=round(spearmanr(Rr.resp,Rr.CS)[0],2),p=round(spearmanr(Rr.resp,Rr.CS)[1],3),
        top_med=round(Rr[Rr.cls==5].resp.median(),2),rest_med=round(Rr[Rr.cls<5].resp.median(),2),mw_p_top_gt=round(mannwhitneyu(Rr[Rr.cls==5].resp,Rr[Rr.cls<5].resp,alternative='greater').pvalue,2))
    if name=='P3>=10':
        sig=[]
        for lab,gg in g.groupby('label'):
            a=gg[gg.st=='wet'].dT.values; b=gg[gg.st=='dry'].dT.values
            if len(a)>=3 and len(b)>=3:
                bs=[np.median(rng.choice(a,len(a)))-np.median(rng.choice(b,len(b))) for _ in range(2000)]
                if np.percentile(bs,2.5)>0: sig.append(lab)
        tur[name]['sig_units']=sig; tur[name]['sig_units_cls']=[int(T.cls[s]) for s in sig]
        Rr.to_csv(f'{OUT}/turbidity_unit_response.csv')
RES['turbidity']=tur
# ---- random points (A+B)
cols=['pid','lon','lat','img','cover','bare','dist','eros','where','old_img','seen','conf','notes']
def pts(p):   # interpretation tables (parsed from the original forms); A = Interpreter 1, B = Interpreter 2
    d=pd.read_csv(p,dtype=str); d[['lon','lat']]=d[['lon','lat']].astype(float); return d[cols]
A=pts(INTERP+'/random_points_interpreter1.csv').assign(who='A'); Bp=pts(INTERP+'/random_points_interpreter2.csv').assign(who='B')
ov=sorted(set(A.pid[:40])); a=A.set_index('pid').loc[ov]; b=Bp.set_index('pid').loc[ov]
ca,cb=a.cover.str[0],b.cover.str[0]; grp=lambda c:c.map({'1':'veg','2':'veg','3':'veg','4':'bare','5':'bare','6':'sealed','7':'rock','8':'water'})
o={'0':0,'น้อยกว่า 10%':1,'10–50%':2,'มากกว่า 50%':3}
RES['kappa']=dict(cover8=round(kap(ca,cb),2),cover_grouped=round(kap(grp(ca),grp(cb)),2),closed=round(kap(ca=='1',cb=='1'),2),sealed=round(kap(ca=='6',cb=='6'),2),
   disturbed=round(kap(a.dist!='ไม่มี',b.dist!='ไม่มี'),2),bare_w=round(kap(a.bare.astype(str).map(o),b.bare.astype(str).map(o),weights='linear'),2),n=len(ov))
P=pd.concat([A,Bp],ignore_index=True); P=P[P.cover.notna()&~P.cover.astype(str).str.startswith('8')]
P=P.join(pd.read_csv(INTERP+'/random_points_key.csv').set_index('pid')[['label']],on='pid')
P['D']=(P.dist!='ไม่มี').astype(int); P['bare10']=P.bare.astype(str).isin(['10–50%','มากกว่า 50%']).astype(int)
P['sealed']=P.cover.astype(str).str.startswith('6').astype(int); P['closed']=P.cover.astype(str).str.startswith('1').astype(int)
P['E']=(~P.eros.isin(['ไม่มี','อื่น ๆ'])&P.eros.notna()).astype(int)
P['cls']=P.label.map(T.cls)
RES['randompts']=dict(n_obs=len(P),n_points=P.pid.nunique(),E_points=int(P.E.sum()),
   by_top={c:dict(top=round(P[P.cls==5][c].mean(),2),rest=round(P[P.cls<5][c].mean(),2)) for c in ['closed','sealed','D','bare10']})
g=P.groupby('label')[['D','bare10','sealed','closed']].mean().join(T.CS); RES['randompts_unit_rho']={c:rci(g[c],g.CS,2000) for c in ['D','bare10','closed']}
RES['randompts_fisher']={}
from scipy.stats import fisher_exact
for c in ['closed','D']:
    tab=pd.crosstab((P.cls==5).values,P[c].values); RES['randompts_fisher'][c]=round(fisher_exact(tab.values)[1],4)
json.dump(RES,open(f'{OUT}/v4_results.json','w'),indent=1,default=lambda o:o.item() if hasattr(o,'item') else str(o))
T.to_csv(f'{OUT}/units_5ha_scored.csv'); P.to_csv(f'{OUT}/randompts_AB.csv',index=False)
print(json.dumps(RES,indent=1,default=str)[:6000])
