import numpy as np, pandas as pd, geopandas as gpd, json, rasterio, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D
from scipy.stats import spearmanr
from common import *
os.makedirs(OUT+'/figs',exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.titlesize':9,'axes.titleweight':'bold','axes.spines.top':False,'axes.spines.right':False})
R=json.load(open(f'{OUT}/v4_results.json')); T=pd.read_csv(f'{OUT}/units_5ha_scored.csv').set_index('label')
U=gpd.read_file(f'{OUT}/units_5ha.shp').set_index('label').join(T[['CS','cls','P_top','Ae_DGcap','IC_C_sealed1']])
f=np.load(f'{OUT}/factors.npz'); ic=np.load(f'{OUT}/ic.npz')
ext=[tr.c,tr.c+tr.a*land.shape[1],tr.f+tr.e*land.shape[0],tr.f]
CLS=['#f7f4ea','#e9d8a6','#ee9b00','#ca6702','#9b2226']
def clean(a): a.set_xticks([]); a.set_yticks([]); a.set_xlabel(''); a.set_ylabel(''); [s.set_visible(False) for s in a.spines.values()]
def sbar(a,x0=None,y0=None):
    x0=x0 or ext[0]+150; y0=y0 or ext[2]+150; a.plot([x0,x0+1000],[y0,y0],'k',lw=2.5); a.text(x0+500,y0+80,'1 km',ha='center',fontsize=7)
def labels(a,col='k',fs=5.5):
    for lab,r in U.iterrows():
        p=r.geometry.representative_point(); a.text(p.x,p.y,lab,fontsize=fs,ha='center',va='center',color=col)
# ---------- Fig 2 framework
fig,a=plt.subplots(figsize=(7.4,4.4)); a.set_xlim(0,100); a.set_ylim(0,64); a.axis('off')
def box(x,y,w,h,t,fc,fs=6.8,bold=False):
    a.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.4,rounding_size=1.2',fc=fc,ec='#333',lw=.9))
    a.text(x+w/2,y+h/2,t,ha='center',va='center',fontsize=fs,weight='bold' if bold else 'normal',linespacing=1.25)
def arr(x1,y1,x2,y2): a.annotate('',xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle='-|>',color='#444',lw=.9))
box(0.5,46,20,15,'Open data\nCopernicus DEM,\nESA WorldCover,\nSentinel-2, soil (9 pts),\nstation rainfall','#eef2f7')
box(25,46,20,15,'Ranking units\n29 complete\nsea-outlet basins\n+ 12 coastal faces\n(whole island)','#eef2f7')
box(50,54.5,22,8,'Supply: RUSLE\n(Desmet–Govers LS)','#fff1d6')
box(50,44,22,8,'Transfer: Index of\nConnectivity (IC)','#fff1d6')
box(77,46,22.5,15,'Compound score\n+ Monte Carlo\n(K, C, LS form,\nIC weighting,\nblock weights)','#fde2c8')
arr(21,53.5,24.6,53.5); arr(45.4,56,49.6,58); arr(45.4,51,49.6,48); arr(72.4,58.5,76.6,55); arr(72.4,48,76.6,52)
a.plot([88,88],[45.5,41],color='#444',lw=.9); a.plot([11,88],[41,41],color='#444',lw=.9)
a.text(50,37.2,'Independent evaluation (no field access, no sediment data)',ha='center',fontsize=8,weight='bold')
ev=[('Sentinel-2 bare-soil\nfrequency\n2019–2024','objective'),('Nearshore turbidity\nafter rain\n(Sentinel-2 + IMERG)','objective'),('Land cover and\ndisturbance at\n193 random points\n(2 interpreters)','reproducible (κ)'),('Erosion-feature\ninventories\n(2 blind interpreters)','interpreter-dependent')]
cols={'objective':'#dcefd8','reproducible (κ)':'#dcefd8','interpreter-dependent':'#f6d5d5'}
for i,(t,k) in enumerate(ev):
    x=1+i*25; box(x,12,23,18,t,cols[k]); arr(x+11.5,41,x+11.5,30.6)
    a.text(x+11.5,8.5,k,ha='center',fontsize=6.8,style='italic',color='#555')
a.text(50,2,'Does the topographic ranking coincide with where disturbance and erosion are observed?',ha='center',fontsize=7.6)
plt.savefig(OUT+'/figs/Fig2_framework.png',dpi=300,bbox_inches='tight'); plt.close()
# ---------- Fig 3 supply / transfer
Ax=f['R']*f['K']*f['LS_DGcap']*f['C']; Ax[~land]=np.nan; icm=ic['IC_C_sealed1'].copy(); icm[~land]=np.nan
fig,ax=plt.subplots(1,4,figsize=(7.4,4.3),gridspec_kw=dict(width_ratios=[1,1,1,1.35]))
U.plot(ax=ax[0],color=np.where(U.kind=='basin','#a9c4e2','#f1c38e'),edgecolor='#555',lw=.4); labels(ax[0],fs=4.8); ax[0].set_title('(a) 41 ranking units',loc='left'); sbar(ax[0],ext[0]+100,ext[2]+100)
ax[0].legend(handles=[plt.Rectangle((0,0),1,1,fc='#a9c4e2',ec='#555'),plt.Rectangle((0,0),1,1,fc='#f1c38e',ec='#555')],labels=['complete basin','coastal face'],fontsize=6.3,loc='upper right',frameon=False)
im=ax[1].imshow(np.log10(np.maximum(Ax,.1)),extent=ext,cmap='YlOrRd',vmin=-1,vmax=3); U.boundary.plot(ax=ax[1],color='#333',lw=.3)
cb=plt.colorbar(im,ax=ax[1],orientation='horizontal',fraction=.04,pad=.02); cb.set_label('log$_{10}$ A (t ha$^{-1}$ yr$^{-1}$)',fontsize=7); ax[1].set_title('(b) Supply',loc='left')
im=ax[2].imshow(icm,extent=ext,cmap='viridis'); U.boundary.plot(ax=ax[2],color='w',lw=.3)
cb=plt.colorbar(im,ax=ax[2],orientation='horizontal',fraction=.04,pad=.02); cb.set_label('IC',fontsize=7); ax[2].set_title('(c) Transfer',loc='left')
for a in ax[:3]: clean(a)
c=np.where(T.kind=='basin','#3a6ea5','#d17a22'); ax[3].scatter(T.Ae_DGcap,T.IC_C_sealed1,c=c,s=22,edgecolor='k',lw=.4)
for lab,r in T.iterrows(): ax[3].annotate(lab,(r.Ae_DGcap,r.IC_C_sealed1),fontsize=5.5,xytext=(2,2),textcoords='offset points')
s=R['supply_vs_IC']['IC_C_sealed1']; ax[3].set_title(f"(d) ρ = {s['rho']:.2f} [{s['lo']:.2f}, {s['hi']:.2f}]",loc='left')
ax[3].set_xlabel('Supply (t ha$^{-1}$ yr$^{-1}$)'); ax[3].set_ylabel('Transfer: mean IC')
ax[3].legend(handles=[Line2D([],[],marker='o',ls='',mfc='#3a6ea5',mec='k',label='complete basin'),Line2D([],[],marker='o',ls='',mfc='#d17a22',mec='k',label='coastal face')],fontsize=6.5,frameon=False,loc='lower right')
plt.tight_layout(); plt.savefig(OUT+'/figs/Fig3_supply_transfer.png',dpi=300); plt.close()
# ---------- Fig 4 ranking
fig,ax=plt.subplots(1,3,figsize=(7.4,3.6),gridspec_kw=dict(width_ratios=[1,1.2,1.2]))
U.plot(ax=ax[0],column='P_top',cmap='magma_r',vmin=0,vmax=1,edgecolor='#444',lw=.3,legend=True,legend_kwds=dict(orientation='horizontal',fraction=.04,pad=.02,label='P(top quintile)'))
for lab_,r_ in U.iterrows():
    p_=r_.geometry.representative_point(); ax[0].text(p_.x,p_.y,lab_,fontsize=4.8,ha='center',va='center',color='w' if r_.P_top>.55 else 'k')
clean(ax[0]); ax[0].set_title('(a) Monte Carlo (n = 1000)',loc='left')
top=T.cls==5; ax[1].scatter(T.slope[~top],T.CS[~top],s=18,c='#bbb',edgecolor='k',lw=.3,label='other units'); ax[1].scatter(T.slope[top],T.CS[top],s=26,c='#9b2226',edgecolor='k',lw=.3,label='top quintile')
for lab,r in T[top].iterrows(): ax[1].annotate(lab,(r.slope,r.CS),fontsize=6,xytext=(2,2),textcoords='offset points')
ax[1].set_xlabel('Mean slope (°)'); ax[1].set_ylabel('Compound score'); ax[1].legend(fontsize=6.5,frameon=False)
ax[1].set_title(f"(b) Top 8 = 8 steepest units",loc='left')
lab=['Interbasins +\nmorphometry\n(previous)','Basins +\nmorphometry','Supply vs IC\nW = C, sealed 1','W = C,\nsealed 0.15','W = roughness','W = 1\n(topographic)']
vals=[(0.20,-0.17,0.52),(R['morph']['M_vs_supply']['rho'],R['morph']['M_vs_supply']['lo'],R['morph']['M_vs_supply']['hi'])]+[(R['supply_vs_IC'][k]['rho'],R['supply_vs_IC'][k]['lo'],R['supply_vs_IC'][k]['hi']) for k in ['IC_C_sealed1','IC_C_sealed0.15','IC_rough','IC_topo']]
y=np.arange(len(vals))[::-1]
for yi,(v,lo,hi),l in zip(y,vals,lab):
    col='#888' if 'previous' in l else '#3a6ea5'; ax[2].plot([lo,hi],[yi,yi],color=col,lw=1.5); ax[2].plot(v,yi,'o',color=col,ms=5)
ax[2].axvline(0,color='k',lw=.6,ls=':'); ax[2].set_yticks(y); ax[2].set_yticklabels(lab,fontsize=6.3); ax[2].set_xlabel('Spearman ρ (95% CI)')
ax[2].set_title('(c) Supply vs transfer',loc='left')
plt.tight_layout(); plt.savefig(OUT+'/figs/Fig4_ranking.png',dpi=300); plt.close()
# ---------- Fig 5 independent objective + reproducible evidence
tur=pd.read_csv(f'{OUT}/turbidity_unit_response.csv').set_index('label'); P=pd.read_csv(f'{OUT}/randompts_AB.csv')
fig,ax=plt.subplots(1,3,figsize=(7.4,2.9))
col=np.where(T.cls==5,'#9b2226','#bbb'); ax[0].scatter(T.CS,T.bare_nobuilt,c=col,s=18,edgecolor='k',lw=.3)
b=R['bare']['CS']; ax[0].set_title(f"(a) Bare soil (ρ = {b['rho']:.2f})",loc='left'); ax[0].set_xlabel('Compound score'); ax[0].set_ylabel('Mean bare frequency\n(built-up masked)')
col=np.where(tur.cls==5,'#9b2226','#bbb'); ax[1].scatter(tur.CS,tur.resp,c=col,s=18,edgecolor='k',lw=.3); ax[1].axhline(0,color='k',lw=.5,ls=':')
t=R['turbidity']['P3>=10']; ax[1].set_title(f"(b) Turbidity (ρ = {t['rho_CS']:.2f})",loc='left'); ax[1].set_xlabel('Compound score'); ax[1].set_ylabel('Δ ln turbidity\n(wet − dry, outlet − control)')
bt=R['randompts']['by_top']; names=['Closed\nforest','Dis-\nturbed','Bare\n≥10%','Sealed']; keys=['closed','D','bare10','sealed']
x=np.arange(4); ax[2].bar(x-.18,[bt[k]['top'] for k in keys],.36,color='#9b2226',label='top quintile'); ax[2].bar(x+.18,[bt[k]['rest'] for k in keys],.36,color='#bbb',label='other units')
kp=R['kappa']; kk=[kp['closed'],kp['disturbed'],kp['bare_w'],kp['sealed']]
for i,v in enumerate(kk): ax[2].text(i,max(bt[keys[i]]['top'],bt[keys[i]]['rest'])+.03,f'κ {v:.2f}',ha='center',fontsize=6.3)
ax[2].set_xticks(x); ax[2].set_xticklabels(names,fontsize=7); ax[2].set_ylabel('Fraction of random points'); ax[2].legend(fontsize=6.5,frameon=False); ax[2].set_title('(c) Random points',loc='left'); ax[2].set_ylim(0,1)
plt.tight_layout(); plt.savefig(OUT+'/figs/Fig5_evidence.png',dpi=300); plt.close()
# ---------- Fig 6 inventories
FA=gpd.read_file(INTERP+'/erosion_features_interpreter1.gpkg'); FB=gpd.read_file(INTERP+'/erosion_features_interpreter2.gpkg'); roc=pd.read_csv(f'{OUT}/roc_features_results.csv')
fig,ax=plt.subplots(1,2,figsize=(7.4,4.2),gridspec_kw=dict(width_ratios=[1,1.5]))
U.plot(ax=ax[0],column='cls',cmap=matplotlib.colors.ListedColormap(CLS),vmin=.5,vmax=5.5,edgecolor='#666',lw=.3)
FA.plot(ax=ax[0],marker='o',color='#1f77b4',markersize=14,edgecolor='k',lw=.4); FB.plot(ax=ax[0],marker='^',color='#2ca02c',markersize=16,edgecolor='k',lw=.4)
ax[0].legend(handles=[Line2D([],[],marker='o',ls='',mfc='#1f77b4',mec='k',label=f'Interpreter 1 (n = {len(FA)})'),Line2D([],[],marker='^',ls='',mfc='#2ca02c',mec='k',label=f'Interpreter 2 (n = {len(FB)})')]+[plt.Rectangle((0,0),1,1,fc=c,ec='#666') for c in CLS[::-1]],
   labels=[f'Interpreter 1 (n = {len(FA)})',f'Interpreter 2 (n = {len(FB)})','Priority 5 (top)','4','3','2','1 (lowest)'],fontsize=6,frameon=False,loc='center left',bbox_to_anchor=(0.62,0.32))
clean(ax[0]); ax[0].set_title('(a) Mapped erosion features',loc='left'); sbar(ax[0],ext[1]-1250)
maps=['slope only','supply A (pixel)','compound score (unit)','IC (pixel)','S2 bare frequency','near built-up (−distance)']
sets=[('A (17)','Interpreter 1','#1f77b4',-.2),('B (12)','Interpreter 2','#2ca02c',0),('A+B (29)','Both pooled','#555',.2)]
for s,lab,colr,off in sets:
    d=roc[roc.features==s].set_index('map').loc[maps]
    yy=np.arange(len(maps))[::-1]+off; ax[1].errorbar(d.AUC,yy,xerr=[d.AUC-d.lo,d.hi-d.AUC],fmt='o',color=colr,ms=4,lw=1,capsize=0,label=lab)
ax[1].axvline(.5,color='k',lw=.6,ls=':'); ax[1].set_yticks(np.arange(len(maps))[::-1]); ax[1].set_yticklabels(['Slope only','Supply (pixel)','Compound score','IC (pixel)','S2 bare frequency','Proximity to built-up'],fontsize=7)
ax[1].set_xlabel('AUC (95% bootstrap CI), features vs 189 random points'); ax[1].set_xlim(0,1); ax[1].legend(fontsize=6.5,frameon=False,loc='lower right')
ax[1].set_title('(b) Discrimination depends on the interpreter',loc='left')
plt.tight_layout(); plt.savefig(OUT+'/figs/Fig6_inventories.png',dpi=300); plt.close()
print('ok')
