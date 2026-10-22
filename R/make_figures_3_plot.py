import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'serif','font.size':8,'axes.linewidth':0.6})
b=pd.read_csv('bkgd.csv').values; c=pd.read_csv('clust.csv').values
gx=pd.read_csv('gx.csv').x.values; gy=pd.read_csv('gy.csv').y.values
ev=pd.read_csv('events.csv'); tg=ev[ev.flag==1]; db=tg[tg.mag>=6]
X,Y=np.meshgrid(gx,gy,indexing='ij'); inside=(Y>=38.35)&(Y<=38.50)&(X>=46.5)&(X<=47.0)
def amap(Z,label,fn,cmap,vmin=None,vmax=None):
    fig,ax=plt.subplots(figsize=(4.2,2.4),dpi=300)
    Zm=np.where(inside,Z,np.nan)
    im=ax.pcolormesh(X,Y,Zm,cmap=cmap,shading='auto',vmin=vmin,vmax=vmax)
    ax.add_patch(plt.Rectangle((46.5,38.35),0.5,0.15,fill=False,ec='k',lw=0.7))
    ax.scatter(tg.long,tg.lat,s=3,c='k',alpha=0.45,lw=0,label='M ≥ 3.2 (251)')
    ax.scatter(db.long,db.lat,marker='*',s=90,fc='none',ec='w',lw=1.4,zorder=5)
    ax.scatter(db.long,db.lat,marker='*',s=55,fc='yellow',ec='k',lw=0.6,zorder=6,label='Mw 6.5 and 6.3')
    ax.set_xlim(46.5,47.0); ax.set_ylim(38.30,38.55); ax.set_aspect(1/np.cos(np.deg2rad(38.42)))
    ax.set_xlabel('Longitude (°E)'); ax.set_ylabel('Latitude (°N)'); ax.legend(loc='upper left',fontsize=6.5,frameon=True,framealpha=.85)
    cb=fig.colorbar(im,ax=ax,fraction=0.046,pad=0.03); cb.set_label(label,fontsize=7.5); cb.ax.tick_params(labelsize=7)
    ax.set_title('(a)',fontsize=9,loc='left'); fig.savefig(fn,dpi=300,bbox_inches='tight',pad_inches=0.03); plt.close(fig)
amap(b,'background rate μ·u(x,y)\n(events per day per deg²)','fig4a.png','viridis')
amap(c,'clustering coefficient','fig5a.png','magma')
# Fig 8
t=pd.read_csv('tau.csv'); tt=t.tt.values; tau=t.tau.values; n=len(tt); N=np.arange(1,n+1); k=224
fig,ax=plt.subplots(2,1,figsize=(3.4,5.8),dpi=300)
a=ax[0]; a.step(tt,N,where='post',c='k',lw=1,label='observed'); a.plot(tt,tau,'--',c='tab:blue',lw=1.2,label='ETAS model (stationary)')
a.axvline(tt[k-1],c='tab:red',ls=':',lw=1); a.set_xlabel('ordinary time (days)'); a.set_ylabel('cumulative number'); a.legend(fontsize=7,loc='lower right'); a.set_title('(a)',fontsize=9,loc='left')
a.text(tt[k-1]+60,25,'day 217.7',color='tab:red',fontsize=7,rotation=90)
a=ax[1]; a.step(tau,N,where='post',c='k',lw=1,label='observed'); a.plot([0,n],[0,n],c='tab:red',lw=1,label='unit-rate Poisson'); a.axvline(tau[k-1],c='tab:red',ls=':',lw=1)
a.set_xlabel('transformed time τ'); a.set_ylabel('cumulative number'); a.legend(fontsize=7,loc='lower right'); a.set_title('(b)',fontsize=9,loc='left')
i=np.argmax(N-tau); a.annotate(f'max excess ≈ {N[i]-tau[i]:.0f}\n(event {i+1})',xy=(tau[i],N[i]),xytext=(tau[i]+25,N[i]-95),fontsize=7,arrowprops=dict(arrowstyle='->',lw=0.6))
fig.tight_layout(); fig.savefig('fig8.png',dpi=300); plt.close(fig); print('max excess',(N-tau).max(), i+1, tt[i])
# Fig 9
g=np.diff(tau); idx=np.arange(1,len(g)+1)
fig,ax=plt.subplots(figsize=(5.2,2.6),dpi=300)
ax.vlines(idx,0,g,color='0.55',lw=0.8); ax.plot(idx,g,'.',c='k',ms=2)
m1=g[:k].mean(); m2=g[k:].mean()
ax.hlines(m1,1,k,color='tab:blue',lw=1.6,label=f'mean {m1:.2f}'); ax.hlines(m2,k+1,len(g),color='tab:orange',lw=1.6,label=f'mean {m2:.2f}')
ax.axhline(1,c='tab:red',ls='--',lw=0.8,label='expected 1'); ax.axvline(k+0.5,c='k',ls='--',lw=0.9)
ax.set_xlabel('event index'); ax.set_ylabel('Δτ (transformed inter-event gap)'); ax.legend(fontsize=7,loc='upper left',frameon=False); ax.set_xlim(0,len(g)+1)
ax.text(k-3,ax.get_ylim()[1]*0.97,'break: event 224 (day 217.7)',fontsize=7,va='top',ha='right')
fig.tight_layout(); fig.savefig('fig9.png',dpi=300); plt.close(fig); print(m1,m2,g.max())
