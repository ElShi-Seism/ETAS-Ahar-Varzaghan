# Completeness check: magnitude-frequency in time periods and before/after the break. Run: python3 R/completeness_check.py (needs ascii_correct.txt). GPL-3.0-or-later.
import pandas as pd,numpy as np
from scipy import stats
d=pd.read_csv('ascii_correct.txt',sep='\t',header=None,names="long lat year month day mag sec hour minute".split())
w=d[(d.lat.between(38.35,38.5))&(d.long.between(46.5,47))].copy()
w['t']=pd.to_datetime(dict(year=w.year,month=w.month,day=w.day,hour=w.hour,minute=w.minute))
t0=pd.Timestamp('2012-08-11'); w=w[w.t>=t0]; w['day']=(w.t-t0).dt.total_seconds()/86400
w=w.sort_values('day'); print(len(w))
def mc_maxc(m):
    bins=np.round(np.arange(2.4,m.max()+0.2,0.1),1); h,_=np.histogram(m,bins=np.append(bins,bins[-1]+0.1)-0.05)
    return bins[np.argmax(h)]
def bmle(m,mc): m=m[m>=mc]; return 1/(np.log(10)*(m.mean()-(mc-0.05))),len(m)
# periods
cuts=[(0,10),(10,100),(100,217.73),(217.73,500),(500,1500),(1500,4490)]
rows=[]
for a,b in cuts:
    s=w[(w.day>=a)&(w.day<b)]; m=s.mag.values
    mc=mc_maxc(m) if len(m)>=30 else np.nan
    bb=bmle(m,3.2) if (m>=3.2).sum()>=10 else (np.nan,0)
    rows.append(dict(a=a,b=b,n25=len(m),n32=int((m>=3.2).sum()),frac_lt32=round((m<3.2).mean(),3),Mc_maxc=mc,b32=round(bb[0],2),mmin=m.min() if len(m) else np.nan))
print(pd.DataFrame(rows).to_string())
# before/after break
pre=w[w.day<217.73]; post=w[w.day>=217.73]
print('M2.5-3.2 share pre/post:',(pre.mag<3.2).mean().round(3),(post.mag<3.2).mean().round(3),len(pre),len(post))
tab=[[ (pre.mag<3.2).sum(),(pre.mag>=3.2).sum()],[(post.mag<3.2).sum(),(post.mag>=3.2).sum()]]
print(tab, stats.fisher_exact(tab))
# the same restricted to day>=30 to avoid early incompleteness
pre2=w[(w.day>=30)&(w.day<217.73)]
tab2=[[ (pre2.mag<3.2).sum(),(pre2.mag>=3.2).sum()],[(post.mag<3.2).sum(),(post.mag>=3.2).sum()]]
print('pre from day30',tab2, stats.fisher_exact(tab2))
# magnitude distribution of M>=3.2 pre vs post
a=pre[pre.mag>=3.2].mag; b=post[post.mag>=3.2].mag
print('KS mags>=3.2',stats.ks_2samp(a,b), 'b pre',bmle(a.values,3.2),'b post',bmle(b.values,3.2))
# rolling: share of M<3.2 per 25-event windows of events (ordered)
ev=(w.mag<3.2).astype(int).values; 
import numpy as np
for k in range(0,len(ev),150): print(int(w.day.values[k]), round(ev[k:k+150].mean(),2))
# yearly Mc estimates (max curvature) using ≥2.5
for y in range(2012,2025):
    m=w[w.year==y].mag.values
    if len(m)>=25: print(y,len(m),mc_maxc(m),round(np.median(m),2))
w[['day','mag','long','lat']].to_csv('results/events_2p5.csv',index=False)
