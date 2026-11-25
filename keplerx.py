"""Educational preliminary design; km, seconds, days. J2000 ecliptic."""
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from datetime import date,timedelta
import csv,json
from pathlib import Path
MU=1.32712440018e11; AU=149597870.7; DAY=86400.; EPOCH=date(2000,1,1)
ELEMENTS={"earth":([1.00000261,.01671123,-.00001531,100.46457166,102.93768193,0],[.00000562,-.00004392,-.01294668,35999.37244981,.32327364,0]),"mars":([1.52371034,.09339410,1.84969142,-4.55343205,-23.94362959,49.55953891],[.00001847,.00007882,-.00813131,19140.30268499,.44441088,-.29257343])}
def pos(planet,d):
    base,rate=ELEMENTS[planet]; a,e,i,L,p,n=np.array(base)+np.array(rate)*(d-.5)/36525
    i,L,p,n=np.deg2rad([i,L,p,n]); m=(L-p+np.pi)%(2*np.pi)-np.pi
    E=brentq(lambda x:x-e*np.sin(x)-m,-np.pi,np.pi); w=p-n
    x=a*AU*(np.cos(E)-e); y=a*AU*np.sqrt(1-e*e)*np.sin(E)
    c,s=np.cos,np.sin
    return np.array([(c(w)*c(n)-s(w)*s(n)*c(i))*x+(-s(w)*c(n)-c(w)*s(n)*c(i))*y,(c(w)*s(n)+s(w)*c(n)*c(i))*x+(-s(w)*s(n)+c(w)*c(n)*c(i))*y,s(w)*s(i)*x+c(w)*s(i)*y])
def state(p,d): return pos(p,d),(pos(p,d+.01)-pos(p,d-.01))/(.02*DAY)
def stump(z):
    if abs(z)<1e-5:return .5-z/24+z*z/720,1/6-z/120+z*z/5040
    if z>0:
        q=np.sqrt(z);return (1-np.cos(q))/z,(q-np.sin(q))/q**3
    q=np.sqrt(-z);return (np.cosh(q)-1)/(-z),(np.sinh(q)-q)/q**3
def lambert(r1,r2,dt):
    r1=np.array(r1);r2=np.array(r2);a,b=np.linalg.norm(r1),np.linalg.norm(r2)
    cos=np.clip(r1@r2/a/b,-1,1); sign=1 if np.cross(r1,r2)[2]>=0 else -1
    A=sign*np.sqrt(a*b*(1+cos))
    if abs(A)<1e-5:raise ValueError('degenerate geometry')
    def ytime(z):
        C,S=stump(z); y=a+b+A*(z*S-1)/np.sqrt(C)
        if y<=0:return y,-dt
        return y,((y/C)**1.5*S+A*np.sqrt(y))/np.sqrt(MU)-dt
    lo=-1.;hi=4*np.pi**2-1e-4
    while ytime(lo)[1]>0 and lo>-1000:lo*=2
    z=brentq(lambda z:ytime(z)[1],lo,hi,xtol=1e-11); y,_=ytime(z)
    f=1-y/a;g=A*np.sqrt(y/MU);gd=1-y/b
    return (r2-f*r1)/g,(gd*r2-r1)/g

def metrics(d,t,offset=None):
    r1,u1=state('earth',d);r2,u2=state('mars',d+t)
    if offset is not None:r1=r1+offset[:3];r2=r2+offset[3:]
    v1,v2=lambert(r1,r2,t*DAY); vin=np.linalg.norm(v1-u1);vout=np.linalg.norm(v2-u2)
    # impulsive departure from 300 km circular Earth orbit and capture into 300 km circular Mars orbit
    dep=np.sqrt(vin**2+2*398600.435/6678.137)-np.sqrt(398600.435/6678.137)
    cap=np.sqrt(vout**2+2*42828.375/3696.19)-np.sqrt(42828.375/3696.19)
    return np.array([vin**2,vout,dep,cap,dep+cap])

def run():
    out=Path(__file__).parent; start=(date(2028,1,1)-EPOCH).days
    rng=np.random.default_rng(42); jitter=rng.normal(0,7,128); offsets=rng.normal(size=(128,6))*np.array([6000]*3+[25000]*3)
    rows=[]
    for d in range(start,start+366,7):
        for t in range(140,341,10):
            try:
                m=metrics(d,t); rows.append([d,t,*m])
            except (ValueError,RuntimeError):pass
    rows=np.array(rows); candidates=rows[(rows[:,2]<=25)&(rows[:,1]<=300)]
    if not len(candidates):raise RuntimeError('no feasible grid points')
    # evaluate common random scenarios on the 45 lowest nominal delta-v candidates
    candidates=candidates[np.argsort(candidates[:,-1])[:45]]; scored=[]
    for row in candidates:
        d,t=row[:2]; sims=np.array([metrics(d+j,t,o) for j,o in zip(jitter,offsets)])
        dv=sims[:,-1]; scored.append([*row,dv.mean(),np.quantile(dv,.95),np.mean(sims[:,0]>25)])
    scored=np.array(scored); nominal=scored[0]; robust=scored[np.argmin(scored[:,-2])]
    def record(x):return dict(launch=str(EPOCH+timedelta(days=int(x[0]))),tof_days=int(x[1]),c3=float(x[2]),arrival_vinf=float(x[3]),departure_dv=float(x[4]),capture_dv=float(x[5]),total_dv=float(x[6]),mean_dv=float(x[7]),p95_dv=float(x[8]),c3_exceedance=float(x[9]))
    # independent numerical integration of the nominal Lambert arc
    d,t=nominal[:2]; r1,_=state('earth',d);r2,_=state('mars',d+t);v1,v2=lambert(r1,r2,t*DAY)
    def ode(_,y):return np.r_[y[3:],-MU*y[:3]/np.linalg.norm(y[:3])**3]
    sol=solve_ivp(ode,[0,t*DAY],np.r_[r1,v1],rtol=1e-11,atol=1e-9); residual=float(np.linalg.norm(sol.y[:3,-1]-r2))
    # Pareto set across nominal cost, duration, and 95th percentile cost, restricted to shortlist
    objectives=scored[:,[6,1,8]]; mask=np.array([not np.any(np.all(objectives<=x,axis=1)&np.any(objectives<x,axis=1)) for x in objectives])
    summary=dict(grid_points=len(rows),shortlist=len(scored),scenarios_per_candidate=len(jitter),nominal=record(nominal),robust=record(robust),integration_endpoint_error_km=residual,pareto_count=int(mask.sum()))
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
    np.savetxt(out/'grid.csv',rows,delimiter=',',header='day_since_2000_01_01,tof_days,c3,arrival_vinf,departure_dv,capture_dv,total_dv',comments='')
    np.savetxt(out/'robustness.csv',scored,delimiter=',',header='day_since_2000_01_01,tof_days,c3,arrival_vinf,departure_dv,capture_dv,total_dv,mean_dv,p95_dv,c3_exceedance',comments='')
    np.savetxt(out/'pareto.csv',scored[mask],delimiter=',',header='day_since_2000_01_01,tof_days,c3,arrival_vinf,departure_dv,capture_dv,total_dv,mean_dv,p95_dv,c3_exceedance',comments='')
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(10,6));dates=[EPOCH+timedelta(days=int(x)) for x in rows[:,0]]; c=ax.tricontourf(matplotlib.dates.date2num(dates),rows[:,1],np.clip(rows[:,6],0,15),levels=25,cmap='viridis_r');ax.xaxis_date();fig.autofmt_xdate();ax.set(ylabel='Flight duration (days)',xlabel='Launch date (approximate)',title='KEPLER-X | Earth–Mars total impulsive delta-v');fig.colorbar(c,ax=ax,label='km/s (clipped at 15)');
    for x,l,mk in [(nominal,'Nominal','o'),(robust,'Robust','x')]:ax.scatter(EPOCH+timedelta(days=int(x[0])),x[1],marker=mk,s=100,label=l,color='red')
    ax.legend();fig.tight_layout();fig.savefig(out/'porkchop.png',dpi=180);plt.close(fig)
if __name__=='__main__':run()
