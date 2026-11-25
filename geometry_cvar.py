"""Exploratory fixed/flexible arrival CVaR benchmark. Not flight software."""
import argparse, csv, json, math, hashlib
from datetime import date,timedelta
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from horizons_experiment import fetch, metrics, daynum, iso
from keplerx import lambert, MU, DAY
from multiwindow import WINDOWS
ROOT=Path(__file__).resolve().parent
C3=25.; DV=6.2; MIN_TOF=140; MAX_TOF=321; MAX_PENALTY=.1
DESIGN=np.arange(15); TEST=np.arange(15,31)

def cvar(x, alpha=.90):
    # Upper tail mean of discrete equally weighted outcomes; fractional boundary mass.
    x=np.sort(np.asarray(x,dtype=float))[::-1]; mass=(1-alpha)*len(x)
    weights=np.clip(mass-np.arange(len(x)),0,1)
    return float(x@weights/mass)

def write(path, rows):
    if rows:
        with path.open('w',newline='',encoding='utf8') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def experiment(year,earth,mars):
    first,last=map(daynum,WINDOWS[year]); cache={}; failures={}
    def calc(d,a):
        key=(d,a)
        if key not in cache and key not in failures:
            try:
                if not MIN_TOF<=a-d<=MAX_TOF:raise ValueError('flight duration out of range')
                x=metrics(earth[d],mars[a],a-d)
                if not np.isfinite(x).all():raise ValueError('nonfinite metrics')
                cache[key]=x
            except (ValueError,RuntimeError,KeyError,OverflowError,FloatingPointError) as ex:failures[key]=str(ex)
        return cache.get(key)
    def eligible(x):return x is not None and x[0]<=C3 and x[-1]<=DV
    def score(x):return (float(x[2][-1]),x[0],x[1])
    # Retain infeasible finite points; refine regardless of whether coarse search finds a viable point.
    coarse=[(d,d+t,calc(d,d+t)) for d in range(first,last+1,3) for t in sorted(set(range(MIN_TOF,MAX_TOF+1,5))|{MAX_TOF})]
    finite=[x for x in coarse if x[2] is not None]; cap=[x for x in finite if x[2][0]<=C3]
    if not finite:return {'year':year,'status':'no finite solutions'},[],[],[]
    seeds=sorted(cap or finite,key=score)[:3]; local=set()
    for d,a,_ in seeds:
        for dd in range(max(first,d-4),min(last,d+4)+1):
            for aa in range(a-6,a+7):
                if MIN_TOF<=aa-dd<=MAX_TOF:local.add((dd,aa))
    fine=[(d,a,calc(d,a)) for d,a in sorted(local)]; points=finite+[x for x in fine if x[2] is not None]
    capped=[x for x in points if x[2][0]<=C3]; valid=[x for x in capped if eligible(x[2])]
    def compact(x):return {'launch':iso(x[0]),'arrival':iso(x[1]),'c3':float(x[2][0]),'dv':float(x[2][-1])} if x else None
    info={'year':year,'coarse_evaluations':len(coarse),'fine_evaluations':len(fine),'min_dv_unconstrained':compact(min(points,key=score)),'min_dv_under_c3':compact(min(capped,key=score)) if capped else None}
    if not valid:return dict(info,status='no feasible point in scanned grid',solver_or_data_failures=len(failures)),[],[],[]
    d0,a0,_=min(valid,key=score)
    # Arrival slack is an explicit mission assumption, not an advantage over a same-arrival baseline.
    rows=[]; profiles={}; invalid=0
    for a in range(a0-7,a0+8):
        for d in range(max(first,d0-28,a-MAX_TOF),min(last,d0+14,a-MIN_TOF)+1):
            v=[calc(d+int(k),a) for k in DESIGN]
            if any(x is None for x in v):invalid+=1;continue
            v=np.array(v); ok=(v[:,0]<=C3)&(v[:,-1]<=DV)
            if not eligible(v[0]):continue
            r={'year':year,'launch':iso(d),'arrival':iso(a),'tof':a-d,'c3':float(v[0,0]),'nominal_dv':float(v[0,-1]),'design_cvar90':cvar(v[:,-1]),'design_feasible':float(ok.mean())}
            rows.append(r);profiles[(r['launch'],r['arrival'])]=v
    if not rows:return dict(info,status='no nominally feasible complete design profile'),[],[],[]
    fixed=[r for r in rows if r['arrival']==iso(a0)]
    nominal=min(fixed,key=lambda r:r['nominal_dv']) if fixed else min(rows,key=lambda r:r['nominal_dv'])
    # Compare objective and cost limits fairly; the nominal baseline is also refined on each domain.
    flexnom=min(rows,key=lambda r:r['nominal_dv'])
    def choose(pool,ref):
        q=[r for r in pool if r['design_feasible']==1 and r['nominal_dv']<=ref['nominal_dv']+MAX_PENALTY]
        return min(q,key=lambda r:(r['design_cvar90'],r['nominal_dv'],r['launch'],r['arrival'])) if q else None
    selected={'fixed_nominal':nominal,'fixed_cvar90':choose(fixed,nominal),'flex_nominal':flexnom,'flex_cvar90':choose(rows,flexnom)}
    for b in (7,14):
        target=date.fromisoformat(nominal['launch'])-timedelta(days=b)
        selected[f'fixed_buffer_{b}']=next((r for r in fixed if r['launch']==str(target)),None)
    records=[]; delays=[]
    for policy,r in selected.items():
        if r is None:continue
        d=daynum(date.fromisoformat(r['launch'])); a=daynum(date.fromisoformat(r['arrival']))
        outcomes=[]; max_inc=0.
        for k in range(31):
            v=calc(d+k,a); angle=None; inc=None
            if v is not None:
                rr=earth[d+k][:3]; ss=mars[a][:3]; cross=np.cross(rr,ss)
                angle=float(np.degrees(np.arctan2(np.linalg.norm(cross),rr@ss)))
                vel,_=lambert(rr,ss,(a-d-k)*DAY); h=np.cross(rr,vel)
                inc=float(np.degrees(np.arctan2(np.linalg.norm(h[:2]),h[2])))
                max_inc=max(max_inc,inc)
            ok=eligible(v);outcomes.append(ok)
            delays.append({'year':year,'policy':policy,'delay':k,'launch':iso(d+k),'arrival':iso(a),'c3':float(v[0]) if v is not None else None,'dv':float(v[-1]) if v is not None else None,'feasible':int(ok),'angle_deg':angle,'inclination_deg':inc})
        records.append(dict(r,policy=policy,stress_feasible=float(np.mean(outcomes[15:22])),extreme_feasible=float(np.mean(outcomes[22:])),uniform_feasible=float(np.mean(outcomes)),max_inclination_deg=max_inc))
    return dict(info,status='ok',reference_arrival=iso(a0),candidate_count=len(rows),incomplete_design_profiles=invalid,solver_or_data_failures=len(failures),missing_policies=[k for k,v in selected.items() if v is None],policies=records),rows,records,delays

def main():
    p=argparse.ArgumentParser();p.add_argument('--years',nargs='+',type=int,choices=list(WINDOWS),default=list(WINDOWS));p.add_argument('--offline',action='store_true');args=p.parse_args()
    out=ROOT/'results_v05';out.mkdir(exist_ok=True);raw=ROOT/'horizons_raw';raw.mkdir(exist_ok=True)
    summary={'model':'Sun-only prograde zero-revolution Lambert; ideal impulsive patched conics','method':'Exploratory analysis after seeing v0.4 outcomes; discrete uniform CVaR90 design days 0-14; evaluation days 15-30; arrival slack +/-7 days; max nominal penalty 0.1 km/s relative to each domain nominal','limits':{'c3':C3,'dv':DV,'tof':[MIN_TOF,MAX_TOF]},'sources':[],'windows':[]};allrows=[];allpol=[];alld=[]
    for y in dict.fromkeys(args.years):
        start,stop=(date(2028,1,1),date(2029,11,20)) if y==2028 else (WINDOWS[y][0],WINDOWS[y][1]+timedelta(days=MAX_TOF))
        if args.offline and not all((raw/f'{i}_{start}_{stop}.json').exists() for i in (399,499)):raise FileNotFoundError(f'Missing cache for {y}')
        e,es=fetch(399,start,stop,raw);m,ms=fetch(499,start,stop,raw);summary['sources'] += [es,ms]
        print(f'Analyzing {y}...',flush=True);s,r,p,d=experiment(y,e,m);summary['windows'].append(s);allrows+=r;allpol+=p;alld+=d
        print(json.dumps(s,indent=2),flush=True);(out/'summary_v05.json').write_text(json.dumps(summary,indent=2),encoding='utf8');write(out/'candidates.csv',allrows);write(out/'policies.csv',allpol);write(out/'geometry.csv',alld)
    print('Saved results_v05/summary_v05.json',flush=True)
if __name__=='__main__':
    assert abs(cvar(np.arange(15))-(14+13*.5)/1.5)<1e-10
    main()
