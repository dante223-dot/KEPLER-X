"""Preliminary Earth-Mars fixed-arrival delay benchmark; not flight navigation."""
import argparse, json, csv, math, platform
from pathlib import Path
from datetime import date, timedelta
import numpy as np
from horizons_experiment import fetch, metrics, daynum, iso

ROOT=Path(__file__).resolve().parent
# Windows are broad search intervals, not guaranteed optimal launch periods.
WINDOWS={2028:(date(2028,7,1),date(2028,12,31)),2031:(date(2030,10,1),date(2031,4,1)),2033:(date(2032,12,1),date(2033,6,1)),2035:(date(2035,2,1),date(2035,8,1))}
C3=25.; DV=6.2; MIN_TOF=140; MAX_TOF=321
DESIGN=np.arange(15); STRESS=np.arange(15,22); EXTREME=np.arange(22,31)
# Hypothetical delay distributions, not empirical launch reliability.
DISTRIBUTIONS={'uniform_0_30':np.ones(31),'early_geometric':.85**np.arange(31),'late_weighted':np.arange(1,32,dtype=float)}

def weighted_quantile(values,weights,q):
    order=np.argsort(values); x=np.asarray(values)[order]; w=np.asarray(weights)[order]
    return float(x[min(np.searchsorted(np.cumsum(w)/sum(w),q,side='left'),len(x)-1)])

def pareto(costs):
    x=np.asarray(costs);return np.array([not np.any(np.all(x<=v,axis=1)&np.any(x<v,axis=1)) for v in x])

def write_csv(path, rows):
    if not rows:return
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def run_window(year,earth,mars):
    first,last=map(daynum,WINDOWS[year]); cache={}; failed={}
    def calc(d,t):
        key=(d,t)
        if key not in cache and key not in failed:
            try:
                if t<=0:raise ValueError('nonpositive flight time')
                m=metrics(earth[d],mars[d+t],t)
                if not np.isfinite(m).all():raise ValueError('nonfinite metric')
                cache[key]=m
            except (ValueError,RuntimeError,KeyError,FloatingPointError,OverflowError) as exc:failed[key]=str(exc)
        return cache.get(key)
    def eligible(m):return m is not None and m[0]<=C3 and m[-1]<=DV
    # Coarse search, then refine both dates in a neighborhood of the best point.
    coarse=[(d,t,calc(d,t)) for d in range(first,last+1,7) for t in range(MIN_TOF,MAX_TOF+1,10)]
    valid=[x for x in coarse if eligible(x[2])]
    if not valid:return {'year':year,'status':'no feasible coarse candidate','failed':len(failed)},[],[],[]
    d0,t0,_=min(valid,key=lambda x:x[2][-1]); fine=[(d,t,calc(d,t)) for d in range(max(first,d0-7),min(last,d0+7)+1) for t in range(max(MIN_TOF,t0-10),min(MAX_TOF,t0+10)+1)]
    d0,t0,_=min([x for x in fine if eligible(x[2])]+valid,key=lambda x:x[2][-1]);arrival=d0+t0
    profiles={}; candidates=[]; failures=[]
    # Shared arrival and candidate set for all methods. Allow up to 42 days advance.
    for d in range(max(first,d0-42,arrival-MAX_TOF),min(last,d0+42,arrival-MIN_TOF)+1):
        t=arrival-d
        arr=[calc(d+k,t-k) if t-k>=MIN_TOF else None for k in range(31)]
        if any(x is None for x in arr):
            failures.append({'year':year,'launch':iso(d),'reason':'missing or out-of-range delay profile'});continue
        a=np.asarray(arr); v=a[:,-1]; ok=(a[:,0]<=C3)&(v<=DV)
        row={'year':year,'launch':iso(d),'arrival':iso(arrival),'tof_days':t,'c3':float(a[0,0]),'nominal_dv':float(v[0]),'design_p95':weighted_quantile(v[DESIGN],np.ones(15),.95),'design_max':float(max(v[DESIGN])),'design_mean':float(np.mean(v[DESIGN])),'design_feasible':float(ok[DESIGN].mean())}
        # Linear interpolated percentile for comparison with v03, separately from empirical quantile.
        row['design_p95_linear']=float(np.quantile(v[DESIGN],.95))
        candidates.append(row);profiles[row['launch']]=(a,ok)
    feasible=[r for r in candidates if r['c3']<=C3 and r['nominal_dv']<=DV]
    if not feasible:return {'year':year,'status':'no feasible complete profiles','failed':len(failed),'profile_failures':len(failures)},candidates,[],[]
    nom=min(feasible,key=lambda r:r['nominal_dv'])
    design=[r for r in feasible if r['design_feasible']==1.0]
    # Design-only selection; never reads test scores before choosing policies.
    selected={'nominal':nom,'buffer_7':next((r for r in candidates if daynum(date.fromisoformat(nom['launch']))-daynum(date.fromisoformat(r['launch']))==7),None),'buffer_14':next((r for r in candidates if daynum(date.fromisoformat(nom['launch']))-daynum(date.fromisoformat(r['launch']))==14),None),'robust_p95':min(design,key=lambda r:r['design_p95_linear']) if design else None,'robust_minimax':min(design,key=lambda r:r['design_max']) if design else None,'robust_margin':min([r for r in design if r['design_max']<=DV-.2 and max(profiles[r['launch']][0][DESIGN,0])<=C3-2],key=lambda r:r['nominal_dv'],default=None)}
    # Pareto front using design-only objectives.
    front=[r for r,flag in zip(design,pareto([[r['nominal_dv'],r['design_p95_linear'],r['design_max']] for r in design])) if flag]
    records=[];samples=[]
    for name,r in selected.items():
        if r is None:continue
        a,ok=profiles[r['launch']];v=a[:,-1]
        rec=dict(r,policy=name,nominal_penalty_m_s=(r['nominal_dv']-nom['nominal_dv'])*1000,stress_feasible=float(ok[STRESS].mean()),extreme_feasible=float(ok[EXTREME].mean()),stress_max=float(max(v[STRESS])),extreme_max=float(max(v[EXTREME])))
        for k,w in DISTRIBUTIONS.items():
            rec[k+'_feasible']=float(np.dot(ok,w)/sum(w));rec[k+'_p95']=weighted_quantile(v,w,.95)
        for isp in (320,450):rec[f'ideal_propellant_fraction_isp{isp}']=float(1-math.exp(-r['nominal_dv']*1000/(9.80665*isp)))
        records.append(rec)
        for k in range(31):samples.append({'year':year,'policy':name,'launch':r['launch'],'arrival':r['arrival'],'delay_days':k,'c3':float(a[k,0]),'dv':float(v[k]),'feasible':int(ok[k]),'split':'design' if k<=14 else 'stress' if k<=21 else 'extreme'})
    # Test data are scored for every selected policy; absent policies explicitly reported.
    summary={'year':year,'status':'ok','coarse_candidates':len(coarse),'fine_candidates':len(fine),'unique_solved':len(cache),'failed':len(failed),'profile_failures':len(failures),'candidate_count':len(candidates),'design_feasible_count':len(design),'pareto_count':len(front),'arrival':iso(arrival),'missing_policies':[k for k,v in selected.items() if v is None],'policies':records}
    return summary,candidates,front,samples

def main():
    p=argparse.ArgumentParser();p.add_argument('--years',nargs='+',type=int,choices=list(WINDOWS),default=list(WINDOWS));p.add_argument('--offline',action='store_true');args=p.parse_args()
    if args.offline and any(y!=2028 for y in args.years):p.error('--offline supports only bundled 2028 data; use --years 2028')
    out=ROOT/'results_v04';out.mkdir(exist_ok=True);raw=ROOT/'horizons_raw';raw.mkdir(exist_ok=True)
    summaries=[];sources=[];all_candidates=[];all_front=[];all_samples=[]
    for y in dict.fromkeys(args.years):
        print(f'Window {y}: loading reference states...',flush=True)
        if y==2028:start,stop=date(2028,1,1),date(2029,11,20)
        else:start,stop=WINDOWS[y][0],WINDOWS[y][1]+timedelta(days=MAX_TOF)
        if args.offline and not all((raw/f'{code}_{start}_{stop}.json').exists() for code in (399,499)):raise FileNotFoundError('Bundled reference data missing')
        earth,e=fetch(399,start,stop,raw);mars,m=fetch(499,start,stop,raw);sources.extend([e,m]);s,c,f,a=run_window(y,earth,mars);summaries.append(s);all_candidates+=c;all_front+=f;all_samples+=a
        print(json.dumps({'year':y,'status':s['status'],'policies':[{'policy':r['policy'],'launch':r['launch'],'nominal_dv':r['nominal_dv'],'stress_feasible':r['stress_feasible'],'extreme_feasible':r['extreme_feasible']} for r in s.get('policies',[])]},indent=2),flush=True)
        report={'model':'Sun-only Lambert, zero revolution prograde, ideal patched conics, 300 km circular parking orbits','development_window':2028,'later_windows':'prospective benchmark relative to 2028, not statistical independent random samples','limits':{'c3':C3,'dv':DV,'tof_days':[MIN_TOF,MAX_TOF],'design_delay':[0,14],'stress_delay':[15,21],'extreme_delay':[22,30],'margin_dv':.2,'margin_c3':2},'environment':{'python':platform.python_version(),'numpy':np.__version__},'sources':sources,'windows':summaries}
        (out/'summary_v04.json').write_text(json.dumps(report,indent=2),encoding='utf-8');write_csv(out/'candidates.csv',all_candidates);write_csv(out/'pareto.csv',all_front);write_csv(out/'delay_samples.csv',all_samples);write_csv(out/'policy_comparison.csv',[r for s in summaries for r in s.get('policies',[])])
    if all_samples:
        import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
        for y in dict.fromkeys(r['year'] for r in all_samples):
            fig,ax=plt.subplots(figsize=(9,5));rr=[r for r in all_samples if r['year']==y]
            for policy in dict.fromkeys(r['policy'] for r in rr):
                pp=[r for r in rr if r['policy']==policy];ax.plot([r['delay_days'] for r in pp],[r['dv'] for r in pp],label=policy)
            ax.axhline(DV,ls='--',label='Illustrative budget');ax.axvspan(14.5,21.5,alpha=.12,label='Stress test');ax.set(xlabel='Launch delay (days; arrival fixed)',ylabel='Ideal total delta-v (km/s)',title=f'KEPLER-X v0.4 | {y}');ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/f'robustness_{y}.png',dpi=160);plt.close(fig)
    print('Saved results_v04/summary_v04.json',flush=True)

if __name__=='__main__':main()
