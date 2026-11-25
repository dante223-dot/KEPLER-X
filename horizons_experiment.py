"""KEPLER-X v0.3: preliminary two-body transfer study with daily Horizons states.
Run online on the user's computer. --self-test uses synthetic parser fixtures and
Keplerian states only; never labels them as real Horizons observations.
"""
import argparse, csv, json, urllib.parse, urllib.request, hashlib
from pathlib import Path
from datetime import date,timedelta
import numpy as np
from keplerx import EPOCH, DAY, state, lambert
OUT=Path(__file__).resolve().parent
START=date(2028,1,1); STOP=date(2029,11,20)
C3_LIMIT=25.; DV_LIMIT=6.2

def daynum(d): return (d-EPOCH).days

def iso(d): return str(EPOCH+timedelta(days=int(d)))

def parse(raw,start,stop):
    if '$$SOE' not in raw or '$$EOE' not in raw: raise ValueError('Missing Horizons data delimiters: '+raw[:300])
    out={}
    for row in csv.reader(raw.split('$$SOE',1)[1].split('$$EOE',1)[0].strip().splitlines()):
        jd=float(row[0]); d=round(jd-2451544.5)
        if abs(jd-2451544.5-d)>1e-6: raise ValueError('Not midnight TDB')
        vec=np.array([float(x) for x in row[2:8]])
        if d in out or vec.shape!=(6,) or not np.isfinite(vec).all(): raise ValueError('Invalid or duplicate vector')
        out[d]=vec
    if set(out)!=set(range(daynum(start),daynum(stop)+1)): raise ValueError('Incomplete daily coverage')
    return out

def fetch(code,start,stop,cache):
    args=dict(format='json',COMMAND=f"'{code}'",CENTER="'500@10'",EPHEM_TYPE='VECTORS',START_TIME=f"'{start}'",STOP_TIME=f"'{stop}'",STEP_SIZE="'1 d'",TIME_TYPE='TDB',OUT_UNITS='KM-S',REF_PLANE='ECLIPTIC',REF_SYSTEM='ICRF',VEC_CORR='NONE',VEC_TABLE='2',CSV_FORMAT='YES',OBJ_DATA='NO')
    url='https://ssd.jpl.nasa.gov/api/horizons.api?'+urllib.parse.urlencode(args)
    path=cache/f'{code}_{start}_{stop}.json'
    if path.exists(): payload=json.loads(path.read_text(encoding='utf-8')); assert payload['url']==url
    else:
        print(f'Downloading JPL Horizons target {code}...',flush=True)
        with urllib.request.urlopen(url,timeout=90) as f: response=json.load(f)
        if response.get('error'): raise RuntimeError(response['error'])
        payload={'url':url,'response':response}
    raw=payload['response'].get('result',''); vectors=parse(raw,start,stop)
    # Archive complete response, including metadata; reject implausible state magnitudes.
    if not all(5e7<np.linalg.norm(x[:3])<4e8 and 5<np.linalg.norm(x[3:])<60 for x in vectors.values()): raise ValueError('Unexpected heliocentric state magnitude')
    if not path.exists(): path.write_bytes(json.dumps(payload,indent=2).encode('utf-8'))
    return vectors,{'target':code,'file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'url':url}

def metrics(a,b,t):
    v1,v2=lambert(a[:3],b[:3],t*DAY); vin=np.linalg.norm(v1-a[3:]); vout=np.linalg.norm(v2-b[3:])
    dep=np.sqrt(vin**2+2*398600.435/6678.137)-np.sqrt(398600.435/6678.137)
    cap=np.sqrt(vout**2+2*42828.375/3696.19)-np.sqrt(42828.375/3696.19)
    return np.array([vin**2,vout,dep,cap,dep+cap])

def approx(p,d): return np.concatenate(state(p,d))

def csv_write(path,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader();w.writerows(rows)

def evaluate(earth,mars):
    rows=[]; profiles={}; failures=[]
    for d in range(daynum(START),daynum(date(2028,12,31))+1,7):
        for t in range(140,301,10):
            a=d+t
            try:
                old=metrics(approx('earth',d),approx('mars',a),t)
                arr=np.array([metrics(earth[d+k],mars[a],t-k) for k in range(22)])
            except (ValueError,RuntimeError,FloatingPointError) as e:
                failures.append({'launch':iso(d),'tof_days':t,'error':str(e)}); continue
            ok=(arr[:,0]<=C3_LIMIT)&(arr[:,-1]<=DV_LIMIT); m=arr[0]
            row=dict(launch=iso(d),arrival=iso(a),tof_days=t,approx_c3=float(old[0]),horizons_c3=float(m[0]),approx_dv=float(old[-1]),horizons_dv=float(m[-1]),dv_change_m_s=float((m[-1]-old[-1])*1000),departure_dv=float(m[2]),capture_dv=float(m[3]),arrival_vinf=float(m[1]),delay_p95_dv=float(np.quantile(arr[:15,-1],.95)),design_feasible_fraction=float(ok[:15].mean()),stress_feasible_fraction=float(ok[15:].mean()))
            rows.append(row);profiles[(row['launch'],t)]=arr
    return rows,profiles,failures

def select(rows):
    feasible=[r for r in rows if r['horizons_c3']<=C3_LIMIT and r['horizons_dv']<=DV_LIMIT]
    robust=[r for r in feasible if r['design_feasible_fraction']>=.95]
    if not robust: raise RuntimeError('No eligible robust candidate')
    nom=min(feasible,key=lambda r:r['horizons_dv']); best=min(robust,key=lambda r:r['delay_p95_dv'])
    same=[r for r in robust if r['arrival']==nom['arrival']]
    matched=min(same,key=lambda r:r['delay_p95_dv']) if same else None
    return nom,best,matched

def fixed_arrival(earth,mars,nom,profiles):
    """Daily departure refinement with common arrival; allow up to 321 days explicitly."""
    a=daynum(date.fromisoformat(nom['arrival'])); d0=daynum(date.fromisoformat(nom['launch'])); rows=[]
    for d in range(max(daynum(START),d0-21),min(daynum(date(2028,12,31)),d0+21)+1):
        t=a-d
        if not 140<=t<=321:continue
        arr=np.array([metrics(earth[d+k],mars[a],t-k) for k in range(22)])
        ok=(arr[:,0]<=C3_LIMIT)&(arr[:,-1]<=DV_LIMIT); m=arr[0]
        r=dict(launch=iso(d),arrival=iso(a),tof_days=t,horizons_c3=float(m[0]),horizons_dv=float(m[-1]),delay_p95_dv=float(np.quantile(arr[:15,-1],.95)),design_feasible_fraction=float(ok[:15].mean()),stress_feasible_fraction=float(ok[15:].mean()))
        rows.append(r);profiles[(r['launch'],t)]=arr
    eligible=[r for r in rows if r['horizons_c3']<=C3_LIMIT and r['horizons_dv']<=DV_LIMIT and r['design_feasible_fraction']>=.95]
    return rows,min(eligible,key=lambda r:r['delay_p95_dv']) if eligible else None

def run():
    cache=OUT/'horizons_raw';cache.mkdir(exist_ok=True)
    earth,me=fetch(399,START,STOP,cache);mars,mm=fetch(499,START,STOP,cache);emb,mb=fetch(3,START,STOP,cache)
    rows,profiles,failures=evaluate(earth,mars);nom,best,_=select(rows)
    fixed,matched=fixed_arrival(earth,mars,nom,profiles);csv_write(OUT/'fixed_arrival_comparison.csv',fixed)
    csv_write(OUT/'trajectory_comparison.csv',rows)
    comparisons=[]
    for p,ref,code in [('earth_moon_barycenter',emb,3),('earth_center',earth,399),('mars',mars,499)]:
        for d,x in sorted(ref.items()):
            y=approx('mars' if code==499 else 'earth',d)
            comparisons.append(dict(target=p,date=iso(d),position_difference_km=float(np.linalg.norm(x[:3]-y[:3])),velocity_difference_m_s=float(1000*np.linalg.norm(x[3:]-y[3:])),comparison='like_for_like' if code!=399 else 'includes_EMB_vs_Earth_center_offset'))
    csv_write(OUT/'ephemeris_comparison.csv',comparisons)
    summary=dict(model='JPL Horizons geometric heliocentric ecliptic states; two-body Lambert and ideal patched-conic burns',time_scale='TDB',sources=[me,mm,mb],grid_success=len(rows),grid_failed=len(failures),failures=failures,nominal=nom,robust=best,same_arrival_robust=matched,fixed_arrival_max_tof_days=321,median_absolute_dv_change_m_s=float(np.median([abs(r['dv_change_m_s']) for r in rows])),p95_absolute_dv_change_m_s=float(np.quantile([abs(r['dv_change_m_s']) for r in rows],.95)),limits={'c3_km2_s2':C3_LIMIT,'dv_km_s':DV_LIMIT,'design_delay_days':[0,14],'stress_delay_days':[15,21]})
    (OUT/'summary_v03.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(9,5))
    for label,r in [('Nominal',nom),('Robust',best),('Same arrival',matched)]:
        if r:ax.plot(range(22),profiles[(r['launch'],r['tof_days'])][:,-1],marker='.',label=label)
    ax.axhline(DV_LIMIT,ls='--',label='Illustrative budget');ax.set(xlabel='Launch delay (days; arrival fixed)',ylabel='Ideal total delta-v (km/s)',title='KEPLER-X v0.3 | JPL Horizons inputs');ax.legend();fig.tight_layout();fig.savefig(OUT/'horizons_robustness.png',dpi=180);plt.close(fig)
    print(json.dumps({k:v for k,v in summary.items() if k not in ('sources','failures')},indent=2));print('Saved summary_v03.json',flush=True)

def self_test():
    d=daynum(START); v=approx('earth',d)
    raw='header\n$$SOE\n'+str(2451544.5+d)+', A.D. 2028-Jan-01 00:00, '+','.join(map(str,v))+',\n$$EOE\n'
    assert np.allclose(parse(raw,START,START)[d],v)
    try:parse(raw,START,START+timedelta(days=1))
    except ValueError:pass
    else:raise AssertionError('Missing day not rejected')
    a={i:approx('earth',i) for i in range(daynum(START),daynum(STOP)+1)};b={i:approx('mars',i) for i in a}
    rows,profiles,fail=evaluate(a,b);nom,best,_=select(rows)
    fixed,matched=fixed_arrival(a,b,nom,profiles);assert fixed and all(r['arrival']==nom['arrival'] for r in fixed)
    assert max(abs(r['dv_change_m_s']) for r in rows)<1e-8
    assert not fail
    print(f'PASS: parser, missing-day guard, {len(rows)} grid points, identical-input regression, policy selection. Synthetic/approximate inputs only; no live JPL validation.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--self-test',action='store_true');args=p.parse_args();self_test() if args.self_test else run()
