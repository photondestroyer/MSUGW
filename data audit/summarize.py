import json, sys
r=json.load(open("results/scan_all.json"))
skip=sys.argv[1:] 
def f(x): return "-" if x is None else (f"{x:.4g}" if isinstance(x,(int,float)) else str(x))
for x in sorted(r,key=lambda z:(z['file'],z['var'])):
    if any(s in x['file'] for s in skip): continue
    if 'fatal' in x: print(x['file'],x['var'],"FATAL",x['fatal']); continue
    if 'valid_frac' not in x: print(x['file'],x['var'],"(char)",x.get('char_sample')); continue
    flags=[]
    if x['errors']: flags.append(f"READERR={len(x['errors'])}")
    if x['n_inf']: flags.append(f"INF={x['n_inf']}")
    if x['sentinels']: flags.append(f"sent={x['sentinels']}")
    if x['n_outside_phys']: flags.append(f"outPhys={x['n_outside_phys']}")
    if x['n_outside_valid_range_attr']: flags.append(f"outVR={x['n_outside_valid_range_attr']}")
    if x['n_outside_valid_minmax_attr']: flags.append(f"outVminmax={x['n_outside_valid_minmax_attr']}")
    nanundecl = x['n_nan'] if (x['fill'] is not None and x['fill']==x['fill'] and x['n_nan']) else 0
    if nanundecl: flags.append(f"NaN-with-nonNaN-fill={nanundecl}")
    p=x.get('pcts',{})
    print(f"{x['file'][:26]:26s} {x['var'][:30]:30s} {x['dtype']:7s} valid={x['valid_frac']:.4f} min={f(x['min'])} p1={f(p.get('1'))} med={f(p.get('50'))} p99={f(p.get('99'))} max={f(x['max'])} mean={f(x['mean'])} {' '.join(flags)}")
