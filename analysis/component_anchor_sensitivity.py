from component_recovery import *
from independent_baseline import window_stats

def anchors_audit():
 checks=[];results=[]
 for name,spec in SPECS.items():
  q,Y,E=load(name)
  for side,(a,b) in zip(['pre','post'],spec['anchors']):checks.append(dict(sample=name,side=side,**window_stats(q,Y,E,a,b)))
  for da,db in itertools.product([-5,0,5],repeat=2):
   anchors=[[a+dt,b+dt] for (a,b),dt in zip(spec['anchors'],[da,db])]
   tests=[window_stats(q,Y,E,a,b) for a,b in anchors]
   r=dict(sample=name,anchors=anchors,anchor_stationarity_pass=all(t['passes'] for t in tests),tests=tests)
   if r['anchor_stationarity_pass']:
    info,_=evaluate(name,'linear',*spec['fid'],anchors=anchors);r['recovery']=info
   results.append(r)
 (OUT/'baseline_anchor_checks.json').write_text(json.dumps(checks,indent=2));(OUT/'anchor_sensitivity.json').write_text(json.dumps(results,indent=2))
 return checks,results
if __name__=='__main__':anchors_audit()
