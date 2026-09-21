from run_bift import *
rows=[]
for name,a,b in [('WT_092025',157,175),('WT_092025',161,175),('WT_092025',159,173),('WT_092025',159,177),('AA_052026',167,188),('AA_052026',171,188),('AA_052026',169,186),('AA_052026',169,190)]:
 rows.append(fit_bift(name,a,b,suffix=f'window_{a}_{b}',mc=80))
(PLOT/'bift_window_sensitivity.json').write_text(json.dumps(rows,indent=2))
