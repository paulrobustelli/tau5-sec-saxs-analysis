from run_efa import *
rows=[]
for endearly in [170,175,180,185,195,205]:
 for startlate in [130,140,145,150]:
  r=rotation('WT_092025',120,205,[[120,endearly],[startlate,205]])
  rows.append(r)
print('BEST',sorted(rows,key=lambda r:r.get('residual_ms',999))[:5])
(OUT/'rotation_extended_windows.json').write_text(json.dumps(rows,indent=2))
# Examine Guinier range dependence instead of forcing the shoulder to an IDP size.
for k in [1,2]:
 d=np.loadtxt(OUT/f'WT_092025_120-205_rotation_component{k}.dat')
 print('Rg scans component',k,[(qm,guinier(*d.T,qmin=qm)) for qm in [.006,.008,.010,.012,.015,.018,.02,.025]])
