from pymol import cmd
import json
from pathlib import Path
OUT=Path(__file__).resolve().parent
settings=json.loads((OUT/'render_settings.json').read_text())
for s in ['WT','AA']:
 ramp=[]
 for f,c,a in [(0.015,'blue',0),(.04,'blue',.012),(.12,'cyan',.015),(.24,'green',.025),(.38,'yellow',.05),(.52,'orange',.08),(.68,'red',.12),(1.1,'red',.15)]:ramp.extend([settings[s]['peak']*f,c,a])
 cmd.volume_ramp_new(s+'_rich_ramp',ramp);cmd.volume_color(s+'_rainbow',s+'_rich_ramp');settings[s]['ramp']=ramp
for s in ['WT','AA']:
 for style in ['density_beads','blue','rainbow']:
  for v in [1,2]:
   cmd.do(f'show_denss {s}, {style}, {v}')
   cmd.png(str(OUT/f'{s}_{style}_view{v}.png'),width=1600,height=1600,dpi=300,ray=0 if style=='rainbow' else 1)
   cmd.scene(f'{s}_{style}_{v}','store')
(OUT/'render_settings.json').write_text(json.dumps(settings,indent=2))
cmd.save(str(OUT/'DENSS_three_styles.pze'))
print('FINAL RENDERS COMPLETE')
