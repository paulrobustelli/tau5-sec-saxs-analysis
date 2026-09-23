from pathlib import Path
from PIL import Image,ImageDraw,ImageFont
from matplotlib import font_manager
p=Path(__file__).resolve().parent
fontpath=font_manager.findfont('DejaVu Sans')
def f(n):return ImageFont.truetype(fontpath,n)
pages=[]
for i in [1,2]:
 out=Image.new('RGB',(2800,1540),'white');d=ImageDraw.Draw(out)
 for c,s in enumerate(['WT','AA']):
  out.paste(Image.open(p/f'{s}_density_view{i}.png').convert('RGB'),(c*1400,75));d.text((c*1400+700,25),s,font=f(55),fill='navy' if s=='WT' else 'orangered',anchor='mt')
 d.text((1400,1490),'PyMOL / DENSS • Blue: 90% density   Cyan: 60%   Red: 30% • Same spatial scale',font=f(31),fill='#444444',anchor='mm')
 out.save(p/f'DENSS_PyMOL_view{i}.png',dpi=(300,300));pages.append(out)
pages[0].save(p/'DENSS_PyMOL_figures.pdf',resolution=300,save_all=True,append_images=pages[1:])
