from pathlib import Path
from PIL import Image,ImageDraw,ImageFont
P=Path(__file__).resolve().parent
font='/System/Library/Fonts/Helvetica.ttc'
for view in [1,2]:
 canvas=Image.new('RGB',(2100,1550),'white');d=ImageDraw.Draw(canvas)
 for col,(style,label) in enumerate([('density_beads','Space-filling density beads'),('blue','Blue envelope'),('rainbow','Rainbow density')]):
  d.text((col*700+350,45),label,font=ImageFont.truetype(font,32),fill='black',anchor='mm')
  for row,sample in enumerate(['WT','AA']):
   im=Image.open(P/f'{sample}_{style}_view{view}.png').convert('RGB').resize((700,700),Image.Resampling.LANCZOS)
   canvas.paste(im,(col*700,90+row*700));d.text((col*700+30,110+row*700),sample,font=ImageFont.truetype(font,36),fill='black')
 d.text((1050,1520),'Same spatial scale • Beads represent density, not atoms • Blue → red: low → high relative density',font=ImageFont.truetype(font,25),fill='#444444',anchor='mm')
 canvas.save(P/f'DENSS_three_styles_view{view}.png',dpi=(300,300));canvas.save(P/f'DENSS_three_styles_view{view}.pdf',resolution=300)
