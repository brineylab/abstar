# Median lines of ggplot box plots in a raster figure (Mikocziova et al. 2021 iScience Figs 3-4).
# usage: boxmedians.py image ytop_px ytop_val ybot_px ybot_val half_width x1 x2 ...
import sys
from PIL import Image
im=Image.open(sys.argv[1]).convert('RGB'); px=im.load()
yt,vt,yb,vb,hw=float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5]),int(sys.argv[6])
dark=lambda p: p[0]<80 and p[1]<80
for x in map(int,sys.argv[7:]):
    rows=[y for y in range(int(yt)-5,int(yb)+6) if sum(dark(px[c,y]) for c in range(x-hw,x+hw+1))>=0.85*(2*hw+1)]
    runs=[];
    for y in rows:
        if runs and y==runs[-1][1]+1: runs[-1][1]=y
        else: runs.append([y,y])
    # median = thickest horizontal run (ggplot draws it at 2x line width)
    best=max(runs,key=lambda r:(r[1]-r[0],-abs((r[0]+r[1])/2-(runs[0][0]+runs[-1][1])/2))) if runs else None
    val=None if best is None else vb+((best[0]+best[1])/2-yb)*(vt-vb)/(yt-yb)
    print(x,[tuple(r) for r in runs],None if val is None else round(val,5))
