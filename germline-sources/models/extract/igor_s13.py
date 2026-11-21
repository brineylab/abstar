# Digitize Marcou et al. 2018 Supplementary Figure 13 (IGoR series, orange) from the
# vector SVG of supplement page 25. Prints JSON.
import json,re,subprocess,sys
svg=sys.argv[1]
rows=[l.split('|') for l in subprocess.run([sys.executable,'-I',sys.argv[2],svg,'--pts'],capture_output=True,text=True).stdout.splitlines()]
OR='rgb(92.155457%, 39.215088%, 12.940979%)'
def bars(x0,x1,y0,y1):
    out=[]
    for r in rows:
        if r[0]==OR:
            a,b,c,d=map(float,r[2:6])
            if a>=x0 and b<=x1 and c>=y0-0.01 and d<=y1+0.01 and abs(d-y1)<0.01: out.append(((a+b)/2,c))
    return sorted(out)
def line(x0,x1,y0,y1):
    best=None
    for r in rows:
        if r[1]==OR and r[0]=='none' and len(r)>7:
            P=[tuple(map(float,p.split(','))) for p in r[7].split()]
            P=[p for p in P if x0-0.5<=p[0]<=x1+0.5 and y0-0.5<=p[1]<=y1+0.5]
            if len(P)>5 and (best is None or len(P)>len(best)): best=P
    return best
VG1=['1-18','1-2','1-24','1-3','1-38-4','1-45','1-46','1-58','1-69','1-8','1-NL1','2-26','2-5','2-70','3-11','3-13','3-15','3-16','3-19','3-20','3-21','3-22','3-23','3-30','3-30-3','3-33','3-35','3-38-3','3-43','3-43D']
VG2=['3-47','3-48','3-49','3-52','3-53','3-64','3-66','3-69-1','3-7','3-71','3-72','3-73','3-74','3-9','3-NL1','4-28','4-30-2','4-34','4-38-2','4-39','4-4','4-55','4-59','4-61','5-10-1','5-51','6-1','7-34-1','7-4-1','7-81']
DG=['1-1','1-14','1-20','1-26','1-7','2-15','2-2','2-21','2-8','3-10','3-16','3-22','3-3','3-9','4-17','4-23','4-4','5-12','5-24','5-5','6-13','6-19','6-25','6-6','7-27']
res={}
def usage(names,x0,x1,y0,y1,per):
    b=bars(x0,x1,y0,y1); assert len(b)==len(names),(len(b),len(names))
    return {n:round((y1-top)/per,5) for n,(x,top) in zip(names,b)}
v=usage(VG1,44.09,230.70,119.008,188.852,12.6995/0.02); v.update(usage(VG2,265.43,452.03,119.008,188.852,12.6995/0.02))
res['V']={'IGHV'+k:x for k,x in v.items()}
res['D']={'IGHD'+k:x for k,x in usage(DG,44.09,230.70,223.410,293.254,17.461/0.05).items()}
res['J']={'IGHJ%d'%(i+1):x for i,x in enumerate(usage(list(range(6)),265.43,452.03,223.410,293.254,13.969/0.1).values())}
def curve(x0,x1,y0,y1,vmin,vmax,ymaxval):
    P=line(x0,x1,y0,y1); per=(x1-x0)/(vmax-vmin); out={}
    pts=sorted((vmin+(x-x0)/per,(y1-y)/(y1-y0)*ymaxval) for x,y in P)
    for ki in range(vmin,vmax+1):
        # matplotlib drops collinear vertices, so linearly interpolate the plotted polyline
        for (xa,ya),(xb,yb) in zip(pts,pts[1:]):
            if xa-0.05<=ki<=xb+0.05 and xb>xa:
                out[ki]=round(max(0.0,ya+(yb-ya)*(ki-xa)/(xb-xa)),5); break
    return out
res['ins_vd']=curve(44.094,230.699,14.606,84.449,0,40,0.20)
res['ins_dj']=curve(265.426,452.031,14.606,84.449,0,40,0.20)
res['del_v3']=curve(44.094,230.699,327.816,397.660,-5,15,0.5)
res['del_j5']=curve(265.426,452.031,327.816,397.660,-5,15,0.5)
res['del_d5']=curve(44.094,230.699,432.219,502.063,-5,15,0.5)
res['del_d3']=curve(265.426,452.031,432.219,502.063,-5,15,0.5)
for k,d in res.items(): print(k,len(d),'sum=%.4f'%sum(d.values()),file=sys.stderr)
print(json.dumps(res,indent=1))
