# Marcou et al. 2018 Supplementary Fig. 4: IGH P(D|J), one panel per D gene (supplement page 16).
# Each panel is 0.30 tall (ticks every 0.05); the blue line has one vertex per J (IGHJ1..IGHJ6).
import json,subprocess,sys
rows=[l.split('|') for l in subprocess.run([sys.executable,'-I',sys.argv[2],sys.argv[1],'--pts'],capture_output=True,text=True).stdout.splitlines()]
D=['1-1','1-14','1-20','1-26','1-7','2-15','2-2','2-21','2-8','3-10','3-16','3-22','3-3','3-9','4-11','4-17','4-23','4-4',
   '5-12','5-18','5-24','5-5','6-13','6-19','6-25','6-6','7-27']
boxes=sorted((float(r[2]),float(r[5])) for r in rows if r[0]=='rgb(100%, 100%, 100%)' and abs(float(r[5])-float(r[4])-47.141)<0.01)
assert len(boxes)==27,len(boxes)
lines=[r for r in rows if r[1]=='rgb(0%, 0%, 100%)' and len(r)>7 and abs(float(r[3])-float(r[2])-119.371)<0.01]
res={}
for name,(x0,bot) in zip(D,boxes):   # boxes sorted by column x then row y: column-major, as labelled
    c=[r for r in lines if abs(float(r[2])-x0)<0.01 and bot-47.2<=float(r[5])<=bot+0.01]
    assert len(c)==1,(name,len(c))
    pts=sorted(tuple(map(float,p.split(','))) for p in c[0][7].split())
    # matplotlib drops collinear vertices: interpolate at the six J positions
    def at(x): return next(ya+(yb-ya)*(x-xa)/(xb-xa) for (xa,ya),(xb,yb) in zip(pts,pts[1:]) if xa-0.01<=x<=xb+0.01)
    res['IGHD'+name]=[round(max(0.0,(bot-at(x0+k*119.371/5))/(47.141/0.30)),5) for k in range(6)]
J=['IGHJ%d'%i for i in range(1,7)]
out={j:{d:v[i] for d,v in res.items()} for i,j in enumerate(J)}
for j in J: print(j,'sum=%.4f'%sum(out[j].values()),file=sys.stderr)
print(json.dumps(out,indent=1))
