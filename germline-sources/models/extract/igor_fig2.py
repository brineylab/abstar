# Digitize Marcou et al. 2018 Fig. 2a-c (vector, article PDF page 3) for TRB and TRA. Prints JSON.
import json,subprocess,sys
rows=[l.split('|') for l in subprocess.run([sys.executable,'-I',sys.argv[2],sys.argv[1],'--pts'],capture_output=True,text=True).stdout.splitlines()]
C={'blue':'rgb(6.666565%, 0%, 100%)','greyblue':'rgb(26.272583%, 41.960144%, 67.0578%)','cyan':'rgb(0%, 81.959534%, 100%)',
   'teal':'rgb(8.235168%, 79.214478%, 85.488892%)','orange':'rgb(100%, 38.430786%, 0%)'}
def curve(col,x0,x1,v0,v1,ybot,ytop,yval):
    out={}
    for r in rows:
        if r[1]==C[col] and len(r)>7 and x0-1<=float(r[2])<=x1 and float(r[3])<=x1+1:
            pts=sorted((v0+(x-x0)*(v1-v0)/(x1-x0),(ybot-y)/(ybot-ytop)*yval) for x,y in (map(float,p.split(',')) for p in r[7].split()))
            # matplotlib drops collinear vertices: interpolate the plotted polyline at integers
            for ki in range(round(pts[0][0]),round(pts[-1][0])+1):
                for (xa,ya),(xb,yb) in zip(pts,pts[1:]):
                    if xa-0.05<=ki<=xb+0.05 and xb>xa:
                        out[ki]=round(max(0.0,ya+(yb-ya)*(ki-xa)/(xb-xa)),5); break
            return out
    raise SystemExit('missing '+col)
A=(107.36,221.63,0,60,203.82,63.87,0.16)
B=(252.23,366.50,-5,20,203.82,63.87,0.25)
Cc=(397.10,511.37,-5,20,203.82,63.87,0.14)
res={'trb_vd_dna':curve('blue',*A),'trb_dj_dna':curve('greyblue',*A),'trb_vd_rna':curve('cyan',*A),'trb_dj_rna':curve('teal',*A),
     'tra_vj_rna':curve('orange',*A),
     'trb_v3_dna':curve('blue',*B),'trb_v3_rna':curve('cyan',*B),'tra_v3_rna':curve('orange',*B),
     'trb_j5_dna':curve('blue',*Cc),'trb_j5_rna':curve('cyan',*Cc),'tra_j5_rna':curve('orange',*Cc)}
for k,d in res.items(): print(k,min(d),max(d),'sum=%.4f'%sum(d.values()),file=sys.stderr)
print(json.dumps(res,indent=1))
