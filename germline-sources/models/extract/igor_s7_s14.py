# Marcou et al. 2018: Supplementary Fig. 7 (TRB "true model" D5'/D3' deletions, black lines, page 19)
# and Supplementary Fig. 14 (IGoR P(D2|J), orange line, page 26). Prints JSON.
import json,subprocess,sys
def rows(svg): return [l.split('|') for l in subprocess.run([sys.executable,'-I',sys.argv[3],svg,'--pts'],capture_output=True,text=True).stdout.splitlines()]
def interp(P,v0,v1):
    out={}
    P=sorted(P)
    for k in range(round(P[0][0]),round(P[-1][0])+1):
        for (xa,ya),(xb,yb) in zip(P,P[1:]):
            if xa-0.05<=k<=xb+0.05 and xb>xa: out[k]=round(max(0.0,ya+(yb-ya)*(k-xa)/(xb-xa)),5); break
    return out
r7=rows(sys.argv[1]); res={}
for key,(xlo,ybot,per) in {'trb_d5':(44.09,193.44,9.327/0.02),'trb_d3':(266.96,193.44,14.92/0.05)}.items():
    c=[r for r in r7 if r[1]=='rgb(0%, 0%, 0%)' and r[0]=='none' and len(r)>7 and abs(float(r[2])-(xlo+7.402))<1 and 190<float(r[5])<193.5]
    assert len(c)==1,(key,len(c))
    P=[((float(x)-xlo)/7.402-5,(ybot-float(y))/per) for x,y in (p.split(',') for p in c[0][7].split())]
    res[key]=interp(P,-5,20)
r14=rows(sys.argv[2])
c=[r for r in r14 if r[1]=='rgb(92.155457%, 39.215088%, 12.940979%)' and len(r)>7][0]
J=['TRBJ1-1','TRBJ1-2','TRBJ1-3','TRBJ1-4','TRBJ1-5','TRBJ1-6','TRBJ2-1','TRBJ2-2','TRBJ2-3','TRBJ2-4','TRBJ2-5','TRBJ2-6','TRBJ2-7']
res['trb_p_d2_given_j']={j:round((150.53-float(p.split(',')[1]))/(26.54/0.2),4) for j,p in zip(J,c[7].split())}
for k,d in res.items(): print(k,'sum=%.4f'%sum(d.values()),file=sys.stderr)
print(json.dumps(res,indent=1))
