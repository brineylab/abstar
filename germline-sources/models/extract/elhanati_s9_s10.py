# Elhanati et al. 2015 Fig. S9 (mononucleotide) and S10 (dinucleotide) "VD sense" bars (dark blue),
# individuals A and B, supplement pages 10-11. Each panel's four bars are normalised to sum to 1
# (they are probabilities over A/C/G/T), so no y-axis calibration is needed. Prints JSON.
import json,subprocess,sys
def bars(svg):
    rows=[l.split('|') for l in subprocess.run([sys.executable,'-I',sys.argv[3],svg],capture_output=True,text=True).stdout.splitlines()]
    b=[tuple(map(float,r[2:6])) for r in rows if r[0].startswith('rgb(20.800781%') and r[3]!='363.14' and float(r[3])-float(r[2])<14]
    panels={}
    for x0,x1,y0,y1 in b: panels.setdefault(round(y1,1),[]).append((x0,y1-y0))
    out=[]
    for bot in sorted(panels):
        h=[v for _,v in sorted(panels[bot])]
        if len(h)==4: out.append([round(v/sum(h),4) for v in h])
    return out
mono=bars(sys.argv[1]); di=bars(sys.argv[2])
assert len(mono)==2 and len(di)==8,(len(mono),len(di))
res={'start_A':mono[0],'start_B':mono[1],
     'trans_A':dict(zip('ACGT',di[:4])),'trans_B':dict(zip('ACGT',di[4:]))}
print(json.dumps(res,indent=1))
