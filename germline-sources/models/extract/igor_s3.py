# Digitize Marcou et al. 2018 Supplementary Figure 3 (TRB allele usage; blue = DNA, red = mRNA
# individuals) from the vector SVG of supplement page 15. Prints JSON {allele: {dna:[..], rna:[..]}}.
import json,subprocess,sys,collections
rows=[l.split('|') for l in subprocess.run([sys.executable,'-I',sys.argv[2],sys.argv[1]],capture_output=True,text=True).stdout.splitlines()]
labels=[l.split('] ',1)[1].split() for l in open(sys.argv[3]).read().strip().splitlines() if l.strip()]
X0,X1=44.094,450.750
# (bottom y, y per unit, number of slots, label row, prefix)
panels=[(61.570,(61.570-55.70)/0.01,30,0,'TRBV'),(148.727,(148.727-140.90)/0.02,30,1,'TRBV'),(235.887,(235.887-226.49)/0.02,30,2,'TRBV'),
        (367.934,(367.934-358.54)/0.05,16,3,'TRBJ'),(499.980,(499.980-494.11)/0.1,3,4,'TRBD')]
COL={'rgb(0%, 0%, 100%)':'dna','rgb(100%, 0%, 0%)':'rgb'}
out={}
for bot,per,n,li,pre in panels:
    g=collections.defaultdict(lambda:{'dna':[],'rna':[]})
    for r in rows:
        if r[0] in ('rgb(0%, 0%, 100%)','rgb(100%, 0%, 0%)'):
            a,b,c,d=map(float,r[2:6])
            if abs(d-bot)<0.05 and X0<=a and b<=X1:
                i=int(((a+b)/2-X0)/((X1-X0)/n))
                g[i]['dna' if r[0]=='rgb(0%, 0%, 100%)' else 'rna'].append(round((bot-c)/per,5))
    for i,name in enumerate(labels[li]): out[pre+name]=g[i]
print(json.dumps(out))
