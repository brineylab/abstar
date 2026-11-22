# Dump styled paths from a pdftocairo SVG in page coordinates (y down).
# usage: svgpaths.py file.svg [--pts]
import re,sys
t=open(sys.argv[1]).read()
def pts(d):
    nums=[float(x) for x in re.findall(r'-?\d+\.?\d*(?:e-?\d+)?',d)]
    return list(zip(nums[0::2],nums[1::2]))
for m in re.finditer(r'<path ([^>]*)/?>',t):
    a=m.group(1)
    if 'clip-rule' in a or ('fill="' not in a and 'stroke="' not in a): continue
    d=re.search(r'\bd="([^"]+)"',a).group(1)
    tr=re.search(r'transform="matrix\(([^)]+)\)"',a)
    P=pts(d)
    if tr:
        A,B,C,D,E,F=[float(x) for x in tr.group(1).split(',')]
        P=[(A*x+C*y+E,B*x+D*y+F) for x,y in P]
    fill=re.search(r'\bfill="([^"]+)"',a); fill=fill.group(1) if fill else None
    st=re.search(r'\bstroke="([^"]+)"',a); st=st.group(1) if st else None
    xs=[p[0] for p in P]; ys=[p[1] for p in P]
    s='%s|%s|%.3f|%.3f|%.3f|%.3f|%d'%(fill,st,min(xs),max(xs),min(ys),max(ys),len(P))
    if '--pts' in sys.argv and len(P)>2: s+='|'+' '.join('%.3f,%.3f'%p for p in P)
    print(s)
