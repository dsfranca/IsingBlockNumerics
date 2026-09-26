#!/usr/bin/env python3
"""Independent direct enumeration of internal and exterior Ising spins."""
import argparse
import csv
import subprocess
import tempfile
from decimal import Decimal, localcontext
import math
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent

def direct(ell,x):
    ns=1<<(ell*ell)
    s=2*((np.arange(ns,dtype=np.int64)[:,None] >> np.arange(ell*ell))&1)-1
    edges=[(r*ell+c,(r+1)*ell+c) for r in range(ell-1) for c in range(ell)]
    edges += [(r*ell+c,r*ell+c+1) for r in range(ell) for c in range(ell-1)]
    e=np.zeros(ns)
    for i,j in edges:e+=s[:,i]*s[:,j]
    # Exterior bits: top, bottom, left, right, each ordered by column/row.
    sites=list(range(ell))+list(range(ell*(ell-1),ell*ell))+list(range(0,ell*ell,ell))+list(range(ell-1,ell*ell,ell))
    count=(s.sum(axis=1)+ell*ell)/2
    means=[]
    for a in range(0,1<<(4*ell),128):
        b=2*((np.arange(a,min(a+128,1<<(4*ell)))[:,None] >> np.arange(4*ell))&1)-1
        weights=np.exp(x*(e[None,:]+b@s[:,sites].T))
        means.extend((weights@count)/weights.sum(axis=1))
    means=np.asarray(means)
    maxima=[]
    for p in range(ell):
        mask=1<<(2*ell+p)
        minus=np.arange(1<<(4*ell));minus=minus[(minus&mask)==0]
        maxima.append(float(np.max(means[minus|mask]-means[minus])))
    return maxima

for text in ('0.30','0.35','0.40','0.42'):
    rows=list(csv.DictReader((ROOT/f'results_{text}.csv').open()))
    for ell in (1,2,3):
        cs=direct(ell,float(text))
        recorded=[float(r['c_estimate']) for r in rows if int(r['ell'])==ell]
        assert max(abs(a-b) for a,b in zip(cs,recorded)) < 2e-13,(text,ell,cs,recorded)
        if ell==1:assert abs(cs[0]-.5*math.tanh(2*float(text)))<2e-15
        print(f'x={text} ell={ell}: exhaustive independent internal-spin enumeration agrees')
for ell in (1,2,3):
    assert max(abs(c) for c in direct(ell,0))<1e-14
    cs=direct(ell,1e-6)
    assert max(cs)<1.00001e-6
print('Zero-coupling and small-coupling checks passed.')

parser=argparse.ArgumentParser()
parser.add_argument('--search',help='Optional compiled search executable for transfer zero/small-coupling checks')
args=parser.parse_args()
if args.search:
    with tempfile.TemporaryDirectory(prefix='ising-check-') as td, localcontext() as ctx:
        ctx.prec=100
        for x in ('0','0.000001'):
            path=Path(td)/'weights.txt';out=Path(td)/'result.csv'
            path.write_text(x+' 9\n'+'\n'.join(format(float((Decimal(x)*k).exp()),'.17g') for k in range(-9,10))+'\n')
            subprocess.run([args.search,str(path),'1','3',str(out)],check=True,capture_output=True)
            for row in csv.DictReader(out.open()):
                c=float(row['c_estimate'])
                assert c==0 if x=='0' else 0<c<1.00001e-6
        print('Compiled transfer implementation also passes zero/small-coupling checks.')
