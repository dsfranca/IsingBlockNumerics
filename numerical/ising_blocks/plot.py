#!/usr/bin/env python3
"""Build the comparison table and a publication-ready PGFPlots figure."""
import csv
import subprocess
import tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent
xs=[.30,.35,.40,.42]
grid={}
for x in xs:
    for r in csv.DictReader((ROOT/f'results_{x:.2f}.csv').open()):
        grid[(x,int(r['ell']))]=r
with (ROOT/'margin_table.csv').open('w',newline='') as f:
    writer=csv.writer(f)
    writer.writerow(['ell']+[f'x={x:.2f}' for x in xs])
    for ell in range(1,7):writer.writerow([ell]+[grid[(x,ell)]['kappa_estimate'] for x in xs])
colors=['444444','0C7C86','D28C17','7851A9','CB4B56','246EB9']
lines=[r'\documentclass[tikz,border=3pt]{standalone}',r'\usepackage{pgfplots}',r'\pgfplotsset{compat=1.18}']
for ell,color in enumerate(colors,1):lines.append(r'\definecolor{line%d}{HTML}{%s}'%(ell,color))
lines += [r'\begin{document}',r'\begin{tikzpicture}',r'\begin{axis}[width=14.1cm,height=8.2cm,',
          r'xlabel={Ising coupling $x=\beta J$},ylabel={Block margin $\kappa(x,\ell)$},',
          r'xmin=.295,xmax=.425,ymin=-.9,ymax=.48,xtick={.30,.35,.40,.42},',
          r'xticklabel style={/pgf/number format/fixed,/pgf/number format/precision=2},',
          r'axis lines=left,ymajorgrids=true,grid style={black!12},',
          r'legend columns=3,legend style={at={(.025,.025)},anchor=south west,draw=none,fill=none,font=\small},',
          r'ticklabel style={font=\small}]',r'\addplot[black!65,dashed,forget plot] coordinates {(.295,0) (.425,0)};']
for ell in range(1,7):
    coords=' '.join('(%s,%s)'%(x,grid[(x,ell)]['kappa_estimate']) for x in xs)
    lines.append(r'\addplot[color=line%d,mark=*,mark size=1.7pt,line width=.9pt] coordinates {%s};'%(ell,coords))
    lines.append(r'\addlegendentry{$\ell=%d$}'%ell)
lines += [r'\end{axis}',r'\end{tikzpicture}',r'\end{document}']
tex=ROOT/'ising_block_margins.tex';tex.write_text('\n'.join(lines)+'\n')
with tempfile.TemporaryDirectory(prefix='ising-plot-') as td:
    out=Path(td)
    proc=subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error',f'-output-directory={out}',str(tex)],capture_output=True,text=True)
    if proc.returncode:
        print(proc.stdout)
        raise RuntimeError('PGFPlots build failed')
    (ROOT/'ising_block_margins.pdf').write_bytes((out/'ising_block_margins.pdf').read_bytes())
subprocess.run(['pdftoppm','-f','1','-l','1','-scale-to','1500','-png','-singlefile',str(ROOT/'ising_block_margins.pdf'),str(ROOT/'ising_block_margins')],check=True)
print((ROOT/'margin_table.csv').read_text())
