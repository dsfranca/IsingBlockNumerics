#!/usr/bin/env python3
"""Certify the double weights exp(k*x) using decimal's correctly rounded exp.
All energies are rational, x in {3/10,7/20,9/25,37/100,2/5,21/50}, |k| <= 9.
The 100-digit correctly rounded Decimal exp has absolute error < 1e-98 here.
We check that converting it to binary64 incurs relative error < 2^-52,
including that 1e-98 enclosure. No runtime transcendental functions are used.
"""
from decimal import Decimal, localcontext
from pathlib import Path
ROOT = Path(__file__).resolve().parent
with localcontext() as ctx:
    ctx.prec = 100
    for text in ('0.30', '0.35', '0.36', '0.37', '0.40', '0.42'):
        vals=[]
        for k in range(-9,10):
            q=(Decimal(text)*k).exp()
            a=float(q)
            assert abs(Decimal.from_float(a)-q)+Decimal('1e-98') < q*(Decimal(2)**-52)
            vals.append(a.hex())
        # Decimal representations read correctly by C++ stream extraction.
        (ROOT / f'weights_{text}.txt').write_text(text+' 9\n'+'\n'.join(format(float.fromhex(v), '.17g') for v in vals)+'\n')
