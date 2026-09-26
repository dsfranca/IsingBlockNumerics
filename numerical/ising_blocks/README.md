# Finite Ising block-influence certificates

This directory contains the full exhaustive boundary search used in the 2D
transverse-field Ising example. No strong-spatial-mixing assumption enters this
finite computation. The reference is the ferromagnetic square-lattice Ising
model at zero longitudinal field. A square block has side length `ell`; its
`4 ell` exterior neighbors are all distinct (on a torus this holds for
`L >= ell + 2`).

The manuscript keeps only the mathematical certificate in Appendix E.
Detailed derivations, the full comparison table and plot, roundoff analysis,
and explicit shell estimates are retained in
[`certification_details.tex`](certification_details.tex), a LaTeX source
fragment using the manuscript notation and macros. A typeset standalone companion is available as [certification_note.pdf](../../output/pdf/certification_note.pdf), with [source](../../certification_note.tex). The scripts and data below
provide the reproducible numerical supplement.

## Reproduce

Run all commands below from the repository root. Use Python 3.10 or later and a C++17 GCC or Clang compiler. The independent verification/plot scripts
use NumPy and a LaTeX installation with PGFPlots, respectively.

```sh
python3 numerical/ising_blocks/weights.py
c++ -O3 -std=c++17 -ffp-contract=off numerical/ising_blocks/search.cpp -o /tmp/ising-block-search
/tmp/ising-block-search numerical/ising_blocks/weights_0.30.txt 1 6 numerical/ising_blocks/results_0.30.csv
/tmp/ising-block-search numerical/ising_blocks/weights_0.35.txt 1 6 numerical/ising_blocks/results_0.35.csv
/tmp/ising-block-search numerical/ising_blocks/weights_0.40.txt 1 6 numerical/ising_blocks/results_0.40.csv
/tmp/ising-block-search numerical/ising_blocks/weights_0.42.txt 1 6 numerical/ising_blocks/results_0.42.csv
/tmp/ising-block-search numerical/ising_blocks/weights_0.36.txt 7 8 numerical/ising_blocks/extra_0.36.csv
python3 numerical/ising_blocks/verify.py --search /tmp/ising-block-search
python3 numerical/ising_blocks/plot.py
```

The additional `extra_*.csv` files extend selected searches to sides 7 and 8. To reproduce the remaining supplementary tables (these searches are substantially more expensive):

```sh
/tmp/ising-block-search numerical/ising_blocks/weights_0.35.txt 7 8 numerical/ising_blocks/extra_0.35.csv
/tmp/ising-block-search numerical/ising_blocks/weights_0.37.txt 7 8 numerical/ising_blocks/extra_0.37.csv
/tmp/ising-block-search numerical/ising_blocks/weights_0.40.txt 7 8 numerical/ising_blocks/extra_0.40.csv
```

The comparison table and figure in this numerical supplement use sides 1–6.
Do not enable `-ffast-math` or reassociation. IEEE binary64, round-to-nearest
arithmetic is checked at startup; fused contraction is disabled in the build
command. Results were produced on Apple silicon using Apple clang.

## Quantity and exhaustive search

For exterior site `j` at position `p` on the left side, FKG monotonicity and
Hamming-Wasserstein duality give

```
c_p = max_{other exterior spins} (E_+[number of plus spins]
                                - E_-[number of plus spins]).
```

A column `s` has `2^ell` states. The vertical and top/bottom factor is
`V_c(s) = exp(x [sum_r s_r s_(r+1) + top_c s_1 + bottom_c s_ell])`.
The intercolumn factor is `K(s,t) = exp(x sum_r s_r t_r)`, the tensor product
of `ell` copies of `[[exp(x),exp(-x)],[exp(-x),exp(x)]]`.
For each assignment of the top/bottom exterior spins, the code propagates
matrices indexed by the entire left boundary assignment and the current
column state. It propagates both the partition function `Z` and its derivative
`U` with respect to a uniform field coupled to the number of plus spins:

```
Z <- (Z K) diag(V_c)
U <- [(U K) + (Z K) diag(n_plus)] diag(V_c).
```

The final multiplication by `K` indexes all right boundary assignments.
Thus `U/Z` is the exact mean plus-spin count before numerical roundoff.
All values entering this recurrence are nonnegative. No finite difference is
used. The tensor-product transform computes each matrix multiplication by
`K` without a dense matrix product.

Global spin flip preserves the positive response difference. We therefore
fix the top-left exterior spin to plus and enumerate all `2^(2ell-1)`
remaining top/bottom patterns; all `2^(2ell)` left/right patterns are handled
in each matrix. For every `p`, all `2^(4ell-2)` representatives of the boundary
assignments other than `j` are checked. No conjecture about the maximizing
boundary, reflection reduction, or unproved monotonicity in `x` is used.
The all-translate block family with weight `ell^-2` has margin

```
kappa(x,ell) = 1 - (4/ell^2) sum_(p=1)^ell c_p.
```

## Roundoff certificate

Every CSV contains a numerical maximizer and certified enclosing bounds.
Here are the details of the arithmetic certificate; these bounds apply to
all enumerated boundary assignments, including nonmaximizers.

1. All couplings `x` are exact terminating rationals. `weights.py` computes
   `exp(k*x)`, `|k| <= 9`, at decimal precision 100. Python's `Decimal.exp`
   is correctly rounded. On the range used here, its absolute error is less
   than `10^-98`. The script checks, using the exact decimal representation
   of each binary64 result, that its relative error (including that
   `10^-98` enclosure) is less than `2^-52 = 2 u`, where `u = 2^-53`.
   The 17-digit decimal serialization round-trips these binary64 inputs.
   No runtime exponential is evaluated by the C++ search.
2. For `ell <= 8`, each operation path through the nonnegative recurrence
   uses at most `4 ell^2 + 5 ell + 5 <= 301` unit-roundoff factors. Each
   tensor-transform coordinate contributes at most four: input-weight
   error contributes two, multiplication one, and addition one. Each
   column update for the derivative contributes at most five, including
   multiplication by the integer plus-spin count, its addition to the
   prior derivative, and the column weight. Initial inputs are covered
   by the remaining five factors. Counting 4096 factors instead of 301
   is therefore conservative. Standard induction for nonnegative
   sums/products bounds the relative error of both `Z` and nonzero `U`
   by `gamma = 4096 u / (1 - 4096 u) < 4.55 * 10^-13`. A zero `U`
   remains exactly zero.
3. All positive intermediates are normal and finite. A partial path has
   at most `2 ell (ell-1) + 4 ell <= 144` interactions, at most `ell^2`
   plus-spin factors, and at most `2^(ell^2)` configurations. For the
   tested `0 <= x <= .42`, all positive intermediates lie between
   `exp(-61)` and `64 * 2^64 * exp(61)`, well inside binary64's normal
   range. Factors and all nonzero integer counts lie in this range too.
4. Dividing `U` by `Z` gives a mean in `[0,ell^2]`. Its absolute error
   is at most
   `ell^2 * [((1+gamma)/(1-gamma))*(1+u)-1]`.
   Subtracting two such means adds at most `2 ell^2 u` plus a negligible
   product term. For `ell <= 8`, the total error is less than
   `1.3 * 10^-10`. We use the looser common bound `10^-9` for each
   boundary response. Taking a maximum preserves this absolute-error
   bound. Endpoints are rounded outward with `nextafter`.
5. The displayed margin estimate is enclosed with error
   `4*10^-9/ell + 10^-13`, where the second term covers the final
   summation and affine evaluation. All signs asserted in the
   manuscript have margins many orders of magnitude larger than
   this enclosure.

The computational certificate assumes the documented correctly rounded
Decimal exponential and IEEE binary64 arithmetic; it does not treat a
floating-point optimum as an exact value.

## Data conventions

CSV column `p` is one-based from the top to the bottom along the left side.
For all boundary encodings, bit 0 represents minus and bit 1 plus.
`top_bottom_code` stores the top boundary in the lowest `ell` bits and
the bottom boundary in the next `ell` bits; each is ordered left to right.
`left_code_minus` and `right_code` order spins top to bottom. In the former,
bit `p-1` is zero; setting it to one gives the compared plus configuration.
`patterns` counts top/bottom assignments after global-spin-flip reduction;
the full number of exterior configurations represented is
`patterns * 2^(2ell)`.

The recorded boundary maximizes the computed response; ties/near-ties can
select symmetry-related representatives differently on different platforms.
The interval enclosure of the exact maximum is independent of that choice.
No uniqueness or structured form of the exact maximizing boundary is claimed.

`margin_table.csv` contains the requested grid. `ising_block_margins.pdf`
and its PNG version plot the four computed values for each side length;
line segments are visual guides, not certified interpolation bounds.

## Independent validation and conclusions

`verify.py` enumerates all `2^(ell^2)` internal spin configurations and all
`2^(4ell)` exterior configurations directly for `ell = 1,2,3` at all four
tabulated couplings. This independent calculation agrees with the transfer
search to less than `2*10^-13`. It also checks the exact single-site identity
`c = tanh(2x)/2`, zero influence at `x=0`, and the small-coupling limit.
Reflection-related coefficients agree within the certified error.

Among sides 1–6, the smallest positive margins occur at side 2 for `x=.30`
and side 6 for `x=.35`. At `x=.40` and `.42`, none of these six block sizes
passes. No pair in this grid has margin at least one half. These statements
only concern the finite grid searched; they are not impossibility claims for
larger blocks or other dynamics. The supplementary search at `x=.36`, side 7,
has positive margin approximately `.01255115638`.

## Explicit perturbation windows

`perturbation_window.py` evaluates the constants derived in the
"Explicit shell constants" section of
[`certification_details.tex`](certification_details.tex),
using upward-rounded bounds for positive constants and downward-rounded
bounds for the sufficient perturbation threshold. It bounds the entire
infinite shell tail by a geometric majorant. Run:

```sh
python3 numerical/ising_blocks/perturbation_window.py
```

The default cases use the certified lower margins `0.0432` at
`(x,ell)=(0.35,6)` and `0.0125` at `(0.36,7)`. With `u0=1`, the script
certifies the sufficient windows `beta*abs(g) <= 1e-68` and `1e-76`,
respectively. The script prints JSON to standard output. To refresh the stored `perturbation_windows.json`, redirect that output to `numerical/ising_blocks/perturbation_windows.json`. The full outward-rounded constants are recorded in
`perturbation_windows.json`. These conservative windows include both
the QBP and thermally dressed pinching contributions.
