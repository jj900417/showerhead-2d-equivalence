# showerhead-equiv

3D-to-2D equivalence of multi-hole showerhead and blocker-plate rings, as introduced in

> J. H. Yoon and J.-H. Lee, "A 3D-to-2D modelling framework of multi-hole showerheads for
> process screening of remote plasma reactors," *Nature Chemical Engineering* (submitted, 2026).
> https://doi.org/`<ARTICLE_DOI>`

For each ring of *N* identical cylindrical holes (radius *r*, length *H*, ring radius *L*) it
returns the axisymmetric (2D) annular slot that reproduces the ring's flow conductance, and the
area-corrected surface reaction probability for the slot walls:

1. **Flow-conductance equivalence (article Eq. 3).** The equivalent annular half-width *R′* solves

   *N r*⁴ = *R*<sub>o</sub>⁴ [ (1 − *k*⁴) − (1 − *k*²)² / ln(1/*k*) ],
   *R*<sub>o</sub> = *L* + *R′*, *R*<sub>i</sub> = *L* − *R′*, *k* = *R*<sub>i</sub> / *R*<sub>o</sub>,

   i.e. the summed Poiseuille conductance of the holes equals that of the annulus.
2. **Wall-loss equivalence (article Eq. 4).** *γ*<sub>eff</sub> = *γ* · *A*<sub>wall,3D</sub> / *A*<sub>wall,2D</sub>
   = *γ N r* / (2*L*), applied to the slot walls only; all other surfaces keep the original *γ*.

Assumptions: fully developed laminar flow in every hole (*H* cancels), identical equally spaced
holes, the same inlet and outlet conditions around the ring, negligible entrance and exit losses,
and a reaction-limited wall loss (surface Damköhler number below 1; Supplementary Note 1 of the
article).

A browser version runs at https://junghwanyoon.dev/tools/showerhead-calculator/v1.0 and uses the
same algorithm.

## Install

```bash
pip install .            # no dependencies beyond the standard library
pip install ".[test]"    # + pytest
```

## Use

```python
from showerhead_equiv import Ring, equivalent_2d

eq = equivalent_2d(Ring(r=0.5, N=18, L=24.0, H=5.0))   # lengths in mm (synthetic example)
eq.r_prime          # equivalent annular half-width R' [mm]
eq.area_ratio       # A_wall,3D / A_wall,2D = N r / (2 L)
eq.gamma_eff(0.015823)
```

Batch, from CSV (columns `label, r_mm, N, L_mm, H_mm [, ratio_override]` and `name, gamma`):

```bash
showerhead-equiv examples/geometry_example.csv \
    --reactions examples/reactions_example.csv -o result.csv
```

| Input | Meaning | Unit |
|---|---|---|
| `r` | hole radius | mm |
| `N` | number of holes in the ring | – |
| `L` | radial position of the hole centres | mm |
| `H` | hole length (plate thickness); not needed by the equations, kept for the record | mm |
| `ratio_override` | optional user-supplied *A*<sub>wall,3D</sub> / *A*<sub>wall,2D</sub> | – |

## Example (`examples/`)

`geometry_example.csv` is a **synthetic** six-ring plate; `reactions_example.csv` holds the
article's two surface reactions (F + SiO₂ etching, pre-exponential *γ*₀ = 1.5823 × 10⁻²;
F recombination, *γ* = 7 × 10⁻⁵). `expected_output.csv` is the reference output, checked by
`pytest`. For the Arrhenius reaction R1 the pre-exponential coefficient *γ*₀ is entered; the area
correction is linear, so it applies to *γ*₀ directly.

The hole layouts of the reactors studied in the article are proprietary and are not included.
The per-ring wall-area ratios and *γ*<sub>eff</sub> values used in the article are listed in its
Supplementary Tables 1 and 2; this code and an independent implementation reproduce them to at
least six significant figures from the (undisclosed) ring geometry.

## Numerics

*R′* is found by 200-step bisection on (0, *L*); *g*(*R′*) is monotonic, so the root is unique
when it exists. For 1 − *k* < 5×10⁻⁴ the bracket is evaluated with its series (4/3) ε³ (1 − ε/2)
to avoid cancellation. Each result reports the relative residual of Eq. (3) and a first-order
error estimate of *R′*.

## Citation

Please cite the article and the archived software version (see `CITATION.cff`; Zenodo DOI
`<VERSION_DOI>`).

## License

MIT — see `LICENSE`.
