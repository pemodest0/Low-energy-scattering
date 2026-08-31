# Low-energy two-body scattering

Numerical laboratory for the low-energy limit of two-body scattering, built as
the first stage of my MSc at the São Carlos Institute of Physics (IFSC-USP).

**Universality.** At low energy only two numbers survive, the scattering length
`a` and the effective range `r0`. The same two-parameter description then holds
for systems with nothing else in common: the deuteron, bound by the strong
force at the MeV scale, and the helium-4 dimer, bound by van der Waals at the
mK scale. Converted to the same units, that is **13 orders of magnitude in
energy** and 6 in length.

**What is measured here.** Both systems, from `(a, r0)` alone, with four
unrelated potentials each. The binding energies come out within 1% of the
published values, and the four potentials agree with one another to the same
1% — the shape-independent approximation, Eq. (55) of the reference. Their
tuned parameters are of course all different; only `(a, r0)` is shared.

This work reproduces Macêdo-Lima & Madeira, *Rev. Bras. Ensino Fís.* **45**,
e20230079 (2023), and checks itself against the closed forms published there.

## Contents

```
notebooks/1_um_corpo/        one particle in a central potential
    schrodinger.py           harmonic oscillator, hydrogen, box
    one_body.ipynb

notebooks/2_dois_corpos/     the two-body laboratory
    lab.py                   potentials, Numerov solver, tabulated data
    test_lab.py              18 checks against closed forms and published tables
    two_body_scattering.ipynb

notas_teoria/
    Theory_and_Implementation.pdf every function and every constant explained
```

Each notebook needs only the `.py` file sitting next to it. Nothing else.

## Running it

```bash
pip install -r requirements.txt
jupyter lab notebooks/2_dois_corpos/two_body_scattering.ipynb   # about 2 min, Run All
python notebooks/2_dois_corpos/test_lab.py                      # the checks
```

## How it is verified

Every check compares against something the code did not compute: a closed form,
an analytic condition, or a published table.

| check | anchor | agreement |
|---|---|---|
| `a` of the square well | Eq. (80) of the article | 2e-10 |
| `r0` of the square well | Eq. (92) of the article | 6e-12 |
| bound state | `k cot(kR) = -kappa` | 3e-10 |
| `r0/R` at the poles of `a` | exact prediction, Fig. 6 | 1.0000000000 |
| node counts, 12 cases | Table 2 | exact |
| tuning, 12 cases | Tables 3 and 4 | 9 within 0.2% |
| helium dimer, 4 potentials | Motovilov *et al.*, Table I | see below |
| each Aziz minimum | its own `(rm, -eps)` | 1e-9 |
| the grid ends exactly on `R` | 120 reduced masses | exact |
| physics under `r -> r/L` | four rescalings | 1e-12 |
| every number since the last run | recorded values | 1e-12 |
| grid convergence | halving the points | no change beyond 1e-6 |

The helium row is the one that costs nothing and proves the most, because no
parameter in it was fitted here. Published parameters go in, published numbers
come out:

| potential | `a` (Å) | published | `E` (mK) | published |
|---|---|---|---|---|
| HFDHE2 | 124.65 | 124.65 | −0.83012 | −0.83012 |
| HFD-B | 88.60 | 88.50 | −1.68541 | −1.68541 |
| LM2M2 | 100.23 | 100.23 | −1.30348 | −1.30348 |
| TTY | 100.01 | 100.01 | −1.30962 | −1.30962 |

Every energy agrees to six figures. HFD-B is the exception on `a`, 0.11% out
while its energy agrees like the others; it is left visible rather than tuned
away. Reference: Motovilov, Sandhas, Sofianos & Kolganova, *Eur. Phys. J. D*
**13**, 33 (2001), Table I for the results and the Appendix for the parameters.

Two of these checks earned their place by catching something. Grid convergence
found the one real bug in the physics: on a uniform grid the Lennard-Jones `r0`
drifted instead of converging, because its hard core and its tail differ by five
orders of magnitude in scale, and a logarithmic grid fixed it. The endpoint
check found a silent one: `exp(log(R))` does not always return `R`, so for 14%
of reduced masses the last grid point fell outside the potential and the answer
was wrong in the seventh figure instead of the tenth.

## Author

Pedro Henrique Gesualdo Modesto — São Carlos Institute of Physics, USP.
Advisor: Lucas Madeira. Co-advisor: Patrícia C. M. Castilho. Supported by CAPES.
