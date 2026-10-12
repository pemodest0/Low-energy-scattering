![The electron cloud of a hydrogen atom: the computed probability density of an eigenstate, with its nodal surfaces](docs/cover.png)

# Low-energy scattering

Numerical laboratory for cold atoms at low energy, built during my MSc at the
São Carlos Institute of Physics (IFSC-USP). Two stages so far: the two-body
problem, where only two numbers survive, and the harmonic trap that holds the
atoms while you look at them.

## Stage one — two-body universality

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

This part reproduces Macêdo-Lima & Madeira, *Rev. Bras. Ensino Fís.* **45**,
e20230079 (2023), and checks itself against the closed forms published there.

## Stage two — the harmonic trap and the road to two dimensions

A trapped gas never sees the free-space problem alone: it sees the trap too.
Three frequencies set the shape of the cloud, and two different means of those
three show up in different places — the geometric mean in the critical
temperature and the counting of states, the quadratic mean in the symmetric
sector of the collective excitations. Squeezing one direction hard enough with
an optical lattice turns the gas two-dimensional, which is a different physical
world: no condensate in the usual sense, and superfluidity arriving through the
Berezinskii-Kosterlitz-Thouless transition instead.

`harmonic_trap/accordion.ipynb` is the design calculation for that squeeze. You
give it the beam parameters and it returns the lattice period, the axial
frequency from the exact band structure (and how far the usual harmonic formula
is off), the three trap frequencies with the resulting cloud shape, whether the
gas is actually two-dimensional, and the effective in-plane interaction. It is
self-contained — `numpy`, `scipy`, `matplotlib`, nothing else — and only one
cell needs editing.

It also says which resonance to look for when calibrating the lattice depth by
amplitude modulation, and why that peak is not the trap frequency: modulating
the depth keeps the potential symmetric, so the drive only connects states of
the same parity and the resonance is the transition to the second band, not the
first. Reading that peak as the trap frequency is wrong by nearly a factor of
two.

Everything inside that needs confirming against a source is flagged where it is
used: the modelled polarizability (the notebook prints its own static limit so
you can compare it with the published value for potassium), the Feshbach
parameters, the constant in the Berezinskii-Kosterlitz-Thouless criterion, and
the range in which the weak-confinement form of the two-dimensional coupling
still holds.

## Contents

```
notebooks/1_one_body/        one particle in a central potential
    schrodinger.py           harmonic oscillator, hydrogen, box
    one_body.ipynb

notebooks/2_two_body/        the two-body laboratory
    lab.py                   potentials, Numerov solver, tabulated data
    test_lab.py              18 checks against closed forms and published tables
    two_body_scattering.ipynb

harmonic_trap/
    accordion.ipynb          accordion lattice design: period, frequencies,
                             cloud shape, two-dimensionality, interaction

theory_notes/
    Theory_and_Implementation.pdf every function and every constant explained
```

Each notebook needs only the `.py` file sitting next to it, or nothing at all.

## Running it

```bash
pip install -r requirements.txt
jupyter lab notebooks/2_two_body/two_body_scattering.ipynb   # about 2 min, Run All
python notebooks/2_two_body/test_lab.py                      # the checks
jupyter lab harmonic_trap/accordion.ipynb                       # seconds, Run All
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
| lattice band structure | `2 sqrt(s)` as `s` grows | asymptotic, see notebook |

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

The lattice row is of the same kind. The harmonic formula for the axial
frequency is an approximation, and the notebook plots how far off it is against
the exact band structure: a few percent for a deep lattice, tens of percent
below about ten recoil energies. Calibrating a lattice from that formula without
checking is a quiet way to be wrong.

## References

Beyond the two already named above, the trap and excitation part rests on:

- S. Stringari, *Collective excitations of a trapped Bose-condensed gas*,
  Phys. Rev. Lett. **77**, 2360 (1996); arXiv:cond-mat/9603126.
- F. Zambelli and S. Stringari, *Finite size effects on the collective
  oscillations of a trapped Bose gas*; arXiv:cond-mat/0107464.

## The cover

The cover is computed, not drawn and not photographed, so the image carries no
third-party rights. It is the electron cloud of a hydrogen atom: the
probability density of a single eigenstate in the plane containing the symmetry
axis, from the exact radial function and the exact angular function. The dark
surfaces crossing it are the real nodes of that state, where the wavefunction
changes sign and the electron is never found.

## Author

Pedro Henrique Gesualdo Modesto — São Carlos Institute of Physics, USP.
Advisor: Lucas Madeira. Co-advisor: Patrícia C. M. Castilho. Supported by CAPES.
