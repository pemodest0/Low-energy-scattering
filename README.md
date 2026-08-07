# Low-energy quantum scattering — a validated numerical laboratory

**Run it in your browser right now, nothing to install:**

| Notebook | Language | Level | |
|---|---|---|---|
| **The Quantum Lab from scratch** — Schrödinger to Efimov, with the algebra | EN | beginner → advanced | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pemodest0/Low-energy-scattering/blob/main/notebooks/colab_quantum_lab_EN.ipynb) |
| Cold atoms from zero (short intro) | PT | beginner | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pemodest0/Low-energy-scattering/blob/main/notebooks/colab_atomos_frios_do_zero.ipynb) |
| 00 — From zero, slowly | PT | beginner | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pemodest0/Low-energy-scattering/blob/main/notebooks/00_do_zero.ipynb) |
| 01 — Guided tour of the lab | PT | intermediate | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pemodest0/Low-energy-scattering/blob/main/notebooks/01_visita_guiada.ipynb) |
| 02 — The formulas, hands on | PT | intermediate | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pemodest0/Low-energy-scattering/blob/main/notebooks/02_formulas_na_pratica.ipynb) |
| 03 — Towards the Efimov effect | PT | intermediate | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pemodest0/Low-energy-scattering/blob/main/notebooks/03_rumo_ao_efimov.ipynb) |
| 04 — Canonical Schrödinger systems | PT | intermediate | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pemodest0/Low-energy-scattering/blob/main/notebooks/04_schrodinger_canonica.ipynb) |

## What this is

A small, heavily validated Python laboratory for **low-energy two-body
quantum scattering**: it computes the s-wave scattering length $a$ and
effective range $r_0$ of microscopic potentials, reproducing every
number and figure of Macêdo-Lima & Madeira, *Rev. Bras. Ens. Fís.* **45**,
e20230079 (2023) — and then goes beyond it (Feshbach resonances of
$^{39}$K, the hyperradial Efimov ladder, canonical bound-state solvers).

Built as part of an MSc project at the São Carlos Institute of Physics
(IFSC-USP), advised by Lucas Madeira, on **Efimov physics under
dimensional crossover** in ultracold gases.

## Highlights

- **Two integrators** (central difference and Numerov) for
  $u''(r) = 2V(r)u(r)$, with $a$ from the log-derivative match and
  $r_0$ from the effective-range integral (trapezoid *and* Simpson).
- **Four potentials**: spherical well, modified Pöschl-Teller, Gaussian,
  Lennard-Jones — plus Aziz HFD-B (He–He), soft/hard spheres.
- **Parameter tuning** (nested loops on $1/a$ and $r_0$) reproducing the
  article's Tables 3–4 to ~0.1–1%.
- **25 automated tests**, each anchored to an analytical or published
  result — never to "what the code gave yesterday".
- **Four-layer validation**: closed formulas (errors $10^{-5}$–$10^{-10}$),
  the article's tables, cross-methods (variable-phase equation agrees to
  $10^{-13}$), and external benchmarks (Jeszenszki *et al.* Gaussian
  $a(v)$ to ~$10^{-6}$).
- **$^{39}$K Feshbach module** (D'Errico 2007 / Zaccanti 2009 data):
  $a(B)$, van der Waals length, predicted first Efimov loss feature at
  $B \approx 403.5$ G.
- **Efimov ladder from scratch**: $s_0 = 1.0062378$ solved from the
  transcendental equation; energy ratios converge to $e^{2\pi/s_0} = 515.03$.
- **Interactive app** (Streamlit): 7 stations, each with
  Simulate / Theory / Code / Results tabs.

### Findings worth knowing

1. **Factor-of-2 convention in the article's Eq. (121)**: the published
   Table-4 Lennard-Jones constants only reproduce $(a, r_0)$ with
   $V = \hbar^2/(2m_r)\,[C_{12}/r^{12} - C_6/r^6]$.
2. **Numerov loses its high order on discontinuous potentials** (the
   spherical well's edge): central difference with an edge-inside
   convention is $O(\Delta r^2)$ there, Numerov degrades to $O(\Delta r)$.
3. **The helium dimer lives on a knife's edge**: the classic de Boer
   Lennard-Jones does *not* bind it ($a = -178$ Å) while Aziz HFD-B does
   ($a = +88.4$ Å, $E = -1.69$ mK) — same atom, opposite answers.

## Quickstart (local)

```bash
pip install -r requirements.txt
python -m pytest tests/ -q          # 25 tests, ~15 s
python main.py --sem-ajustes        # regenerate figures + tables (~1 min)
streamlit run app/Inicio.py         # the interactive lab
```

## Layout

```
src/          physics core (potentials, solvers, a & r0, tuning, Feshbach, Efimov)
tests/        the 25-test safety net
notebooks/    the teaching track (badges above)
app/          Streamlit lab (7 stations)
results/      CSV outputs vs. the article
references/   external benchmark CSVs
```

Code comments are currently in Portuguese (being translated
progressively); all physics, tests and this documentation are in English.

## Cite

See `CITATION.cff`. Please also cite the article this laboratory
reproduces: M. Macêdo-Lima and L. Madeira, *Rev. Bras. Ens. Fís.* **45**,
e20230079 (2023), [doi:10.1590/1806-9126-RBEF-2023-0079](https://doi.org/10.1590/1806-9126-RBEF-2023-0079).

## License

MIT — see `LICENSE`.
