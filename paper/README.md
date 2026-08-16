# paper/ — manuscript draft

## Build

```bash
make          # -> main.pdf
```

## The one rule

**No number is ever typed in `main.tex`.**

```
notebooks/2_dois_corpos/two_body_scattering.ipynb
        │
        ├── writes ──> paper/tables/*.tex     ──┐
        └── writes ──> paper/figures/*.pdf    ──┤
                                                │
                       main.tex \input{}s them ─┘
```

Change the physics, re-run the notebook, type `make`. The paper updates itself.
If a number in the PDF looks wrong, the fix is in the notebook — never here.

## Files

| File | What it is | Edit by hand? |
|---|---|---|
| `main.tex` | the text | **yes — this is the only file you write** |
| `references.bib` | bibliography | yes, add entries |
| `tables/*.tex` | generated tables | no |
| `tables/*.csv` | same data, for inspection | no |
| `figures/*.pdf` | generated figures | no |
| `Makefile` | build commands | rarely |

## The five LaTeX things that actually matter

1. `\input{tables/deuteron.tex}` drops a whole table in. That is the whole trick.
2. `\label{tab:deuteron}` + `\ref{tab:deuteron}` — never write "Table 3" by hand;
   LaTeX renumbers when you reorder.
3. `\cite{hackenburg2006}` pulls from `references.bib`. `latexmk` runs BibTeX for you.
4. Math goes between `$...$` inline or `\[...\]` displayed.
5. **Read the FIRST error in the log, not the last.** One missing brace produces
   fifty errors; only the first one is real.
