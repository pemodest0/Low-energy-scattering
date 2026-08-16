# Build the manuscript.  Just type:  make
#
# The tables and figures are NOT edited by hand — they are written by
# notebooks/2_dois_corpos/two_body_scattering.ipynb. Run the notebook first
# if the numbers changed.

pdf:            ## build main.pdf (runs BibTeX automatically)
	latexmk -pdf -interaction=nonstopmode main.tex

clean:          ## delete intermediates, keep the PDF
	latexmk -c

wipe:           ## delete everything including the PDF
	latexmk -C

.PHONY: pdf clean wipe
