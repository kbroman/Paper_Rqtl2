R_OPTS=--no-save --no-restore --no-init-file --no-site-file

all: rqtl2_paper.pdf

rqtl2_paper.pdf: rqtl2_paper.tex rqtl2_paper.bib Figs/fig1.png Figs/fig2.pdf
	pdflatex rqtl2_paper
	bibtex rqtl2_paper
	pdflatex rqtl2_paper
	pdflatex rqtl2_paper

Figs/fig1.png: R/fig1.R R/Perms/operm_add.rds R/Perms/operm_full.rds R/Perms/operm_snps.rds
	cd R;R $(R_OPTS) -e "source('$(<F)', echo=TRUE)"

Figs/fig2.pdf: R/fig2.R GnanResults/gnan_results.rds
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"

GnanResults/gnan_results.rds: GnanResults/combine_results.R
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"
