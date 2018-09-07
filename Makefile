R_OPTS=--no-save --no-restore --no-init-file --no-site-file

rqtl2_paper.pdf: rqtl2_paper.tex rqtl2_paper.bib Figs/fig1.png
	pdflatex rqtl2_paper
	bibtex rqtl2_paper
	pdflatex rqtl2_paper
	pdflatex rqtl2_paper

Figs/fig1.png: R/fig1.R R/Perms/operm_add.rds R/Perms/operm_full.rds R/Perms/operm_snps.rds
	cd R;R $(R_OPTS) -e "source('$(<F)', echo=TRUE)"

clean:
	rm rqtl2_paper.aux rqtl2_paper.bbl rqtl2_paper.blg rqtl2_paper.log rqtl2_paper.out rqtl2_paper.pdf Figs/fig1.png

clean_cache: clean
	rm R/cache/*.rds R/cache/*.zip R/cache/*.RData

web: rqtl2_paper.pdf
	scp rqtl2_paper.pdf broman-10.biostat.wisc.edu:Website/publications/
