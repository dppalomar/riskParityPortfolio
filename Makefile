build:
	Rscript .roxygenize.R

install:
	R CMD INSTALL ../riskParityPortfolio

test:
	Rscript -e "devtools::test()"
