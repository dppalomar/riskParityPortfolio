clean:
	rm -v src/*.so src/*.o
	rm -v R/RcppExports.R
	rm -v src/RcppExports.cpp

build:
	Rscript .compileAttributes.R
	Rscript .roxygenize.R

install:
	R CMD INSTALL ../riskParityPortfolio

test:
	Rscript -e "devtools::test()"

all:
	make build && make install && make test
