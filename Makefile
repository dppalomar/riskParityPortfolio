clean:
	rm -v src/*.so src/*.o
	rm -v R/RcppExports.R
	rm -v src/RcppExports.cpp

build:
	Rscript .compileAttributes.R
	Rscript .roxygenize.R
	echo "useDynLib(riskParityPortfolio)" >> NAMESPACE
	echo "importFrom(Rcpp, sourceCpp)" >> NAMESPACE

install:
	R CMD INSTALL ../riskParityPortfolio

test:
	Rscript -e "devtools::test()"

all:
	make build && make install && make test
