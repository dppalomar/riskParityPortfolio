RUN apt-get -y install r-base
RUN Rscript -e "devtools::install_github('dppalomar/riskParityPortfolio')"
