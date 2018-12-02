FROM rocker/verse:latest
RUN Rscript -e "devtools::install_github('dppalomar/riskParityPortfolio')"
