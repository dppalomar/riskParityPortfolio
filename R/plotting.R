#' @title Create portfolio barplots with the capital allocation and the risk allocation
#' 
#' @description Creates a barplot on top with the portfolio capital allocation and another
#' at the bottom with the risk contribution allocation whose profile is the target of the
#' risk parity portfolio design with \code{\link{riskParityPortfolio}}.
#' By default the plot is based on the package \code{ggplot2}, but the user
#' can also specify a simple base plot.
#' 
#' @param w Vector or matrix containing the portfolio(s) weights. 
#'          For multiple portfolios, they should be columnwise and named for the legend.
#' @param Sigma Covariance matrix of the assets.
#' @param type Type of plot. Valid options: \code{"ggplot2", "simple"}. Default is 
#'             \code{"ggplot2"} (the packages \code{ggplot2} and \code{gridExtra} must be installed).
#' @param colors Vector of colors for the portfolios (default uses grDevices::topo.colors()).
#' @examples
#' library(riskParityPortfolio)
#' 
#' # generate random covariance matrix
#' set.seed(42)
#' N <- 10
#' V <- matrix(rnorm(N^2), nrow = N)
#' Sigma <- cov(V)
#' 
#' # generate random portfolio vectors
#' w_single <- runif(N)
#' w_single <- w_single/sum(w_single)  # normalize
#' names(w_single) <- LETTERS[1:N]
#' 
#' w_multiple <- matrix(runif(4*N), ncol = 4)
#' w_multiple <- sweep(w_multiple, MARGIN = 2, STATS = colSums(w_multiple), FUN = "/")  # normalize each column
#' rownames(w_multiple) <- LETTERS[1:N]
#' 
#' # plot
#' barplotPortfolioRisk(w_single, Sigma)
#' barplotPortfolioRisk(w_multiple, Sigma)
#' barplotPortfolioRisk(w_multiple, Sigma, colors = viridisLite::viridis(ncol(w)))
#' 
#' @import ggplot2
#' @import gridExtra
#' @export
barplotPortfolioRisk <- function(w, Sigma, type = c("ggplot2", "simple"), colors = topo.colors(ncol(w))) {
  w <- as.matrix(w)
  if (is.null(colnames(w)))
    colnames(w) <- paste0("portf-", 1:ncol(w))
  if (is.null(rownames(w)))
    rownames(w) <- paste0("stock", 1:nrow(w))  
  RRC <- w * (Sigma %*% w)
  RRC <- sweep(RRC, MARGIN = 2, STATS = colSums(RRC), FUN = "/")  # normalize each column
  
  # plot
  switch(match.arg(type),
         "simple" = {
           old_par <- par(mfrow=c(2,1))
           barplot(t(w), col = colors,
                   beside = ncol(w) > 1, legend = colnames(w),
                   main = "Portfolio weight allocation", ylab = "capital")
           barplot(t(RRC), col = colors,
                   beside = ncol(w) > 1,
                   main = "Relative risk contribution", xlab = "stocks", ylab = "risk")
           par(old_par)
         },
         "ggplot2" = {
           # stock <- portfolio <- value <- NULL  # ugly hack to deal with CRAN note
           p1 <- ggplot(data = melt_portf_matrix(w), aes(x = stock, y = value)) + 
             geom_bar(aes(fill = portfolio), color = "black", stat = "identity", position = "dodge", width = 0.8) + 
             scale_fill_manual(values = colors) +
             labs(title = "Portfolio weight allocation", x = "stocks", y = "capital") + 
             theme(legend.title = element_blank())
           p2 <- ggplot(data = melt_portf_matrix(RRC), aes(x = stock, y = value)) + 
             geom_bar(aes(fill = portfolio), color = "black", stat = "identity", position = "dodge", width = 0.8) + 
             scale_fill_manual(values = colors) +
             labs(title = "Relative risk contribution", x = "stocks", y = "risk") + 
             theme(legend.title = element_blank())
           # theme(plot.title = element_text(size = 16),
           #       axis.title = element_text(size = 14, color = '#555555'),
           #       axis.text = element_text(size = 11),
           #       legend.position = c(0.9, 0.86), #c(0.905, 0.88), 
           #       legend.title = element_blank(), 
           #       legend.spacing.x = unit(4, 'pt'),
           #       legend.text = element_text(size = 13))
           if (ncol(w) == 1) {  # remove legend if only one portfolio
             p1 <- p1 + theme(legend.position = "none")
             p2 <- p2 + theme(legend.position = "none")
           }           
           grid.arrange(p1, p2, nrow = 2)
         },
         stop("Barplot type unknown."))
}


# Equivalent to reshape2::melt(data, varnames = c("stock", "portfolio"))
melt_portf_matrix <- function (data, na.rm = FALSE, as.is = FALSE) 
{
  varnames = c("stock", "portfolio")
  value.name = "value"
  
  var.convert <- function(x) {
    if (!is.character(x)) 
      return(x)
    x <- type.convert(x, as.is = TRUE)
    if (!is.character(x)) 
      return(x)
    factor(x, levels = unique(x))
  }
  dn <- dimnames(data)
  names(dn) <- varnames
  if (!as.is) {
    dn <- lapply(dn, var.convert)
  }
  labels <- expand.grid(dn, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  if (na.rm) {
    missing <- is.na(data)
    data <- data[!missing]
    labels <- labels[!missing, ]
  }
  value_df <- setNames(data.frame(as.vector(data)), value.name)
  cbind(labels, value_df)
}
