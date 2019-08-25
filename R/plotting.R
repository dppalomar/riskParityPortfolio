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
#' @param colors Vector of colors for the portfolios (optional).
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
#' w_multiple <- sweep(w_multiple,  # normalize each column
#'                     MARGIN = 2, 
#'                     STATS = colSums(w_multiple), FUN = "/")
#' rownames(w_multiple) <- LETTERS[1:N]
#' 
#' # plot
#' barplotPortfolioRisk(w_single, Sigma)
#' barplotPortfolioRisk(w_multiple, Sigma)
#' barplotPortfolioRisk(w_multiple, Sigma, colors = viridisLite::viridis(4))
#' barplotPortfolioRisk(w_multiple, Sigma) + ggplot2::scale_fill_viridis_d()
#' 
#' @author Daniel P. Palomar and Ze Vinicius
#' 
#' @export
barplotPortfolioRisk <- function(w, Sigma, type = c("ggplot2", "simple"), colors = NULL) {
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
           if (is.null(colors))
             colors <- grDevices::topo.colors(ncol(w))
           old_par <- graphics::par(mfrow=c(2,1))
           graphics::barplot(t(w), col = colors,
                             beside = ncol(w) > 1, legend = colnames(w),
                             main = "Portfolio weight allocation", ylab = "capital")
           graphics::barplot(t(RRC), col = colors,
                             beside = ncol(w) > 1,
                             main = "Relative risk contribution", xlab = "stocks", ylab = "risk")
           graphics::par(old_par)
         },
         "ggplot2" = {
           if (!requireNamespace("ggplot2", quietly = TRUE)) 
             stop("Please install package \"ggplot2\" or choose another plot type", call. = FALSE)
           molten_portf_matrix <- rbind(cbind(melt_portf_matrix(w), type = "w"),
                                        cbind(melt_portf_matrix(RRC), type = "RRC"))
           labels <- c(w = "Weight allocation", RRC = "Relative risk contribution")
           p <- ggplot2::ggplot(molten_portf_matrix, ggplot2::aes_string(x = "stock", y = "value")) + 
             ggplot2::geom_bar(ggplot2::aes_string(fill = "portfolio"), color = "black", stat = "identity", position = "dodge", width = 0.8) +
             ggplot2::facet_wrap(~ type, ncol = 1, scales = "free", labeller = ggplot2::labeller(type = labels)) +
             ggplot2::labs(title = "Portfolio capital and risk distribution", x = "stocks", y = NULL) + 
             ggplot2::theme(legend.title = ggplot2::element_blank()) +
             ggplot2::theme(strip.text = ggplot2::element_text(size = 11))
           if (!is.null(colors))
             p <- p + ggplot2::scale_fill_manual(values = colors)
           if (ncol(w) == 1)  # remove legend if only one portfolio
             p <- p + ggplot2::theme(legend.position = "none")
           p
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
    x <- utils::type.convert(x, as.is = TRUE)
    if (!is.character(x)) 
      return(x)
    factor(x, levels = unique(x))
  }
  dn <- dimnames(data)
  names(dn) <- varnames
  if (!as.is)
    dn <- lapply(dn, var.convert)
  labels <- expand.grid(dn, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  if (na.rm) {
    missing <- is.na(data)
    data <- data[!missing]
    labels <- labels[!missing, ]
  }
  value_df <- stats::setNames(data.frame(as.vector(data)), value.name)
  cbind(labels, value_df)
}
