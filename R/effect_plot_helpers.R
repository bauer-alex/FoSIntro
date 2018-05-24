
#' Estimate posterior simulation prediction intervals for Gamma regression
#' 
#' Prediction intervals for a function-on-scalar Gamma regression model
#' are estimated using the 'posterior simulation' approach for the parameter vector,
#' e.g. outlined here: https://stat.ethz.ch/pipermail/r-help/2011-April/275632.html.
#' This approach is more generally applicable to generalized regression models,
#' however the implementation is currently only focused on Gamma regression.
#' 
#' @param model Function-on-scalar regression model fitted with \code{\link[refund]{pffr}}
#' @param newdata The new data.frame for which the predictions are calculated
#' @param alpha Alpha level, defaults to 0.05, i.e. 95\% prediction intervals
#' @param n.rep number of replications for the posterior simulation
#' @param seed Optional seed for reproducibility
#' @return \code{data.frame} with two columns \code{pi_lower} and \code{pi_upper}
#' @export
pi_gamma <- function(model, newdata, alpha, n.rep = 10000, seed = NULL) {
  if (!is.null(seed))
    set.seed(seed)

    ## extract parameter estiamtes and cov matrix
  beta <- model$coefficients
  Vb <- vcov(model)
  
  ## simulate replicate beta vectors from posterior
  Cv <- chol(Vb)
  nb <- length(beta)
  br <- t(Cv) %*% matrix(rnorm(n.rep*nb),nb,n.rep) + beta
  
  ## turn these into replicate linear predictors
  Xp <- refund:::predict.pffr(model, newdata = newdata, type="lpmatrix")
  lp <- Xp%*%br
  fv <- exp(lp) ## ... finally, replicate expected value vectors
  
  ## now simulate from Gamma deviates with mean as in fv and estimated scale
  scale <- model$scale
  yr <- matrix(rgamma(fv*0,shape=1/scale, scale=fv*scale), nrow(fv), ncol(fv))
  
  ## compute the prediction interval
  PI <- apply(yr,1,quantile,prob=c(alpha/2, 1-(alpha/2)))
  data.frame("pi_lower" = PI[1,],
             "pi_upper" = PI[2,])
}


#' Create a grid of multiple ggplots sharing the same legend
#'
#' Function to create a grid of multiple \code{\link[ggplot2]{ggplot}}s sharing
#' the same legend.
#' Note: This function is based on a code example shown on the
#' \href{https://github.com/tidyverse/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs}{tidyverse GitHub page}
#' (last access: 25.10.2017)
#'
#' @param ... Multiple \code{\link[ggplot2]{ggplot}} objects
#' @param nrow,ncol Parameters specifying the grid to arrange the ggplots in.
#' Defaults to one row and multiple columns for the plots.
#' @param legend_position The legend can be plotted at the \code{"bottom"} or on
#' the \code{"right"} of all plots
#' @param main,sub Text over and under the whole graphic, respectively
#' @param outer_fontsize Fontsize of \code{main} and \code{sub}. Defaults to 15.
#' @param newpage If TRUE, \code{\link[grid]{grid.newpage}} is called before
#' plot creation. This is necessary for the plot window in R, but should be
#' suppressed if one exports the plot, e.g. using \code{\link{pdf}}
#' @return Silently returns the plot object
#' @importFrom grid textGrob gpar unit.c grid.newpage grid.draw
#' @importFrom gridExtra arrangeGrob
#' @import ggplot2
grid_arrange_shared_legend <- function(..., nrow = 1, ncol = length(list(...)),
                                       legend_position = c("bottom", "right"), main = NULL, sub = NULL,
                                       outer_fontsize = 15, newpage = TRUE) {
  plots <- list(...)
  legend_position <- match.arg(legend_position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = legend_position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(legend_position,
                     "bottom" = gridExtra::arrangeGrob(do.call(gridExtra::arrangeGrob, gl), legend,
                                                       ncol = 1,
                                                       top = grid::textGrob(main, gp=grid::gpar(fontsize = outer_fontsize)),
                                                       bottom = grid::textGrob(sub, gp=grid::gpar(fontsize = outer_fontsize)),
                                                       heights = grid::unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = gridExtra::arrangeGrob(do.call(gridExtra::arrangeGrob, gl), legend,
                                                      ncol = 2,
                                                      top = grid::textGrob(main, gp=grid::gpar(fontsize = outer_fontsize)),
                                                      bottom = grid::textGrob(sub, gp=grid::gpar(fontsize = outer_fontsize)),
                                                      widths = grid::unit.c(unit(1, "npc") - lwidth, lwidth)))
  if (newpage)
    grid::grid.newpage()
  grid::grid.draw(combined)
  invisible(combined) # silently return the plot object
}


#' Extract returned values of plot.gam() while suppressing creation of the plot
#'
#' Function to extract the values returned of \code{\link[refund]{plot.pffr}} while
#' suppressing creation of the plot.
#' 
#' @param model Function-on-scalar model fitted with \code{\link[refund]{pffr}}
#' @importFrom grDevices png dev.off
#' @importFrom mgcv plot.gam
no.plot <- function (model) {
  png("temp.xyz")
  # plot.gam returns all terms by default, select and rug are set only to decrease evaluation time
  plot.df <- mgcv::plot.gam(model, select = 1, rug = FALSE)
  dev.off()
  unlink("temp.xyz", recursive = TRUE)
  # delete 'raw' elements as they are very large but not necessary for plotting the effects
  plot.df <- lapply(plot.df, function(x) {
    x$raw <- NULL
    x
  })
  return(invisible(plot.df))
}
