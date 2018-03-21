
#' Residual plot vs fitted values for \code{\link[refund]{pffr}} models
#'
#' Plot the residuals of a function-on-scalar model fitted with 
#' \code{\link[refund]{pffr}} against its fitted values.
#'
#' @param model Function-on-scalar model fitted with \code{\link[refund]{pffr}}
#' @param base_size Size of plot elements. see \code{\link[ggplot2]{theme_bw}}.
#' Defaults to 11.
#' @param xlab,ylab,main Optional plot annotation
#' @param legend_limits,legend_breaks Optional arguments, passed as \code{limits}
#' and \code{breaks} to \code{\link[ggplot2]{scale_fill_gradientn}}.
#' @param ... Further arguments passed to \code{\link[ggplot2]{theme}}
#' @import ggplot2
#' @import refund
#' @export
plot_resVSfitted <- function(model, base_size = 11, xlab = "fitted values",
                             ylab = "residuals", main = "Residuals vs fitted values\nHeatmap of binned points",
                             legend_limits = NULL, legend_breaks = waiver(), ...) {
  plot_data <- data.frame("fitted" = as.vector(refund:::fitted.pffr(model)),
                          "residuals" = as.vector(refund:::residuals.pffr(model, type = "response")))
  ggplot(plot_data, aes_string(x="fitted", y="residuals")) + stat_binhex(color = gray(0.7), bins = 30) +
    scale_fill_gradientn(colours=c(gray(0.9),"darkblue"),name = "Frequency",
                         limits = legend_limits, breaks = legend_breaks) +
    xlab(xlab) + ylab(ylab) + ggtitle(main) +
    theme_bw(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5), ...)
}


#' Residual plot vs functional domain for \code{\link[refund]{pffr}} models
#'
#' Plot the residuals of a function-on-scalar model fitted with 
#' \code{\link[refund]{pffr}} against its functional domain (e.g. time).
#'
#' @inheritParams plot_resVSfitted
#' @import ggplot2
#' @import refund
#' @export
plot_resVSyindex <- function(model, base_size = 11, xlab = "yindex",
                             ylab = "residuals", main = "Residuals vs yindex\nHeatmap of binned points",
                             legend_limits = NULL, legend_breaks = waiver(), ...) {
  yindex <- model$pffr$yind
  resid_matrix <- refund:::residuals.pffr(model, type = "response")
  plot_data <- data.frame("yindex" = rep(yindex, each = nrow(resid_matrix)),
                          "residuals" = as.vector(resid_matrix))
  ggplot(plot_data, aes_string(x="yindex", y="residuals")) + stat_binhex(color = gray(0.7), bins = 30) +
    scale_fill_gradientn(colours=c(gray(0.9),"darkblue"), name = "Frequency",
                         limits = legend_limits, breaks = legend_breaks) +
    xlab(xlab) + ylab(ylab) + ggtitle(main) +
    theme_bw(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5), ...)
}

#' Residual autocovariance plot for \code{\link[refund]{pffr}} models
#' 
#' Plot the autocovariance of the residuals over the functional domain (e.g. time) for 
#' a function-on-scalar model fitted with \code{\link[refund]{pffr}}.
#' 
#' @param model Function-on-scalar model fitted with \code{\link[refund]{pffr}}
#' @param base_size Size of plot elements. see \code{\link[ggplot2]{theme_bw}}.
#' Defaults to 11.
#' @param legend.position,legend.key.width Arguments for \code{\link[ggplot2]{theme}}
#' @param cov.use see \code{use} argument in \code{\link[stats]{cov}}
#' @param ... Additional arguments passed to \code{\link[ggplot2]{scale_fill_gradient}}
#' @importFrom stats cov
#' @importFrom tidyr gather
#' @import ggplot2
#' @import refund
#' @export
plot_residAutocov <- function(model, base_size = 11, legend.position = "bottom",
                              legend.key.width = unit(2,"lines"), cov.use = "everything", ...) {
  yindex <- model$pffr$yind
  dat <- data.frame(time = yindex, cov(refund:::residuals.pffr(model), use = cov.use))
  colnames(dat)[-1] <- yindex
  dat <-  tidyr::gather(dat, key = "time2", value = "covariance", 2:ncol(dat))
  dat$time2 <- as.numeric(dat$time2)
  ggplot(dat, aes_string(x="time", y="time2")) +
    geom_tile(aes_string(fill = "covariance")) +
    xlab("time [s]") + ylab("time [s]") +
    scale_fill_gradient(low = "grey90", high = "darkblue", ...) +
    theme_bw(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = legend.position,
          legend.key.width = legend.key.width) +
    ggtitle("Autocovariance of residuals")
}


#' Data preparation for \code{\link{plot_residsVSxy}}
#' 
#' Prepare the data set containing the mean residuals per measurement location
#' for the following call of \code{\link{plot_residsVSxy}}.
#' Note: This function can only be used for situations, where measurements
#' took place at a set of fixed locations in 2D, which all have multiple measurements.
#' Situations where each observations has a unique metric value on the x and/or y 
#' axis are not supported. Neither are other spatial dimensions than 2D.
#' 
#' @param model Function-on-scalar model fitted with \code{\link[refund]{pffr}}
#' @param data Data for which the residuals of the model should be looked at
#' @param yvar name of the response variable in \code{model}
#' @param xCoord_var,yCoord_var name of the variables containing location information
#' for the x and y axis of the measurements
#' @param ... Additional arguments passed to \code{\link[refund]{predict.pffr}}
#' @import dplyr
#' @import refund
#' @export
prepareData_residsVSxy <- function(model, data, yvar, xCoord_var, yCoord_var, ...) {
  colnames(data)[colnames(data) == xCoord_var] <- "x"
  colnames(data)[colnames(data) == yCoord_var] <- "y"
  data$location <- factor(paste0(data$x,"_",data$y))
  y <- data[,yvar]
  dat <- data[,colnames(data) != yvar] # the y variable has to be excluded to call predict()
  dat$fitted <- refund:::predict.pffr(model, newdata = dat, type = "response", ...)
  dat$residuals <- y - dat$fitted
  # Create a dataset with 1 row per location
  dat_res <- dat[,!(colnames(dat) %in% c("residuals","fitted"))] %>%
    distinct(location, .keep_all = TRUE)
  # Calculate means of residuals per location
  dat_res$resids_mean <- unlist(sapply(dat_res$location,
                                       FUN = function(l) mean(unlist(dat$residuals[dat$location == l,]))))
  # Sort residuals close to zero to top so that they get plotted in the background
  ord <- order(abs(dat_res$resids_mean))
  dat_res <- dat_res[ord,]
  
  dat_res
}



#' Residual plot vs space for \code{\link[refund]{pffr}} models
#'
#' Plot the residuals of a function-on-scalar model fitted with 
#' \code{\link[refund]{pffr}}, which includes a spatial component, against the latter.
#' Note: This function can only be used for situations, where measurements
#' took place at a set of fixed locations in 2D, which all have multiple measurements.
#' Situations where each observations has a unique metric value on the x and/or y 
#' axis are not supported. Neither are other spatial dimensions than 2D.
#' 
#' @param model Function-on-scalar model fitted with \code{\link[refund]{pffr}}
#' @param dat_xy Data set precalculated with \code{\link{prepareData_residsVSxy}},
#' containing information on the mean residuals per measurement location.
#' @param breaks breaks to be used for the colors of the residual categories.
#' E.g. \code{c(-10,0,10)} to show the residuals categorized in \code{c("(-10,0]","(0,10]")}
#' @param labels_cut Optional labels for the categories defined using \code{breaks}.
#' @param base_size see \code{\link[ggplot2]{theme_bw}}
#' @param legend_rows,legend_cols Optional parameters specifying the grid of
#' legend elements in \code{\link[ggplot2]{guide_legend}}
#' @param mark_location Optional row number \code{dat_xy} for which the location
#' should be marked by a special black dot in the plot
#' @param mark_location_size Size of the black dot for \code{mark_location}
#' @param ... Further arguments passed to \code{\link[ggplot2]{theme}}
#' @importFrom grDevices gray
#' @importFrom RColorBrewer brewer.pal
#' @import ggplot2
#' @export
plot_residsVSxy <- function(model, dat_xy, breaks, labels_cut = NULL, base_size = 11,
                            legend_rows = NULL, legend_cols = NULL, mark_location = NULL,
                            mark_location_size = 5, ...) {
  if (missing(breaks))
    stop("Please specify 'breaks' for the residual legend!")
  if (any(abs(dat_xy$resids_mean) > 1))
    stop("At least one mean of a seismogram has an absolute value bigger than one -> Please edit the breaks accordingly!")
  # Plot
  resids_kat <- cut(dat_xy$resids_mean, breaks = breaks, labels = labels_cut)
  colors <- rev(RColorBrewer::brewer.pal(length(levels(resids_kat)), "RdBu"))
  index_neutralCat <- which(sapply(1:(length(breaks)-1), function(i) (0 > breaks[i]) && (0 < breaks[i+1])))
  colors[index_neutralCat] <- gray(0.9) # change neutral color to make it visible on b/w print
  colors[1:(index_neutralCat-1)] <- paste0(colors[1:(index_neutralCat-1)],"99") # add slight transparency of points
  shapes <- c(rep(15,index_neutralCat), rep(18, length(colors) - index_neutralCat))
  sizes <- rep(2,length(colors))
  sizes[index_neutralCat] <- 3
  
  names(colors) <- labels_cut
  names(shapes) <- labels_cut
  names(sizes) <- labels_cut

  gg <- ggplot(dat_xy, aes_string(x="x", y="y", color="resids_kat")) +
    geom_point(mapping = aes_string(shape="resids_kat", size="resids_kat")) + theme_bw(base_size = base_size) +
    theme(legend.key.height = unit(2,"line"), plot.title = element_text(hjust = 0.5),
          axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), ...) +
    scale_colour_manual(name="Residual mean", values=colors) +
    scale_shape_manual(name="Residual mean", values=shapes) +
    scale_size_manual(name="Residual mean", values=sizes) +
    ggtitle("Mean residuals over space")
  if (!is.null(mark_location))
    gg <- gg + geom_point(data = dat_xy[mark_location,c("x","y")],
                          color = "black", size = mark_location_size)
  if(!is.null(legend_rows))
    gg <- gg + guides(colour = guide_legend(nrow = legend_rows, byrow = TRUE, override.aes = list(size=6)))
  if(!is.null(legend_cols))
    gg <- gg + guides(colour = guide_legend(ncol = legend_cols, byrow = TRUE, override.aes = list(size=6)))
  gg
}
