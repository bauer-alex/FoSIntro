
#' Plot 1D smooth effects for \code{\link[refund]{pffr}} models
#'
#' Plots 1D smooth effects for a function-on-scalar model fitted with 
#' \code{\link[refund]{pffr}}. Additionally, Marra \& Wood (2012) or bootstrapped-based confidence
#' intervals (CIs) can be plotted, specifying the \code{plot_type} argument.
#'
#' @param model Function-on-scalar model fitted with \code{\link[refund]{pffr}}
#' @param plot_type Either 1 (plotting CIs after Marra \& Wood (2012)) or 2
#' (bootstrap-based CIs, by specyfing \code{CIs_list}). Defaults to 1.
#' @param plot_ci If \code{TRUE} CIs are plotted. Only used if \code{plot_type = 1}.
#' @param alpha \code{(1-alpha)} CIs are calculated. The default 0.05 leads to 
#' 95\% CIs. Only used if \code{plot_type = 1}.
#' @param CIs_list Output of function \code{\link{calc_bootstrapCIs}}.
#' Only necessary if \code{plot_type = 2}.
#' @param select Index of smooth term to be plotted. If a vector of length two
#' is passed, then the difference of the second smooth minus the first smooth will
#' be plotted. Note that CIs from \code{plot_type == 1} are suppressed in this case.
#' However, using the \code{CIs_list} argument one can still plot boostrapped
#' CIs for such a difference of smooths, by passing a list with only one CI element.
#' Such a CI element can be obtained by using \code{\link{calc_singleBootstrapCI}},
#' using the \code{select} argument.
#' @param xlab,ylab,ylim,main Optional parameters for overwriting the default
#' plot parameters
#' @param effect_label Used in the main title of the plot as
#' '\code{effect_label} with global 95\% CI'. Only used if \code{main} is not specified.
#' @param shift Numeric shift of the plotted effect.
#' @param base_size Size of plot elements. see \code{\link[ggplot2:ggtheme]{theme_bw}}.
#' Defaults to 11.
#' @importFrom grDevices gray
#' @importFrom stats qnorm
#' @import ggplot2
#' @export
plot_1D <- function(model, plot_type = 1, plot_ci = TRUE, alpha = 0.05, CIs_list, select,
                    xlab, ylab, ylim, main, effect_label, shift = 0, base_size = 11) {
  if(class(model)[1] != "pffr")
    stop("Please specify a fitted pffr model as 'model'!")
  if(plot_type == 2 & missing(CIs_list))
    stop("Please specify 'CIs_list'! Necessary as plot_type = 2.")
  if(missing(select))
    stop("Please specify which model term should be plotted using the 'select' argument!")
  
  plotObject <- no.plot(model)
  if (length(select) == 1) {
    plotObject <- plotObject[[select]]
    plot_data <- data.frame("x" = plotObject$x,
                            "fit" = shift + plotObject$fit)
  } else if (length(select) == 2) {
    plot_data <- data.frame("x" = plotObject[[select[1]]]$x,
                            "fit" = shift + plotObject[[select[2]]]$fit - plotObject[[select[1]]]$fit)
    plotObject <- plotObject[[select[1]]]
    if (missing(main))
      main <- paste0("Difference between smooths ",select[2]," and ",select[1])
    if (plot_type != 2)
      plot_ci <- FALSE
  }
  if (plot_ci) {
    if (plot_type == 1 & length(select) == 1) {
      plot_data$ci_lower <- shift + plotObject$fit - qnorm(1-alpha/2)*plotObject$se
      plot_data$ci_upper <- shift + plotObject$fit + qnorm(1-alpha/2)*plotObject$se
      if(missing(main)) {
        main <- ifelse(missing(effect_label), plotObject$ylab, effect_label)
        if (plot_ci)
          main <- paste0(main, "\nwith ",100-100*alpha,"% Marra & Wood (2012) confidence intervals")
      }
    } else if (plot_type == 2) {
      if (length(select) == 1) {
        ci <- CIs_list[[select]]
      } else if (length(select) == 2)
        ci <- CIs_list[[1]]
      plot_data$ci_lower <- ci$ci_lower
      plot_data$ci_upper <- ci$ci_upper
      if(missing(main)) {
        effect_label <- ifelse(missing(effect_label), plotObject$ylab, effect_label)
        main <- paste0(effect_label, "\nwith ", ci$ci_type, " ", 100-100*ci$alpha, "% confidence intervals")
      }
    }
  }
  if (!plot_ci)
    main <- ifelse(missing(main), NULL, main)
  
  xlab <- ifelse(missing(xlab), plotObject$xlab, xlab)
  ylab <- ifelse(missing(ylab), ifelse(length(select) == 1, plotObject$ylab, "Difference"), ylab)
  
  if (plot_ci)
    poly_dat <- data.frame("x" = c(plot_data$x, rev(plot_data$x)),
                           "y" = c(plot_data$ci_lower, rev(plot_data$ci_upper)))
  gg <- ggplot(plot_data, aes_string(x = "x", y = "fit"))
  if(!(plot_type == 1 & !plot_ci)) {
    gg <- gg + geom_polygon(data = poly_dat, aes_string(x = "x", y = "y"), fill = gray(0.75)) +
      geom_hline(yintercept = 0, lty = 2, col = "darkgray")
  }
  gg <- gg + geom_line() +
    theme_bw(base_size = base_size) +
    xlab(xlab) + ylab(ylab) + ggtitle(main) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
  if (!missing(ylim))
    gg <- gg + scale_y_continuous(limits = ylim)
  gg
}





#' 3D perspective plot for 2D smooth effects for \code{\link[refund]{pffr}} models
#'
#' This function produces 3D plots for 2D smooth effects of function-on-scalar 
#' models fitted with \code{\link[refund]{pffr}}.
#' This function mainly does the same as \code{mgcv::plot.gam(..., pers = TRUE)}.
#' However, \code{\link[mgcv]{plot.gam}} only uses the \code{shade} argument
#' for 1D smooths and doesn't pass it to the \code{\link{persp}} function. Such
#' arguments not passed to \code{persp} can be used with this function.
#'
#' @param model Function-on-scalar model fitted with \code{\link[refund]{pffr}}
#' @param select Index of the smooth effect to be plotted
#' @param theta,phi Some specific arguments for \code{\link{persp}} with
#' changed default values
#' @param ... Further arguments passed to \code{persp}
#' @importFrom graphics persp
#' @export
plot_2D <- function(model, select, theta = 30, phi = 30, ...) {
  p <- no.plot(model)
  P <- p[[select]]
  n <- length(P$x)
  persp(P$x, P$y, matrix(P$fit, n, n), theta = theta, phi = phi, ...)
}



#' Heatmat plot for 2D smooth effects for \code{\link[refund]{pffr}} models
#'
#' Plots heatmap plots for 2D smooth effects for a function-on-scalar model
#' fitted with \code{\link[refund]{pffr}}. Additionally, Marra \& Wood (2012) or
#' bootstrapped-based confidence intervals (CIs) can be plotted, specifying the
#' \code{plot_type} argument.
#'
#' @param plot_ci_type Only used if \code{plot_ci = TRUE}. One of
#' \code{c('all','lower','upper')}. Determines if the effect and the CI
#' boundaries (\code{'all'}), only the lower or only the upper boundary
#' (\code{'lower'/'upper'}) should be plotted
#' @param lowCol,highCol Colors defining the color gradient of the plots
#' @param legend_limits Optional vector specifying the limits of the color
#' scale, e.g. \code{c(-2,2)}
#' @param legend_position Either "bottom", "right" or "none"
#' @param sub Optional subtitle of the plot
#' @param outer_fontsize Fontsize of the main and sub which are outer of the
#' single ggplots. Defaults to 15.
#' @param newpage See \code{\link[FoSIntro]{grid_arrange_shared_legend}} function.
#' Should be set to FALSE when exporting the plot, e.g. using \code{\link{pdf}}.
#' @param ... Further arguments passed to \code{\link[ggplot2]{theme}}
#' @inheritParams plot_1D
#' @return The function silently returns the ggplot object
#' @import ggplot2
#' @export
plot_2Dheatmap <- function(model, plot_type = 1, plot_ci = TRUE, plot_ci_type = "all",
                           alpha = 0.05, CIs_list, select,
                           lowCol = "blue", highCol = "red", legend_limits, legend_position = "right",
                           xlab, ylab, main, sub = NULL, effect_label, base_size = 11,
                           outer_fontsize = 15, newpage = TRUE, ...) {
  if(class(model)[1] != "pffr")
    stop("Please specify a fitted pffr model as 'model'!")
  if(plot_type == 2 & missing(CIs_list))
    stop("Please specify 'CIs_list'! Necessary as plot_type = 2.")
  if(missing(select))
    stop("Please specify which model term should be plotted using the 'select' argument!")
  
  plotObject <- no.plot(model)
  plotObject <- plotObject[[select]]
  plot_data <- expand.grid(plotObject$x,
                           plotObject$y)
  plot_data$fit <- as.vector(plotObject$fit)
  if (plot_type == 1) {
    plot_data$ci_lower <- plot_data$fit - qnorm(1-alpha/2)*as.vector(plotObject$se)
    plot_data$ci_upper <- plot_data$fit + qnorm(1-alpha/2)*as.vector(plotObject$se)
    if(missing(main)) {
      main <- ifelse(missing(effect_label), plotObject$main, effect_label)
      if (plot_ci)
        main <- paste0(main, "\nwith ",100-100*alpha,"% Marra & Wood (2012) confidence intervals")
    }
  } else if (plot_type == 2) {
    ci <- CIs_list[[select]]
    plot_data$ci_lower <- ci$ci_lower
    plot_data$ci_upper <- ci$ci_upper
    if(missing(main)) {
      effect_label <- ifelse(missing(effect_label), plotObject$main, effect_label)
      main <- paste0(effect_label, "\nwith ", ci$ci_type, " ", 100-100*ci$alpha, "% confidence intervals")
    }
  }
  
  xlab <- ifelse(missing(xlab), plotObject$xlab, xlab)
  ylab <- ifelse(missing(ylab), plotObject$ylab, ylab)
  
  if (missing(legend_limits)) {
    if (!plot_ci) {
      legend_limits <- range(plot_data$fit, na.rm = T)
    } else
      legend_limits <- c(min(plot_data$ci_lower, na.rm = T), max(plot_data$ci_upper, na.rm = T))
  }
  gg_fit <- ggplot(plot_data, aes_string(x="Var1", y="Var2", fill = "fit")) + geom_tile() +
    xlab(xlab) + ylab("") + ggtitle("Point estimate") +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    scale_fill_gradient2(limits = legend_limits, low = lowCol, mid = "white", high = highCol, name = "estimate") +
    geom_contour(aes_string(x="Var1", y="Var2", z = "fit"), colour = "darkgrey") +
    theme_bw(base_size) + theme(plot.title = element_text(hjust = 0.5), legend.position = legend_position, ...)
  if(!plot_ci)
    gg_fit <- gg_fit + ylab(ylab)
  if(!(plot_type == 1 & !plot_ci)) {
    gg_lower <- ggplot(plot_data, aes_string(x="Var1", y="Var2", fill = "ci_lower")) + geom_tile() +
      xlab(xlab) + ylab(ylab) + ggtitle("Lower CI boundary") +
      scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
      scale_fill_gradient2(limits = legend_limits, low = lowCol, mid = "white", high = highCol, name = "estimate") +
      geom_contour(aes_string(x="Var1", y="Var2", z = "ci_lower"), colour = "darkgrey") +
      theme_bw(base_size) + theme(plot.title = element_text(hjust = 0.5), legend.position = legend_position, ...)
    ylab_upper <- ifelse(plot_ci_type == "upper", ylab, "")
    gg_upper <- ggplot(plot_data, aes_string(x="Var1", y="Var2", fill = "ci_upper")) + geom_tile() +
      xlab(xlab) + ylab(ylab_upper) + ggtitle("Upper CI boundary") +
      scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
      scale_fill_gradient2(limits = legend_limits, low = lowCol, mid = "white", high = highCol) +
      geom_contour(aes_string(x="Var1", y="Var2", z = "ci_upper"), colour = "darkgrey") +
      theme_bw(base_size) + theme(plot.title = element_text(hjust = 0.5), legend.position = legend_position, ...)
    
    if (plot_ci_type == "lower")
      (gg <- gg_lower)
    else if (plot_ci_type == "upper")
      (gg <- gg_upper)
    else if (plot_ci_type == "all")
      (gg <- grid_arrange_shared_legend(gg_lower, gg_fit, gg_upper,
                                        ncol = 3, legend_position = legend_position,
                                        main = main, sub = sub, outer_fontsize = outer_fontsize, newpage = newpage))
  } else { # plot effect without CIs
    xlab <- ifelse(is.null(xlab), sub, xlab)
    gg <- gg_fit + ggtitle(main) + xlab(xlab)
  }
  if (plot_ci == FALSE | plot_ci_type != "all")
    print(gg)
  invisible(gg) # invisibly return plot object
}




#' Plot model predictions for \code{\link[refund]{pffr}} models
#'
#' Plots the model predictions for a function-on-scalar model fitted with
#' \code{\link[refund]{pffr}}.
#'
#' @param model Function-on-scalar model fitted with \code{\link[refund]{pffr}}
#' @param newdata See \code{\link[refund]{predict.pffr}}
#' @param yvar name of the response variable in \code{model}
#' @param ci_type One of \code{"none"} (default; no confidence bands), \code{"ci"}
#' (pointwise confidence intervals for the predicted mean) or \code{"pi_norm"}/\code{"pi_gamma"}
#' (pointwise prediction intervals), specifying which type of uncertainty bands should be plotted.
#' Prediction intervals are either calculated summing the predicted mean variance and
#' the residual variance (\code{"pi_norm"}) or using 'posterior simulation' of the parameter vector
#' (\code{"pi_gamma"}) as outlined here: https://stat.ethz.ch/pipermail/r-help/2011-April/275632.html.
#' Note that the \code{"pi_norm"} PIs are not generally applicable for generalized models.
#' The posterior simulation approach is more generally applicable, but is currently only implemented
#' for Gamma regression.
#' @param ci_alpha Alpha to be used for confidence or prediction intervals.
#' Defaults to 0.05, i.e. 95\% intervals.
#' @param labels Labels vector with length equal to \code{nrow(newdata)}.
#' Used for the legend.
#' @param legend_title Title of legend
#' @param col_palette Color palette to use in plot, see
#' \code{\link[ggplot2:scale_colour_brewer]{scale_color_brewer}}. Only used if \code{nrow(newdata)>1}
#' and if \code{is.null(col_vector)}.
#' @param col_vector Vector of colors to use
#' @param rev_cols If TRUE, the color palette is turned around.
#' @param base_size Base size of plot elements, see \code{\link[ggplot2:ggtheme]{theme_bw}}.
#' Defaults to 11.
#' @param legend_symbol_size Size of symbols in legend. Default to 2
#' @param type \code{type} argument for \code{\link[refund]{predict.pffr}}
#' @param log10 If TRUE the predictions are shown on a log10 scale
#' @param plot.margin Margin of the ggplot
#' @param hline Vector of \code{\link[ggplot2:geom_abline]{geom_hline}}
#' @param ybreaks Numeric vector of breaks for the y axis
#' @param ylabels Character vector of labels for the y axis. Only used if
#' \code{ybreaks} is specified
#' @param xlab,ylab,main,ylim,lwd Further arguments for the plot function
#' @param ... Additional arguments passed to \code{\link{pi_gamma}}
#' @importFrom stats var
#' @importFrom tidyr gather_
#' @import dplyr
#' @import magrittr
#' @import ggplot2
#' @import refund
#' @export
plot_predictions <- function(model, newdata, yvar, ci_type = "none", ci_alpha = 0.05, labels = NULL, legend_title = "group",
                             col_palette = "Set1", col_vector = NULL, rev_cols = FALSE,
                             base_size = 11, legend_symbol_size = 2, type = "response", log10 = FALSE,
                             plot.margin = unit(rep(5.5,4),"points"), hline, ybreaks = waiver(), ylabels = waiver(),
                             xlab = "yindex", ylab = "prediction", main, ylim = NULL, lwd = 1, ...) {
  if (missing(yvar))
    stop("Please specify the response variable yvar")
  newdata[,yvar] <- NULL # predict.pffr only works if no yvar is in the data.frame
  if (ci_type %in% c("none","pi_gamma"))
    fit <- as.data.frame(refund:::predict.pffr(model, newdata, type = type))
  else if (ci_type %in% c("ci","pi_norm")) {
    fit <- refund:::predict.pffr(model, newdata, type = type, se.fit = TRUE)
    se_fit <- as.data.frame(fit$se.fit)
    colnames(se_fit) <- model$pffr$yind
    fit <- as.data.frame(fit$fit)
  }
  colnames(fit) <- model$pffr$yind
  if (!is.null(labels)) {
    fit$label <- labels
    fit$label <- factor(fit$label, levels = labels)
  } else
    fit$label <- paste("group", 1:nrow(fit))
  plot_data <- fit %>%
    tidyr::gather_(key = "yindex", value = "prediction", gather_cols = colnames(fit)[!(colnames(fit) %in% c("yindex","prediction","label"))])
  if (ci_type %in% c("ci","pi_norm")) {
    se_fit_data <- se_fit %>%
      gather(key = "yindex", value = "se")
    plot_data$se <- se_fit_data$se
  }
  plot_data$yindex <- as.numeric(plot_data$yindex)
  
  if (missing(main)) {
    alpha_label <- as.character(round(100 - 100*ci_alpha,1))
    main <- paste0("predictions", ifelse(ci_type == "none", "",
                                         ifelse(ci_type == "ci", paste0("\nwith ",alpha_label,"% pointwise confidence intervals\nfor the predicted mean"),
                                                paste0("\nwith ",alpha_label,"% pointwise prediction intervals"))))
  }
  if (class(ylabels) == "waiver" & class(ybreaks) != "waiver") # if only ybreaks are specified
    ylabels <- as.character(ybreaks)
  if (is.null(ylim) & class(ybreaks) != "waiver")
    ylim <- range(ybreaks) # if ybreaks are specified, but ylim not
  
  # create dataset with confidence/prediction intervals
  if (ci_type != "none") {
    if (ci_type == "pi_norm") {
      sigma2 <- var(as.vector(refund:::residuals.pffr(model, type = "response")), na.rm = T) # error variance
      plot_data %<>% mutate(se = sqrt(plot_data$se^2 + sigma2))
    }
    if (ci_type %in% c("ci","pi_norm")) {
      plot_data %<>% mutate(interval_lower = plot_data$prediction - qnorm(1-(ci_alpha/2))*plot_data$se,
                            interval_upper = plot_data$prediction + qnorm(1-(ci_alpha/2))*plot_data$se)
    }
    if (ci_type == "pi_gamma") {
      pi_dat <- pi_gamma(model, newdata, ci_alpha, ...)
      plot_data %<>% mutate(interval_lower = pi_dat$pi_lower, interval_upper = pi_dat$pi_upper)
    }
    # if y is plotted on the log10 scale and the intervals are <0, then set those values to 0.00001
    if (log10 & any(plot_data$interval_lower < 0)) {
      warning("Confidence/Prediction interval lower boundaries < 0 are set to 0.00001 in order to plot on log10 scale!")
      plot_data$interval_lower[plot_data$interval_lower < 0] <- 0.00001
    }
    polydat <- data.frame("yindex" = c(plot_data$yindex, rev(plot_data$yindex)),
                          "label" = factor(c(as.character(plot_data$label), rev(as.character(plot_data$label)))),
                          "interval_border" = c(plot_data$interval_lower, rev(plot_data$interval_upper)))
  }
  
  gg <- ggplot(plot_data, aes_string(x="yindex", y="prediction", color="label"))
  if (ci_type != "none") {
    gg <- gg + geom_polygon(data = polydat, aes_string(x = "yindex", y = "interval_border", fill = "label", col = NULL), show.legend = FALSE, alpha = 0.3)
    if (nrow(newdata) == 1)
      gg <- gg + scale_fill_manual(values = "black")
    else if (is.null(col_vector))
      gg <- gg + scale_fill_brewer(name = legend_title, palette = col_palette, direction = ifelse(rev_cols, -1, 1))
    else
      gg <- gg + scale_fill_manual(name = legend_title, values = col_vector)
  }
  if (!missing(hline))
    gg <- gg + geom_hline(yintercept = hline, lty = 2, col = "grey", lwd = 1.1)
  if (log10)
    gg <- gg + scale_y_log10(breaks = ybreaks, labels = ylabels, limits = ylim, name = ylab)
  else
    gg <- gg + scale_y_continuous(breaks = ybreaks, labels = ylabels, limits = ylim, name = ylab)
  gg <- gg +
    geom_line(lwd = lwd) +
    xlab(xlab) + ggtitle(main) +
    guides(colour = guide_legend(override.aes = list(size = legend_symbol_size))) +
    theme_bw(base_size) +
    theme(legend.position = ifelse(nrow(newdata) == 1, "none", "bottom"),
          plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank(),
          plot.margin = plot.margin)
  if (nrow(newdata) == 1)
    gg <- gg + scale_color_manual(values = "black")
  else if (is.null(col_vector))
    gg <- gg + scale_color_brewer(name = legend_title, palette = col_palette, direction = ifelse(rev_cols, -1, 1))
  else
    gg <- gg + scale_color_manual(name = legend_title, values = col_vector)
  gg
}




#' Plot predictions vs real observations for \code{\link[refund]{pffr}} models
#'
#' Plots the predictions o a function-on-scalar model fitted with
#' \code{\link[refund]{pffr}} agains the real observations.
#'
#' @param model Function-on-scalar model fitted with \code{\link[refund]{pffr}}
#' @param data Dataset with one row for which the predictions are to be plotted
#' @param yvar name of the response variable in \code{model}
#' @param yvar_label Optional label of the response variable used for plot annotation
#' @param type \code{"response"} (default) or \code{"link"}.
#' See \code{\link[refund]{predict.pffr}}.
#' @param log10 If TRUE the predictions are shown on a log10 scale
#' @param pi_type One of \code{"none"} (default; no uncertainty bands), \code{"pi_norm"}
#' or \code{"pi_gamma"} (pointwise prediction intervals), specifying which type of uncertainty bands should be plotted.
#' Prediction intervals are either calculated summing the predicted mean variance and
#' the residual variance (\code{"pi_norm"}) or using 'posterior simulation' of the parameter vector
#' (\code{"pi_gamma"}) as outlined here: https://stat.ethz.ch/pipermail/r-help/2011-April/275632.html.
#' Note that the \code{"pi_norm"} PIs are not generally applicable for generalized models.
#' The posterior simulation approach is more generally applicable, but is currently only implemented
#' for Gamma regression.
#' @param pi_alpha If \code{pi_type != "none"}, \code{1-pi_alpha}\% prediction
#' intervals are plotted. Defaults to 0.05.
#' @param base_size Base size of plot elements, see \code{\link[ggplot2:ggtheme]{theme_bw}}.
#' Defaults to 11.
#' @param xlab,ylab,main Optional arguments for plot annotation
#' @param ybreaks,ylabels Optional specification of the tick position and the
#' corresponding labels (character vector) for the y axis
#' @param ylim Optional y axis limits of plot
#' @param lwd Optional line width. Defaults to 2
#' @param legend_symbol_size Optional size of symbols in legend.
#' @param seed Additional arguments passed to \code{\link{pi_gamma}}
#' @importFrom tidyr gather
#' @import refund
#' @import ggplot2
#' @export
plot_predVSobs <- function(model, data, yvar, yvar_label, type = "response", log10 = FALSE,
                           pi_type = "none", pi_alpha = 0.05, base_size = 11,
                           xlab = "yindex", ylab, main = "Prediction vs observation",
                           ybreaks, ylabels, ylim = NULL, lwd = 1, legend_symbol_size = 2, ...) {
  if (missing(yvar))
    stop("Please specify the response variable yvar")
  if (missing(yvar_label))
    yvar_label <- yvar
  if (pi_type != "none" & type != "response")
    stop("Note: Prediction intervals can currently only be plotted for type == 'response'!")
  if (missing(ylab))
    ylab <- ifelse(type == "link", paste0("logarithmized\n",yvar_label),
                   ifelse(type == "response" & !log10, yvar_label,
                          ifelse(type == "response" & log10, paste0(yvar_label,"\non log10-scale"))))
  obs <- data[,yvar]
  data[,yvar] <- NULL
  if (pi_type != "none") {
    p <- refund:::predict.pffr(model, newdata = data, type = type, se.fit = TRUE)
    if (pi_type == "pi_norm") {
      sigma2 <- var(as.vector(refund:::residuals.pffr(model, type = "response"))) # error variance
      se.fit <- sqrt(p$se.fit^2 + sigma2)
      p <- p$fit
      pi_lower <- p - qnorm(1 - pi_alpha/2) * se.fit
      pi_upper <- p + qnorm(1 - pi_alpha/2) * se.fit
    } else if (pi_type == "pi_gamma") {
      pi_dat <- pi_gamma(model, data, pi_alpha, ...)
      pi_lower <- pi_dat$pi_lower
      pi_upper <- pi_dat$pi_upper
      p <- p$fit
    }
    if (log10) {
      if (any(pi_lower <= 0)) {
        warning("Prediction intervals are suppressed as they include negative values on the non-log scale!")
        pi_type <- "none"
      } else {
        pi_lower <- log10(pi_lower)
        pi_upper <- log10(pi_upper)
      }
    }
  } else
    p <- refund:::predict.pffr(model, newdata = data, type = type)
  colnames(p) <- colnames(obs)
  if (type == "link")
    obs <- log(obs)
  if (missing(ylabels) & !missing(ybreaks)) ylabels <- as.character(ybreaks)
  
  plot_data <- rbind(obs, p)
  plot_data$type <- c("observation","prediction")
  plot_data <- tidyr::gather(plot_data, key = "yindex", value = "y", -type)
  plot_data$yindex <- as.numeric(substr(plot_data$yindex, 3,6))
  
  gg <- ggplot(plot_data, aes_string(x="yindex", y="y", color="type"))
  if (pi_type != "none") {
    yindex <- unique(plot_data$yindex)
    poly_data <- data.frame("yindex" = c(yindex, rev(yindex)),
                            "pi_border" = c(pi_lower, rev(pi_upper)))
    gg <- gg + geom_polygon(data = poly_data, aes(x = yindex, y = pi_border), 
                            fill = "dodgerblue3", color = "dodgerblue3", alpha = 0.5)
  }
  gg <- gg + geom_line(lwd = lwd) +
    xlab(xlab) + ylab(ylab) + ggtitle(main) +
    scale_color_manual(name = NULL, values = c("black","dodgerblue3")) +
    guides(colour = guide_legend(override.aes = list(size = legend_symbol_size))) +
    theme_bw(base_size = base_size) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "bottom")
  if (log10) {
    if (!missing(ybreaks))
      gg <- gg + scale_y_log10(breaks = ybreaks, labels = ylabels, limits = ylim, name = ylab)
    else
      gg <- gg + scale_y_log10(limits = ylim, name = ylab)
  } else {
    if (!missing(ybreaks))
      gg <- gg + scale_y_continuous(breaks = ybreaks, labels = ylabels, limits = ylim, name = ylab)
    else
      gg <- gg + scale_y_continuous(limits = ylim, name = ylab)
  }
  gg
}
