
#' Perform bootstrap for smooth effects in \code{\link[refund]{pffr}} models
#'
#' With this function you can perform a parametric or nonparametric
#' bootstrap for a function-on-scalar model fitted with \code{\link[refund]{pffr}}.
#' Based on the results, confidence intervals (CIs) for smooth effects can be
#' calculated using \code{link{calc_bootstrapCIs}}.
#'
#' @param B Number of bootstrap samples, defaults to 1000
#' @param type One of \code{"parametric"} or \code{"nonparametric"}, specifying
#' the type of bootstrap to be used. If \code{type = "parametric"} you also have
#' to specify the \code{model} argument.
#' @param formula See formula argument of \code{\link[refund]{pffr}}
#' @param data See data argument of \code{pffr}
#' @param cores Number of cores to use for parallel processing. Possible for
#' both Linux-based systems and Windows.
#' @param model Function-on-scalar model fitted with \code{pffr} used for
#' parametric bootstrapping. Only used if \code{type=='parametric'}.
#' @param param_yvar Name of response variable used in \code{model} and present
#' in \code{data}. Only used if \code{type=='parametric'}.
#' @param param_simFun Specifies the simulation function to be used for
#' parametric bootstrapping. See \code{\link{simulate_pffr}}.
#' @param param_yMinValue Minimum value to which the bootstrapped response values
#' are set using a parametric bootstrap. See \code{\link{simulate_pffr}}.
#' @param log_file Optional filename of a log file where progress of parallel
#' processing should printed to. Defaults to console output (only on Linux-based
#' systems). Only used if \code{cores>1}.
#' @param ... Further arguments passed to \code{pffr}
#' @importFrom refund pffr
#' @import parallel
#' @export
bootstrap_pffr <- function(type = "parametric", B = 1000, formula, data, model = NULL,
                           param_yvar = NULL, cores = 1, param_simFun = "simulate",
                           param_yMinValue = NULL, log_file = NULL, ...) {
  if (!(type %in% c("parametric","nonparametric")))
    stop("Please specify either 'parametric' or 'nonparametric' for the 'type' argument!")
  if (type == "parametric" & class(model)[1] != "pffr")
    stop("Please specifiy the 'model' argument as a fitted pffr model!")
  if (type == "parametric" & (is.null(param_yvar) || !(param_yvar %in% colnames(data))))
    stop("Please specify the 'param_yvar' argument correctly as a column name of data and the response variable used in model!")

  if (!is.null(log_file) && file.exists(log_file)) # delete log_file for overwriting it
   unlink(log_file, recursive = TRUE)
  bs_results <- vector(mode = "list", length = B)
  singleIteration <- function(b, B, type, formula, data, model, param_yvar, param_simFun, param_yMinValue, log_file = NULL) { # b and B are only used to print progress of the whole bootstrap
    if(!is.null(log_file))
      sink(log_file, append = TRUE)
    print(paste0("Performing bootstrap iteration ",b,"/",B,"..."))
    if (type == "parametric") {
      data_b <- data
      data_b[,param_yvar] <- FoSIntro::simulate_pffr(model, sim_fun = param_simFun, y_minValue = param_yMinValue)
    } else if (type == "nonparametric")
      data_b <- data[sample(seq_len(nrow(data)), size = nrow(data), replace = TRUE),]
    model_b <- pffr(formula, data = data_b, ...)
    if(!is.null(log_file))
      sink()
    return(no.plot(model_b)) # save returned value of plot.gam while suppressing plot creation
  }

  # perform bootstrap
  if (cores == 1) { # no parallel call
    bs_results <- lapply(seq_len(B), function(b) { singleIteration(b,B,type,formula,data,model,param_yvar,param_simFun,param_yMinValue) })
  } else if (Sys.info()["sysname"] != "Windows") { # parallel call for Linux-based systems
    bs_results <- mclapply(seq_len(B), function(b) { singleIteration(b,B,type,formula,data,model,param_yvar,param_simFun,param_yMinValue,log_file) }, mc.cores = cores)
  } else { # parallel call for Windows
    local_cluster <- makePSOCKcluster(rep("localhost", cores)) # create cluster
    # Export functions and packages to the cluster
    clusterExport(cl = local_cluster, c("singleIteration","no.plot","B","type","formula","data","model","param_yvar","param_simFun","param_yMinValue","log_file",
                                        "pffr","yindex"), envir=environment())
    clusterEvalQ(cl = local_cluster, c(library(parallel), library(refund),library(mgcv)))

    bs_results <- parLapply(cl = local_cluster, X = seq_len(B),
                         fun = function(b) { singleIteration(b,B,type,formula,data,model,param_yvar,param_simFun,param_yMinValue,log_file) })
    stopCluster(cl = local_cluster) # close cluster
  }
  if(!is.null(log_file))
    sink(log_file, append = TRUE)
  print("Finished!")
  if(!is.null(log_file))
    sink()
  return(bs_results)
}


#' Calculate bootstrap-based confidence intervals for smooth effects in \code{\link[refund]{pffr}} models
#'
#' Wrapper for \code{\link{calc_singleBootstrapCI}} that calculates CIs for all
#' smooth terms in \code{bs_results}.
#'
#' @inheritParams calc_singleBootstrapCI
#' @export
calc_bootstrapCIs <- function(bs_results, alpha = 0.05) {
  terms <- length(bs_results[[1]])
  CIs_list <- lapply(seq_len(terms), function(i) calc_singleBootstrapCI(bs_results, i, alpha))
  return(CIs_list)
}

#' Calculate a pointwise confidence interval band for a single smooth term in a \code{\link[refund]{pffr}} model
#'
#' Function to calculate a pointwise confidence interval (CI) band for a single
#' smooth term, based on the output of \code{\link{bootstrap_pffr}}.
#' This function is normally called more comfortably by using \code{\link{calc_bootstrapCIs}}.
#' However, compared to the latter function, it also offers the additional
#' possibility to compute CIs for differences of two smooths.
#' NOTE: Only 1D and 2D smooths are currently handled correctly!
#'
#' @param bs_results Output of function \code{\link{bootstrap_pffr}}
#' @param select Index of smooth term to be looked at. If a vector of length two
#' is passed, then the CI will be computed for the difference of the second smooth
#' minus the first smooth.
#' @param alpha \code{(1-alpha)} CIs are calculated. The default 0.05 leads to 95\% CIs
#' @importFrom stats quantile
#' @export
calc_singleBootstrapCI <- function(bs_results, select, alpha = 0.05) {
  x <- lapply(bs_results, function(x) x[[select[1]]])
  if (length(select) == 2)
    x2 <- lapply(bs_results, function(x) x[[select[2]]])
  B <- length(x)
  fit_length <- length(x[[1]]$fit)
  q <- c(alpha/2, 1-alpha/2)
  if(!("y" %in% names(x[[1]]))) # if 1D smooth effect(s)
    y <- NULL
  else # if 2D smooth effect(s)
    y <- x[[1]]$y
  fits <- lapply(seq_len(fit_length), function(i) {
    if (length(select) == 1) {
      return(sapply(seq_len(B), function(b) x[[b]]$fit[i]))
    } else if (length(select) == 2) {
      return(sapply(seq_len(B), function(b) x2[[b]]$fit[i] - x[[b]]$fit[i]))
    }
  })
  # calculate pointwise CIs
  fitCI_lower <- sapply(fits, function(f) unname(stats::quantile(f, probs = q[1])))
  fitCI_upper <- sapply(fits, function(f) unname(stats::quantile(f, probs = q[2])))

  ci_list <- list("x" = x[[1]]$x,
                  "y" = y,
                  "ci_lower" = fitCI_lower,
                  "ci_upper" = fitCI_upper,
                  "ci_type" = "pointwise",
                  "alpha" = alpha)
  return(ci_list)
}

#' Simulate new observations from a \code{pffr} model
#'
#' This function can be used to draw a parametric bootstrap sample from a
#' function-on-scalar model fitted with \code{\link[refund]{pffr}}. Per default,
#' new observations are sampled using \code{\link[stats]{simulate}}. If the
#' distribution family used is not supported by \code{simulate}, one can also
#' define and pass an own function as \code{sim_fun}. One example for this is
#' \code{\link[FoSIntro]{simulate_gamma}}, which supports simulation from a
#' Gamma distribution and can be called using \code{sim_fun = "gamma"}.
#'
#' @param model Function-on-scalar model fitted with \code{\link[refund]{pffr}}
#' @param sim_fun Function that should be used to draw the parametric bootstrap
#' sample based on \code{model}. Default value is \code{'simulate'} which uses
#' the \code{\link[stats]{simulate}} function. For distribution families not
#' supported in \code{\link[stats]{simulate}} one can either use \code{sim_fun = 'gamma'}
#' to perform bootstrapping for Gamma models or one can pass the name of its own function
#' to perform the bootstrap. The function should thereby take a \code{pffr} model
#' as single argument and return a bootstrapped matrix of y values with one row
#' per observation and one column per functional index (e.g. time point), equivalent to
#' \code{\link[FoSIntro]{simulate_gamma}}.
#' @param y_minValue Optional argument. If specified, possible bootstrapped response
#' values smaller than \code{y_minValue} are set to \code{y_minValue}. This can e.g.
#' be used for Gamma models to prevent estimation problems of models based on the
#' bootstrapped data, caused by response values being too close to zero.
#' @return Matrix of reponse values with one row per observation and one column
#' per functional index (e.g. time point).
#' @export
simulate_pffr <- function(model, sim_fun = "simulate", y_minValue = NULL) {
  if (sim_fun == "simulate") {
    pred <- as.data.frame(matrix(stats::simulate(model), byrow = FALSE, ncol = length(model$pffr$yind)))
  } else if (sim_fun == "gamma") {
    pred <- simulate_gamma(model)
    if (is.null(y_minValue))
      message("Note: You didn't specify a <y_minValue> and it could happen that simulated y values too close to zero cause estimation problems in following Gamma models...")
  } else {
    sim_fun <- get(sim_fun)
    pred <- sim_fun(model)
  }
  colnames(pred) <- paste0("y_", model$pffr$yind)
  if (!is.null(y_minValue))
    pred <- replace(pred, pred <= y_minValue, y_minValue)
  pred
}

#' Function to draw one parametric bootstrap sample for a Gamma model
#'
#' This function draws one parametric bootstrap sample from a \code{\link[refund]{pffr}}
#' model using a Gamma response distribution.
#'
#' @param model Function-on-scalar model fitted with \code{\link[refund]{pffr}}
#' @importFrom stats rgamma
#' @import refund
simulate_gamma <- function(model) {
  pred_mean <- refund:::predict.pffr(model, type = "response")
  scale <- model$sig2
  pred <- apply(pred_mean, c(1,2), function(x) rgamma(1, shape = x / scale, scale = scale))
}



#' Plot bootstrapped confidence intervals
#'
#' Function to plot bootstrapped confidence intervals for function-on-scalar
#' models fitted with \code{\link[refund]{pffr}}. The bootstrap confidence intervals
#' have to be precalculated beforehand using \code{\link{calc_bootstrapCIs}}.
#'
#' @param model Function-on-scalar model fitted with \code{\link[refund]{pffr}},
#' which was fit on the raw, unbootstrapped data to extract point estimate.
#' @param CIs_list Output of function \code{\link{calc_bootstrapCIs}}
#' @param select Index of smooth term to be plotted. Default plots all effects
#' @param xlab,ylab,main Optional parameters to overwrite the default
#' plot parameters
#' @param effect_label Used in main as '\code{effect_label} with global 95\% CI'.
#' Only used if \code{main = NULL}.
#' @param ... Further arguments passed to \code{\link{plot_2Dheatmap}}
#' @return Silently returns list with plot objects
#' @export
plot_bootstrapCIs <- function(model, CIs_list, select = NULL,
                              xlab, ylab, main, effect_label, ...) {
  if(is.null(select))
    select <- 1:length(CIs_list)
  plotObjects <- vector(mode = "list", length = length(select))
  for(i in select) { # plot each effect with CIs
    ci <- CIs_list[[i]]
    if (model$smooth[[i]]$dim == 1) { # 1D smooth effect
      plotObjects[[match(i,select)]] <- plot_1D(model, plot_type = 2, CIs_list = CIs_list,
                                                select = i, xlab = xlab, ylab = ylab, main = main, effect_label = effect_label)
      print(plotObjects[[match(i,select)]])
    } else if (model$smooth[[i]]$dim == 2) { # 2D smooth effect
      plotObjects[[match(i,select)]] <- plot_2Dheatmap(model, plot_type = 2, CIs_list = CIs_list,
                                                       select = i, xlab = xlab, ylab = ylab, main = main, effect_label = effect_label, ...)
    } else
      stop("Only 1D and 2D smooth effects are currently handled correctly!")
    if(i < max(select))
      readline(prompt="Press [enter] for next plot: ")
  }
  invisible(plotObjects)
}
