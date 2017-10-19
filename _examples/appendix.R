
# load packages and helper functions --------------------------------------
# Install FoSIntro
devtools::install_github("bauer-alex/FoSIntro")

packages <- c("FoSIntro", "refund", "ggplot2", "grid", "gridExtra",
              "dplyr", "tidyr", "magrittr")
# install.packages(packages)
sapply(packages, require, character.only = TRUE)

# set ggplot2 theme
theme_set(theme_bw())

# data preparation --------------------------------------------------------
data <- FoSIntro::sample_data
yindex <- attr(data, "yindex")

# Manually create dummy variables for factor variables
# (currently necessary as pffr v.0.1-16 doesn't dummy code categorical variables correctly)
data <- cbind(data, with(data, model.matrix(~velocity))[,-1,drop=FALSE])

# Create an index variable that is needed for specifying smooth error terms in the model
# (which are equal to observational-specific functional random intercepts)
data$index <- factor(1:nrow(data))
# inclusion in model with term <s(index, bs = "re")>.
# Not currently feasible with pffr for datasets of this size.

# reshape data for some visualizations
data_wide <- data %>%
  select(-ground_velocity) %>%
  bind_cols(data$ground_velocity) %>%
  gather(key = "time", value = "ground_velocity", starts_with("y")) %>%
  mutate(time = as.numeric(substring(time, 3)))



# Data visualization ------------------------------------------------------
# Plot single observations
gg <- data_wide %>%
  filter(seismometer %in% c(62,2524,5936), simulation == 88) %>%
  ggplot(aes(x=time, y=ground_velocity)) +
  geom_line() +
  scale_y_log10() +
  ggtitle("Three single observations (for simulation 88)") +
  ylab("ground velocity\non log10-scale")
gg + aes(color=seismometer)
gg + facet_wrap(~ seismometer) +
  theme(legend.position = "none")

# Plot all observations at specific seismometers
# to see seismometer-specific variation
# Note: As the data only contains all simulations per seismometer where
#       15 seconds >0 where observed the seismometers have different n's
data_wide %>%
  filter(seismometer %in% c(62,2524,5936)) %>%
  ggplot(aes(x=time, y=ground_velocity, group=simulation)) +
  geom_line(alpha = 0.2) +
  scale_y_log10() +
  facet_wrap(~ seismometer) +
  ggtitle("All observations for specific seismometers") +
  ylab("ground velocity\non log10-scale")

# Plot all observations for one simulation
data_wide %>%
  filter(simulation == 88) %>%
  ggplot(aes(x=time, y=ground_velocity, group=seismometer)) +
  geom_line(alpha = 0.2) +
  scale_y_log10() +
  ggtitle("All observations for simulation 88") +
  ylab("ground velocity\non log10-scale")



# Some simple regression models -------------------------------------------
### 1) Intercept-only model
m1 <- pffr(ground_velocity ~ 1,
           family = Tweedie(p = 2, link = "log"), yind = yindex, data = data)
summary(m1)

### 2) Model with one factor covariate: hypocentral distance high/medium/low
m2_dat <- data
m2_dat$distance_cat <- cut(m2_dat$hypocentral_distance, breaks = c(17,25,40,103),
                         labels = c("large","medium","small"))
# manually create dummy variables for factor variables (necessary as of a bug in pffr v.0.1-16)
m2_dat <- cbind(m2_dat, with(m2_dat, model.matrix(~distance_cat))[,-1,drop=FALSE])
m2 <- pffr(ground_velocity ~ distance_catmedium + distance_catsmall,
           family = Tweedie(p = 2, link = "log"), yind = yindex, data = m2_dat)
summary(m2)


# The main regression model -----------------------------------------------
model <- pffr(ground_velocity ~ s(hypocentral_distance) + s(stat_friction_coef) +
                c(dyn_friction_coef) + c(s(slip_weakening)) +
                dir_background_stress + velocitysediment,
              family = Tweedie(p = 2, link = "log"), yind = yindex, data = data)



# Effect plots ------------------------------------------------------------
# Note: 1D and 2D smooth effects can both be plotted using the default
#       plot.pffr() function (based on plot.gam()). For better customization
#       however we use our own functions based on ggplot2.

### 1) Only-intercept model
# plot time-dependent intercept
plot(m1, rug = FALSE) # plot.pffr()
plot_1D(m1, select = 1, xlab = "time [s]")
# plot prediction of the model
# Note: newdata has to have some entry, even if the content isn't used for the intercept-only model
newdata <- data.frame("intercept" = 1)
ybreaks <- 10^c(-2,-1,0,1)
plot_predictions(m1, newdata, log10 = TRUE,
                 ylim = c(.01,10), ybreaks = ybreaks, hline = ybreaks, lwd = 2,
                 xlab = "time [s]", ylab = "prediction on log10-scale",
                 main = "Overall mean ground velocity\nbased on an intercept-only model")

### 2) Model with one factor covariate: hypocentral distance high/medium/low
plot(m2, pages = 1, rug = FALSE) # plot.pffr()
ylim <- c(-3,0.5)
gg1 <- plot_1D(m2, select = 1, xlab = "time [s]", ylim = ylim)
gg2 <- plot_1D(m2, select = 2, xlab = "time [s]", ylim = ylim)
gg3 <- plot_1D(m2, select = 3, xlab = "time [s]", ylim = ylim)
blank_plot <- grid::grid.rect(gp=gpar(col="white"))
grid.arrange(gg1, blank_plot, gg2, gg3, nrow = 2, ncol = 2)
# plot predictions for the three factor levels
newdata <- data.frame("distance_catmedium" = c(0,1,0),
                      "distance_catsmall" = c(0,0,1))
labels <- c("small","medium","large")
plot_predictions(m2, newdata, log10 = TRUE,
                 xlab = "time [s]", ylab = "predictions on log10-scale", lwd = 2,
                 legend_title = "hypocentral distance", labels = labels,
                 col_vector = c("#08519C","#3182BD","#9ECAE1"), rev_cols = TRUE)

### 3) Full model
ylim <- c(-2.25,3) # use same effect scale for all effects for easier interpretation
# intercept
plot_1D(model, select = 1, ylim = ylim)
# time-constant effect of slip weakening
plot_1D(model, select = 4, ylim = ylim)
# time-varying effect of the direction of the background stress
plot_1D(model, select = 5, xlab = "time [s]", ylim = ylim)
# time-varying effect of hypocentral distance
plot_2Dheatmap(model, select = 2, plot_ci = FALSE,
               xlab = "hypocentral distance", ylab = "time [s]", legend_limits = ylim)
plot_2D(model, select = 2,
        xlab = "hypocentral distance", ylab = "time [s]", zlab = "estimate")
# time-varying effect of the static coefficient of friction
plot_2Dheatmap(model, select = 3, plot_ci = FALSE,
               xlab = "static coefficient of friction", ylab = "time [s]", legend_limits = ylim)
plot_2D(model, select = 3,
        xlab = "static coefficient of friction", ylab = "time [s]", zlab = "estimate")



# Uncertainty quantification ----------------------------------------------
# Note: Currently only pointwise confidence/prediction intervals can be calculated with our code.
#       A method for computing simultaneous/global confidence bands
#       (and analogously intervalwise bands and global/intervalwise prediction bands)
#       is shown in Krivobokova et al. (2010).

### 1) Confidence intervals after Marra & Wood (2012)
# for 1D smooth effects
plot_1D(model, select = 4)
# for 2D smooth effects
plot_2Dheatmap(model, select = 2,
               xlab = "hypocentral distance", ylab = "time [s]")
plot_2Dheatmap(model, select = 3,
               xlab = "static coefficient of friction", ylab = "time [s]")

### 2) Bootstrap confidence intervals
# As an example, we only use B = 10 bootstrap iterations.
B <- 10 # bootstrap iterations
n_cores <- 1 # for parallel computation (possible both for Linux-based systems and Windows)
formula <- ground_velocity ~ s(hypocentral_distance) + s(stat_friction_coef) + c(dyn_friction_coef) +
  c(s(slip_weakening)) + dir_background_stress + velocitysediment
# Possibility 1: Parametric Bootstrap
bs_results <- bootstrap_pffr("parametric", B = B, cores = n_cores,
                             formula, data, model = model,
                             family = Tweedie(p = 2, link = "log"), yind = yindex, # arguments passed to pffr()
                             param_yvar = "ground_velocity", param_simFun = "gamma", param_yMinValue = 0.01,
                             log_file = "bootstrap.log") # useful for parallel computation

# Possibility 2: Nonparametric bootstrap
bs_results <- bootstrap_pffr("nonparametric", B = B, cores = n_cores,
                             formula, data, model = model,
                             family = Tweedie(p = 2, link = "log"), yind = yindex, # arguments passed to pffr()
                             log_file = "bootstrap.log") # useful for parallel computation
# Prepare bs_results for plotting
CIs_list <- calc_bootstrapCIs(bs_results, 0.05)
# Plots
plot_bootstrapCIs(model, CIs_list, select = 6, xlab = "time [s]", ylab = "estimate",
                  effect_label = "Time-varying effect of sediment velocity")
plot_bootstrapCIs(model, CIs_list, select = 2, ylab = "time [s]")


### 2) Confidence intervals for the predicted mean
newdata <- data[1,colnames(data) != "ground_velocity"]
plot_predictions(model, newdata, log10 = FALSE, ci_type = "ci", # ybreaks = ybreaks,
                 xlab = "time [s]", ylab = "prediction on original scale")

### 3) prediction intervals
plot_predictions(model, newdata, log10 = FALSE, ci_type = "pi",
                 xlab = "time [s]", ylab = "prediction on original scale")
# -> the pointwise prediction intervals in our application are much
#    wider as the pointwise confidence intervals for the predicted mean!
#    This makes sense: We can predict the MEAN expected ground velocity very
#    well with our model, but we cannot predict the ground velocities for
#    SINGLE earthquakes equally well as these are highly volatile.



# Hypothesis testing ------------------------------------------------------
# Note: As noted in the text we didn't report p-values as our quite high-dimensional
#       dataset led to all p-values being <0.0001.
#       However, the following code shows how the p-values for all shown
#       hypothesis tests can be extracted/computed.
# Save the model summary as we will extract some p-values from it
sm <- summary(model)

# Likelihood ratio test
# Example: Test if hypocentral_distance has a time-varying effect by comparing
#          the full model to 'model_2', which only contains a time-constant
#          effect of hypocentral distance.
# Note: as of a current bug in pffr v.0.1-16 we have to overwrite the class of a pffr-object
#       in order to make the anova.gam() function work properly
model_2 <- pffr(ground_velocity ~ c(s(hypocentral_distance)) + s(stat_friction_coef) + c(dyn_friction_coef) +
                  c(s(slip_weakening)) + dir_background_stress + velocitysediment,
                family = Tweedie(p = 2, link = "log"), yind = yindex, data = data)
mtest <- model; mtest_2 <- model_2
class(mtest) <- class(mtest)[-1]
class(mtest_2) <- class(mtest_2)[-1]
anova(mtest, mtest_2, test = "F")

### 1) Is the linear effect different from zero?
# 1a) for a metric or binary variable
sm$p.table
sm$p.pv["dyn_friction_coef"] # extract the p-value in 'sm$p.table' for one variable
# 1b) for a categorical variable with >2 categories
# -> Likelihood ratio test (see above)

### 2) Is the smooth effect different from zero?
# 2a) Globally
sm$s.table
sm$s.pv[2] # extract the p-value in 'sm$s.table' for one variable
# 2b) At a specific point
# -> Use bootstraped confidence intervals
# 2c) In a specific interval
# -> Use bootstrapped confidence intervals

### 3) Is at least one of multiple parameters different from zero?
# -> Likelihood ratio test (see above)

### 4) Are two linear effects different from one another?
# see explanation in paper

### 5) Are two smooth effects different from one another?
# see explanation in paper

### 6) Is the linear effect different depending on another variable?
# 6.1) both variables are metric or binary
#      6.1a) x_k is binary or the effect is varying linearly over the metric x_k
#      -> see hypothesis 1a)
#      6.1b) the effect is varying nonlinearly over the metric x_k
#      -> Likelihood ratio test (see above)
# 6.2) at least one of the variables is categorical with >2 categories
# -> Likelihood ratio test (see above)

### 7) Is the smooth effect different depending on another variable?
# 7a) x_k is binary
# -> see 2a)
# 7b) x_k is metric or categorical with >2 categories
# -> Likelihood ratio test (see above)

### 8) Is one model better than another model?
# -> Likelihood ratio test (see above)




# Model evaluation --------------------------------------------------------
### 1) Residual plots
# Residuals vs fitted values
plot_resVSfitted(model)

# Residuals vs the functional domain (here: time)
plot_resVSyindex(model, xlab = "time [s]",
                 main = "Residuals vs time\nHeatmap of binned points")

# Residuals vs space
# Note: the plot is not as highly resoluted as in the paper as we only work with
# a subset of the real data. For the figure in the paper, in order to show the
# visualization method properly, we predicted not just the subset-data with the
# model, but all seismographs from the original dataset of Bauer (2016).
dat_xy <- prepareData_residsVSxy(model, data, "ground_velocity", 
                                 "seismometer_xCoord", "seismometer_yCoord")
breaks <- c(-1,-0.25,-0.05,0.05,0.25,1)
labels_cut <- c("(-1.00,-0,25]","(-0.25,-0.05]","(-0.05,  0.05]","(  0.05,  0.25]","(  0.25,  1.00]")
plot_residsVSxy(model, dat_xy, breaks, labels_cut,
                mark_location = which(dat_xy$seismometer == 62)) # Mark epicenter


# Residual autocovariance
plot_residAutocov(model)

### 2) Compare predictions with observations
ybreaks <- 10^c(-2,-1,0,1)
plot_predVSobs(model, data[data$simulation == 88 & data$seismometer == 62,],
               yvar = "ground_velocity", yvar_label = "ground velocity [m/s]",
               log10 = TRUE, ylim = 10^c(-2,1), ybreaks = ybreaks, xlab = "time [s]")
