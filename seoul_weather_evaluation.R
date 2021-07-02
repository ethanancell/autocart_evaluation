library(gstat)
library(rpart)
library(randomForest)
library(sp)
library(tidyverse)

# Set a seed for reproducibility
set.seed(1234)

rmse <- function(y_pred, y_true) {
  sqrt(mean((y_pred - y_true)^2))
}

# =============================
# == Weather Bias Correction ==
# =============================

bc <- read_csv("../autocart_evaluation/data/Bias_correction_ucl.csv") %>%
  rename(solar_radiation = `Solar radiation`)

# Covariates are long, lat, elevation, slope, solar radiation
# predicting:
# * Present_Tmax
# * Present_Tmin
# Split into all the days
bc <- bc %>% select(Date, Present_Tmax, Present_Tmin, lon, lat, DEM, solar_radiation, Slope) %>%
  na.omit()

# Now we need to split into individual days
all_days <- sort(unique(bc$Date))
start_day <- all_days[1]
test_day <- all_days[20]
end_day <- all_days[length(all_days)]
day <- start_day
num_days <- as.integer(end_day - start_day)
day_index <- 1

# Tune some parameters: These histories allow us to get a sense for what values
# of the parameters work well.
day <- start_day

num_test_days <- as.integer(test_day - start_day)

alpha_wl_history_tmax <- rep(NA, num_test_days)
beta_wl_history_tmax <- rep(NA, num_test_days)
bandwidth_wl_history_tmax <- rep(NA, num_test_days)
power_wl_history_tmax <- rep(NA, num_test_days)
powerRange_wl_history_tmax <- rep(NA, num_test_days)
alpha_wol_history_tmax <- rep(NA, num_test_days)
beta_wol_history_tmax <- rep(NA, num_test_days)
bandwidth_wol_history_tmax <- rep(NA, num_test_days)
power_wol_history_tmax <- rep(NA, num_test_days)
powerRange_wol_history_tmax <- rep(NA, num_test_days)

alpha_wl_history_tmin <- rep(NA, num_test_days)
beta_wl_history_tmin <- rep(NA, num_test_days)
bandwidth_wl_history_tmin <- rep(NA, num_test_days)
power_wl_history_tmin <- rep(NA, num_test_days)
powerRange_wl_history_tmin <- rep(NA, num_test_days)
alpha_wol_history_tmin <- rep(NA, num_test_days)
beta_wol_history_tmin <- rep(NA, num_test_days)
bandwidth_wol_history_tmin <- rep(NA, num_test_days)
power_wol_history_tmin <- rep(NA, num_test_days)
powerRange_wol_history_tmin <- rep(NA, num_test_days)

day_index <- 1
while (day <= test_day) {
  paste(paste0("Tuning on day ", day))
  bc_thisday <- bc[bc$Date == day, ] %>% select(-Date)
  bc_thisday <- SpatialPointsDataFrame(coords = bc_thisday %>% select(lon, lat),
                                       data = bc_thisday)
  num_obs_this_day <- nrow(bc_thisday)

  # Autocart Tuning
  tune_w_loc_tmax <- autotune(Present_Tmax ~ . -Present_Tmin, data = bc_thisday, k=5,
                              alphaVals = seq(0, 1, by = 0.20),
                              betaVals = c(0, 0.05, 0.10),bandwidthVals = c(1.0),
                              powerVals = c(0, 0.3, 0.5, 1.0, 1.5),rangedPowerOffset = c(0, 0.10, 0.25),
                              autocartControl(minsplit = 8),
                              outputProgress = TRUE)
  tune_wo_loc_tmax <- autotune(Present_Tmax ~ . -Present_Tmin -lon -lat, data = bc_thisday, k=5,
                               alphaVals = seq(0, 1, by = 0.20),betaVals = c(0, 0.05, 0.10),
                               bandwidthVals = c(1.0),
                               powerVals = c(0, 0.3, 0.5, 1.0, 1.5),
                               rangedPowerOffset = c(0, 0.10, 0.25),
                               autocartControl(minsplit = 8),
                               outputProgress = TRUE)
  tune_w_loc_tmin <- autotune(Present_Tmin ~ . -Present_Tmax, data = bc_thisday, k=5,
                              alphaVals = seq(0, 1, by = 0.20),
                              betaVals = c(0, 0.05, 0.10),bandwidthVals = c(1.0),
                              powerVals = c(0, 0.3, 0.5, 1.0, 1.5),rangedPowerOffset = c(0, 0.10, 0.25),
                              autocartControl(minsplit = 8),
                              outputProgress = TRUE)
  tune_wo_loc_tmin <- autotune(Present_Tmin ~ . -Present_Tmax -lon -lat, data = bc_thisday, k=5,
                               alphaVals = seq(0, 1, by = 0.20),betaVals = c(0, 0.05, 0.10),
                               bandwidthVals = c(1.0),
                               powerVals = c(0, 0.3, 0.5, 1.0, 1.5),
                               rangedPowerOffset = c(0, 0.10, 0.25),
                               autocartControl(minsplit = 8),
                               outputProgress = TRUE)

  alpha_wl_history_tmax[day_index] <- tune_w_loc_tmax$alpha
  beta_wl_history_tmax[day_index] <- tune_w_loc_tmax$beta
  bandwidth_wl_history_tmax[day_index] <- tune_w_loc_tmax$bandwidth
  power_wl_history_tmax[day_index] <- tune_w_loc_tmax$power
  powerRange_wl_history_tmax[day_index] <- tune_w_loc_tmax$powerOffset
  alpha_wol_history_tmax[day_index] <- tune_wo_loc_tmax$alpha
  beta_wol_history_tmax[day_index] <- tune_wo_loc_tmax$beta
  bandwidth_wol_history_tmax[day_index] <- tune_wo_loc_tmax$bandwidth
  power_wol_history_tmax[day_index] <- tune_wo_loc_tmax$power
  powerRange_wol_history_tmax[day_index] <- tune_wo_loc_tmax$powerOffset

  alpha_wl_history_tmin[day_index] <- tune_w_loc_tmin$alpha
  beta_wl_history_tmin[day_index] <- tune_w_loc_tmin$beta
  bandwidth_wl_history_tmin[day_index] <- tune_w_loc_tmin$bandwidth
  power_wl_history_tmin[day_index] <- tune_w_loc_tmin$power
  powerRange_wl_history_tmin[day_index] <- tune_w_loc_tmin$powerOffset
  alpha_wol_history_tmin[day_index] <- tune_wo_loc_tmin$alpha
  beta_wol_history_tmin[day_index] <- tune_wo_loc_tmin$beta
  bandwidth_wol_history_tmin[day_index] <- tune_wo_loc_tmin$bandwidth
  power_wol_history_tmin[day_index] <- tune_wo_loc_tmin$power
  powerRange_wol_history_tmin[day_index] <- tune_wo_loc_tmin$powerOffset

  day_index <- day_index + 1
  day <- day + 1
}

# Plots to look at histories
plot(alpha_wl_history_tmax)
lines(alpha_wl_history_tmax)
plot(alpha_wol_history_tmax)
lines(alpha_wol_history_tmax)
plot(beta_wl_history_tmax)
lines(beta_wl_history_tmax)
plot(beta_wol_history_tmax)
lines(beta_wol_history_tmax)
plot(bandwidth_wl_history_tmax)
lines(bandwidth_wl_history_tmax)
plot(bandwidth_wol_history_tmax)
lines(bandwidth_wol_history_tmax)
plot(power_wl_history_tmax)
lines(power_wl_history_tmax)
plot(power_wol_history_tmax)
lines(power_wol_history_tmax)
plot(powerRange_wl_history_tmax)
lines(powerRange_wl_history_tmax)
plot(powerRange_wol_history_tmax)
lines(powerRange_wol_history_tmax)

plot(alpha_wl_history_tmin)
lines(alpha_wl_history_tmin)
plot(alpha_wol_history_tmin)
lines(alpha_wol_history_tmin)
plot(beta_wl_history_tmin)
lines(beta_wl_history_tmin)
plot(beta_wol_history_tmin)
lines(beta_wol_history_tmin)
plot(bandwidth_wl_history_tmin)
lines(bandwidth_wl_history_tmin)
plot(bandwidth_wol_history_tmin)
lines(bandwidth_wol_history_tmin)
plot(power_wl_history_tmin)
lines(power_wl_history_tmin)
plot(power_wol_history_tmin)
lines(power_wol_history_tmin)
plot(powerRange_wl_history_tmin)
lines(powerRange_wl_history_tmin)
plot(powerRange_wol_history_tmin)
lines(powerRange_wol_history_tmin)

use_alpha_wl_tmax <- mean(alpha_wl_history_tmax)
use_alpha_wol_tmax <- mean(alpha_wol_history_tmax)
use_beta_wl_tmax <- mean(beta_wl_history_tmax)
use_beta_wol_tmax <- mean(beta_wol_history_tmax)
use_bandwidth_wl_tmax <- mean(bandwidth_wl_history_tmax)
use_bandwidth_wol_tmax <- mean(bandwidth_wol_history_tmax)
use_power_wl_tmax <- mean(power_wl_history_tmax)
use_power_wol_tmax <- mean(power_wol_history_tmax)
use_powerRange_wl_tmax <- mean(powerRange_wl_history_tmax)
use_powerRange_wol_tmax <- mean(powerRange_wol_history_tmax)

use_alpha_wl_tmin <- mean(alpha_wl_history_tmin)
use_alpha_wol_tmin <- mean(alpha_wol_history_tmin)
use_beta_wl_tmin <- mean(beta_wl_history_tmin)
use_beta_wol_tmin <- mean(beta_wol_history_tmin)
use_bandwidth_wl_tmin <- mean(bandwidth_wl_history_tmin)
use_bandwidth_wol_tmin <- mean(bandwidth_wol_history_tmin)
use_power_wl_tmin <- mean(power_wl_history_tmin)
use_power_wol_tmin <- mean(power_wol_history_tmin)
use_powerRange_wl_tmin <- mean(powerRange_wl_history_tmin)
use_powerRange_wol_tmin <- mean(powerRange_wol_history_tmin)

# TEST PARAMETERS ON REST OF DAYS
# -----------------------------------
# Set a seed again: we want a reproducible script below without needing to
# run everything above to get the right randomization state. (this made it
# easier for me personally when writing this code.)
set.seed(1234)

# Test these parameters on some other data
# end_day <- all_days[length(all_days)]
end_day <- all_days[50]
start_day <- test_day + 1
test_days <- all_days[start_day <= all_days & all_days <= end_day]
num_test_days <- length(test_days)

rmse_ac_w_loc_tmax <- rep(NA, num_test_days)
rmse_ac_wo_loc_tmax <- rep(NA, num_test_days)
rmse_ac_w_loc_rangep_tmax <- rep(NA, num_test_days)
rmse_ac_wo_loc_rangep_tmax <- rep(NA, num_test_days)
rmse_af_w_loc_tmax <- rep(NA, num_test_days)
rmse_af_wo_loc_tmax <- rep(NA, num_test_days)
rmse_af_w_loc_rangep_tmax <- rep(NA, num_test_days)
rmse_af_wo_loc_rangep_tmax <- rep(NA, num_test_days)
rmse_rp_w_loc_tmax <- rep(NA, num_test_days)
rmse_rp_wo_loc_tmax <- rep(NA, num_test_days)
rmse_rf_w_loc_tmax <- rep(NA, num_test_days)
rmse_rf_wo_loc_tmax <- rep(NA, num_test_days)
rmse_idw_tmax <- rep(NA, num_test_days)
rmse_idw_drift_tmax <- rep(NA, num_test_days) # IDW with regression
rmse_null_tmax <- rep(NA, num_test_days)

rmse_ac_w_loc_tmin <- rep(NA, num_test_days)
rmse_ac_wo_loc_tmin <- rep(NA, num_test_days)
rmse_ac_w_loc_rangep_tmin <- rep(NA, num_test_days)
rmse_ac_wo_loc_rangep_tmin <- rep(NA, num_test_days)
rmse_af_w_loc_tmin <- rep(NA, num_test_days)
rmse_af_wo_loc_tmin <- rep(NA, num_test_days)
rmse_af_w_loc_rangep_tmin <- rep(NA, num_test_days)
rmse_af_wo_loc_rangep_tmin <- rep(NA, num_test_days)
rmse_rp_w_loc_tmin <- rep(NA, num_test_days)
rmse_rp_wo_loc_tmin <- rep(NA, num_test_days)
rmse_rf_w_loc_tmin <- rep(NA, num_test_days)
rmse_rf_wo_loc_tmin <- rep(NA, num_test_days)
rmse_idw_tmin <- rep(NA, num_test_days)
rmse_idw_drift_tmin <- rep(NA, num_test_days)
rmse_null_tmin <- rep(NA, num_test_days)

day_index <- 1
while (day_index <= num_test_days) {

  # Set an appropriate day based upon the day index
  day <- test_days[day_index]

  print(paste0("Day ", day, ": ", 100 * as.integer(day_index - num_test_days) / num_test_days, "%"))

  bc_thisday <- bc[bc$Date == day, ] %>% select(-Date)
  bc_thisday <- SpatialPointsDataFrame(coords = bc_thisday %>% select(lon, lat),
                                       data = bc_thisday)
  proj4string(bc_thisday) <- "+proj=longlat"
  num_obs_this_day <- nrow(bc_thisday)

  Tmax_guess_ac_w_loc <- rep(NA, num_obs_this_day)
  Tmax_guess_ac_wo_loc <- rep(NA, num_obs_this_day)
  Tmax_guess_ac_w_loc_rangep <- rep(NA, num_obs_this_day)
  Tmax_guess_ac_wo_loc_rangep <- rep(NA, num_obs_this_day)
  Tmax_guess_af_w_loc <- rep(NA, num_obs_this_day)
  Tmax_guess_af_wo_loc <- rep(NA, num_obs_this_day)
  Tmax_guess_af_w_loc_rangep <- rep(NA, num_obs_this_day)
  Tmax_guess_af_wo_loc_rangep <- rep(NA, num_obs_this_day)
  Tmax_guess_rp_w_loc <- rep(NA, num_obs_this_day)
  Tmax_guess_rp_wo_loc <- rep(NA, num_obs_this_day)
  Tmax_guess_rf_w_loc <- rep(NA, num_obs_this_day)
  Tmax_guess_rf_wo_loc <- rep(NA, num_obs_this_day)
  Tmax_guess_idw <- rep(NA, num_obs_this_day)
  Tmax_guess_idw_drift <- rep(NA, num_obs_this_day)
  Tmax_guess_null <- rep(NA, num_obs_this_day)

  Tmin_guess_ac_w_loc <- rep(NA, num_obs_this_day)
  Tmin_guess_ac_wo_loc <- rep(NA, num_obs_this_day)
  Tmin_guess_ac_w_loc_rangep <- rep(NA, num_obs_this_day)
  Tmin_guess_ac_wo_loc_rangep <- rep(NA, num_obs_this_day)
  Tmin_guess_af_w_loc <- rep(NA, num_obs_this_day)
  Tmin_guess_af_wo_loc <- rep(NA, num_obs_this_day)
  Tmin_guess_af_w_loc_rangep <- rep(NA, num_obs_this_day)
  Tmin_guess_af_wo_loc_rangep <- rep(NA, num_obs_this_day)
  Tmin_guess_rp_w_loc <- rep(NA, num_obs_this_day)
  Tmin_guess_rp_wo_loc <- rep(NA, num_obs_this_day)
  Tmin_guess_rf_w_loc <- rep(NA, num_obs_this_day)
  Tmin_guess_rf_wo_loc <- rep(NA, num_obs_this_day)
  Tmin_guess_idw <- rep(NA, num_obs_this_day)
  Tmin_guess_idw_drift <- rep(NA, num_obs_this_day)
  Tmin_guess_null <- rep(NA, num_obs_this_day)

  # Do some LOOCV
  for (i in 1:num_obs_this_day) {
    bc_train <- bc_thisday[-i, ]
    bc_test <- bc_thisday[i, ]

    # IDW and null model
    idw_tmax <- gstat(formula = Present_Tmax ~ 1, data = bc_train)
    idw_tmin <- gstat(formula = Present_Tmin ~ 1, data = bc_train)
    # idw_tmax_drift <- gstat(formula = Present_Tmax ~ lon + lat + solar_radiation + DEM + Slope, data = bc_train)
    # idw_tmin_drift <- gstat(formula = Present_Tmin ~ lon + lat + solar_radiation + DEM + Slope, data = bc_train)

    Tmax_guess_idw[i] <- predict(idw_tmax, bc_test)$var1.pred
    # Tmax_guess_idw_drift[i] <- predict(idw_tmax_drift, bc_test)$var1.pred
    Tmin_guess_idw[i] <- predict(idw_tmin, bc_test)$var1.pred
    # Tmin_guess_idw_drift[i] <- predict(idw_tmin_drift, bc_test)$var1.pred

    Tmax_guess_null[i] <- mean(bc_train$Present_Tmax)
    Tmin_guess_null[i] <- mean(bc_train$Present_Tmin)

    # Machine learning models
    ac_w_loc_tmax <- autocart(Present_Tmax ~ . -Present_Tmin, data = bc_train,
                              use_alpha_wl_tmax, use_beta_wl_tmax, control = autocartControl(minsplit = 15))
    ac_wo_loc_tmax <- autocart(Present_Tmax ~ . -Present_Tmin -lon -lat, data = bc_train,
                               use_alpha_wol_tmax, use_beta_wol_tmax, control = autocartControl(minsplit = 15))
    af_w_loc_tmax <- autoforest(Present_Tmax ~ . -Present_Tmin, data = bc_train,
                                use_alpha_wl_tmax, use_beta_wl_tmax, numtrees = 15, control = autocartControl(minsplit = 15))

    af_wo_loc_tmax <- autoforest(Present_Tmax ~ . -Present_Tmin -lon -lat, data = bc_train,
                                 use_alpha_wol_tmax, use_beta_wol_tmax, numtrees = 15, control = autocartControl(minsplit = 15))

    rp_w_loc_tmax <- rpart(Present_Tmax ~ . -Present_Tmin, data = bc_train, control = rpart.control(minsplit = 15))
    rp_wo_loc_tmax <- rpart(Present_Tmax ~ . -Present_Tmin -lon -lat, data = bc_train, control = rpart.control(minsplit = 15))
    rf_w_loc_tmax <- randomForest(Present_Tmax ~ . -Present_Tmin, data = bc_train, ntree = 15)
    rf_wo_loc_tmax <- randomForest(Present_Tmax ~ . -Present_Tmin -lon -lat, data = bc_train, ntree = 15)

    ac_w_loc_tmin <- autocart(Present_Tmin ~ . -Present_Tmax, data = bc_train,
                              use_alpha_wl_tmin, use_beta_wl_tmin, control = autocartControl(minsplit = 15))
    ac_wo_loc_tmin <- autocart(Present_Tmin ~ . -Present_Tmax -lon -lat, data = bc_train,
                               use_alpha_wol_tmin, use_beta_wol_tmin, control = autocartControl(minsplit = 15))
    af_w_loc_tmin <- autoforest(Present_Tmin ~ . -Present_Tmax, data = bc_train,
                                use_alpha_wl_tmin, use_beta_wl_tmin, numtrees = 15, control = autocartControl(minsplit = 15))

    af_wo_loc_tmin <- autoforest(Present_Tmin ~ . -Present_Tmax -lon -lat, data = bc_train,
                                 use_alpha_wol_tmin, use_beta_wol_tmin, numtrees = 15, control = autocartControl(minsplit = 15))

    rp_w_loc_tmin <- rpart(Present_Tmin ~ . -Present_Tmax, data = bc_train, control = rpart.control(minsplit = 15))
    rp_wo_loc_tmin <- rpart(Present_Tmin ~ . -Present_Tmax -lon -lat, data = bc_train, control = rpart.control(minsplit = 15))
    rf_w_loc_tmin <- randomForest(Present_Tmin ~ . -Present_Tmax, data = bc_train, ntree = 15)
    rf_wo_loc_tmin <- randomForest(Present_Tmin ~ . -Present_Tmax -lon -lat, data = bc_train, ntree = 15)

    # Predict for the left out day
    Tmax_guess_ac_w_loc[i] <- predict(ac_w_loc_tmax, bc_test, spatialNodes = TRUE, p = use_power_wl_tmax)
    Tmax_guess_ac_wo_loc[i] <- predict(ac_wo_loc_tmax, bc_test, spatialNodes = TRUE, p = use_power_wol_tmax)
    Tmax_guess_ac_w_loc_rangep[i] <- predict(ac_w_loc_tmax, bc_test, spatialNodes = TRUE,
                                             pRange = c(use_power_wl_tmax-use_powerRange_wl_tmax, use_power_wl_tmax+use_powerRange_wl_tmax))
    Tmax_guess_ac_wo_loc_rangep[i] <- predict(ac_wo_loc_tmax, bc_test, spatialNodes = TRUE,
                                              pRange = c(use_power_wol_tmax-use_powerRange_wol_tmax, use_power_wol_tmax+use_powerRange_wol_tmax))
    Tmax_guess_af_w_loc[i] <- predict(af_w_loc_tmax, bc_test, spatialNodes = TRUE, p = use_power_wl_tmax)
    Tmax_guess_af_wo_loc[i] <- predict(af_wo_loc_tmax, bc_test, spatialNodes = TRUE, p = use_power_wol_tmax)
    Tmax_guess_af_w_loc_rangep[i] <- predict(af_w_loc_tmax, bc_test, spatialNodes = TRUE,
                                             pRange = c(use_power_wl_tmax-use_powerRange_wl_tmax, use_power_wl_tmax+use_powerRange_wl_tmax))
    Tmax_guess_af_wo_loc_rangep[i] <- predict(af_wo_loc_tmax, bc_test, spatialNodes = TRUE,
                                              pRange = c(use_power_wol_tmax-use_powerRange_wol_tmax, use_power_wol_tmax+use_powerRange_wol_tmax))
    Tmax_guess_rp_w_loc[i] <- predict(rp_w_loc_tmax, bc_test)
    Tmax_guess_rp_wo_loc[i] <- predict(rp_wo_loc_tmax, bc_test)
    Tmax_guess_rf_w_loc[i] <- predict(rf_w_loc_tmax, bc_test)
    Tmax_guess_rf_wo_loc[i] <- predict(rf_wo_loc_tmax, bc_test)

    Tmin_guess_ac_w_loc[i] <- predict(ac_w_loc_tmin, bc_test, spatialNodes = TRUE, p = use_power_wl_tmin)
    Tmin_guess_ac_wo_loc[i] <- predict(ac_wo_loc_tmin, bc_test, spatialNodes = TRUE, p = use_power_wol_tmin)
    Tmin_guess_ac_w_loc_rangep[i] <- predict(ac_w_loc_tmin, bc_test, spatialNodes = TRUE,
                                             pRange = c(use_power_wl_tmin-use_powerRange_wl_tmin, use_power_wl_tmin+use_powerRange_wl_tmin))
    Tmin_guess_ac_wo_loc_rangep[i] <- predict(ac_wo_loc_tmin, bc_test, spatialNodes = TRUE,
                                              pRange = c(use_power_wol_tmin-use_powerRange_wol_tmin, use_power_wol_tmin+use_powerRange_wol_tmin))
    Tmin_guess_af_w_loc[i] <- predict(af_w_loc_tmin, bc_test, spatialNodes = TRUE, p = use_power_wl_tmin)
    Tmin_guess_af_wo_loc[i] <- predict(af_wo_loc_tmin, bc_test, spatialNodes = TRUE, p = use_power_wol_tmin)
    Tmin_guess_af_w_loc_rangep[i] <- predict(af_w_loc_tmin, bc_test, spatialNodes = TRUE,
                                             pRange = c(use_power_wl_tmin-use_powerRange_wl_tmin, use_power_wl_tmin+use_powerRange_wl_tmin))
    Tmin_guess_af_wo_loc_rangep[i] <- predict(af_wo_loc_tmin, bc_test, spatialNodes = TRUE,
                                              pRange = c(use_power_wol_tmin-use_powerRange_wol_tmin, use_power_wol_tmin+use_powerRange_wol_tmin))
    Tmin_guess_rp_w_loc[i] <- predict(rp_w_loc_tmin, bc_test)
    Tmin_guess_rp_wo_loc[i] <- predict(rp_wo_loc_tmin, bc_test)
    Tmin_guess_rf_w_loc[i] <- predict(rf_w_loc_tmin, bc_test)
    Tmin_guess_rf_wo_loc[i] <- predict(rf_wo_loc_tmin, bc_test)
  }

  # Accuracy checks
  rmse_ac_w_loc_tmax[day_index] <- rmse(Tmax_guess_ac_w_loc, bc_thisday$Present_Tmax)
  rmse_ac_wo_loc_tmax[day_index] <- rmse(Tmax_guess_ac_wo_loc, bc_thisday$Present_Tmax)
  rmse_ac_w_loc_rangep_tmax[day_index] <- rmse(Tmax_guess_ac_w_loc_rangep, bc_thisday$Present_Tmax)
  rmse_ac_wo_loc_rangep_tmax[day_index] <- rmse(Tmax_guess_ac_wo_loc_rangep, bc_thisday$Present_Tmax)
  rmse_af_w_loc_tmax[day_index] <- rmse(Tmax_guess_af_w_loc, bc_thisday$Present_Tmax)
  rmse_af_wo_loc_tmax[day_index] <- rmse(Tmax_guess_af_wo_loc, bc_thisday$Present_Tmax)
  rmse_af_w_loc_rangep_tmax[day_index] <- rmse(Tmax_guess_af_w_loc_rangep, bc_thisday$Present_Tmax)
  rmse_af_wo_loc_rangep_tmax[day_index] <- rmse(Tmax_guess_af_wo_loc_rangep, bc_thisday$Present_Tmax)
  rmse_rp_w_loc_tmax[day_index] <- rmse(Tmax_guess_rp_w_loc, bc_thisday$Present_Tmax)
  rmse_rp_wo_loc_tmax[day_index] <- rmse(Tmax_guess_rp_wo_loc, bc_thisday$Present_Tmax)
  rmse_rf_w_loc_tmax[day_index] <- rmse(Tmax_guess_rf_w_loc, bc_thisday$Present_Tmax)
  rmse_rf_wo_loc_tmax[day_index] <- rmse(Tmax_guess_rf_wo_loc, bc_thisday$Present_Tmax)
  rmse_idw_tmax[day_index] <- rmse(Tmax_guess_idw, bc_thisday$Present_Tmax)
  rmse_null_tmax[day_index] <- rmse(Tmax_guess_null, bc_thisday$Present_Tmax)

  rmse_ac_w_loc_tmin[day_index] <- rmse(Tmin_guess_ac_w_loc, bc_thisday$Present_Tmin)
  rmse_ac_wo_loc_tmin[day_index] <- rmse(Tmin_guess_ac_wo_loc, bc_thisday$Present_Tmin)
  rmse_ac_w_loc_rangep_tmin[day_index] <- rmse(Tmin_guess_ac_w_loc_rangep, bc_thisday$Present_Tmin)
  rmse_ac_wo_loc_rangep_tmin[day_index] <- rmse(Tmin_guess_ac_wo_loc_rangep, bc_thisday$Present_Tmin)
  rmse_af_w_loc_tmin[day_index] <- rmse(Tmin_guess_af_w_loc, bc_thisday$Present_Tmin)
  rmse_af_wo_loc_tmin[day_index] <- rmse(Tmin_guess_af_wo_loc, bc_thisday$Present_Tmin)
  rmse_af_w_loc_rangep_tmin[day_index] <- rmse(Tmin_guess_af_w_loc_rangep, bc_thisday$Present_Tmin)
  rmse_af_wo_loc_rangep_tmin[day_index] <- rmse(Tmin_guess_af_wo_loc_rangep, bc_thisday$Present_Tmin)
  rmse_rp_w_loc_tmin[day_index] <- rmse(Tmin_guess_rp_w_loc, bc_thisday$Present_Tmin)
  rmse_rp_wo_loc_tmin[day_index] <- rmse(Tmin_guess_rp_wo_loc, bc_thisday$Present_Tmin)
  rmse_rf_w_loc_tmin[day_index] <- rmse(Tmin_guess_rf_w_loc, bc_thisday$Present_Tmin)
  rmse_rf_wo_loc_tmin[day_index] <- rmse(Tmin_guess_rf_wo_loc, bc_thisday$Present_Tmin)
  rmse_idw_tmin[day_index] <- rmse(Tmin_guess_idw, bc_thisday$Present_Tmin)
  rmse_null_tmin[day_index] <- rmse(Tmin_guess_null, bc_thisday$Present_Tmin)

  day_index <- day_index + 1
}

# Look at the results
mean(rmse_ac_w_loc_tmax, na.rm = TRUE)
mean(rmse_ac_wo_loc_tmax, na.rm = TRUE)
mean(rmse_ac_w_loc_rangep_tmax, na.rm = TRUE)
mean(rmse_ac_wo_loc_rangep_tmax, na.rm = TRUE)
mean(rmse_af_w_loc_tmax, na.rm = TRUE)
mean(rmse_af_wo_loc_tmax, na.rm = TRUE)
mean(rmse_af_w_loc_rangep_tmax, na.rm = TRUE)
mean(rmse_af_wo_loc_rangep_tmax, na.rm = TRUE)
mean(rmse_rp_w_loc_tmax, na.rm = TRUE)
mean(rmse_rp_wo_loc_tmax, na.rm = TRUE)
mean(rmse_rf_w_loc_tmax, na.rm = TRUE)
mean(rmse_rf_wo_loc_tmax, na.rm = TRUE)
mean(rmse_idw_tmax, na.rm = TRUE)
mean(rmse_null_tmax, na.rm = TRUE)

mean(rmse_ac_w_loc_tmin, na.rm = TRUE)
mean(rmse_ac_wo_loc_tmin, na.rm = TRUE)
mean(rmse_ac_w_loc_rangep_tmin, na.rm = TRUE)
mean(rmse_ac_wo_loc_rangep_tmin, na.rm = TRUE)
mean(rmse_af_w_loc_tmin, na.rm = TRUE)
mean(rmse_af_wo_loc_tmin, na.rm = TRUE)
mean(rmse_af_w_loc_rangep_tmin, na.rm = TRUE)
mean(rmse_af_wo_loc_rangep_tmin, na.rm = TRUE)
mean(rmse_rp_w_loc_tmin, na.rm = TRUE)
mean(rmse_rp_wo_loc_tmin, na.rm = TRUE)
mean(rmse_rf_w_loc_tmin, na.rm = TRUE)
mean(rmse_rf_wo_loc_tmin, na.rm = TRUE)
mean(rmse_idw_tmin, na.rm = TRUE)
mean(rmse_null_tmin, na.rm = TRUE)

# Box plots of results
result_df <- data.frame(rmse_ac_w_loc, rmse_ac_wo_loc, rmse_ac_w_loc_rangep,
                        rmse_ac_wo_loc_rangep, rmse_af_w_loc, rmse_af_wo_loc,
                        rmse_af_w_loc_rangep, rmse_af_wo_loc_rangep,
                        rmse_rp_w_loc, rmse_rp_wo_loc, rmse_rf_w_loc, rmse_rf_wo_loc)
ggplot(stack(result_df), aes(x = ind, y = values)) +
  geom_boxplot()
