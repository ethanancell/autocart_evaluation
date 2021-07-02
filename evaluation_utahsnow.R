library(autocart)
library(CAST)
library(rpart)
library(randomForest)
library(sp)
library(tidyverse)

# Set a seed for reproducibility
set.seed(1234)

# Function to evaluate RMSE
rmse <- function(y_pred, y_true) {
  sqrt(mean((y_pred - y_true)^2, na.rm = TRUE))
}

# Spatial cross-validation script
source("../autocart_evaluation/cluster_stations.R")

# =========================
# == UTAH 2017 SNOW DATA ==
# =========================

# Set this variable to TRUE if you want to see updates through the LONGITUDE amounts
# of computing that's required for this script (e.g. when various CVs are done)
# Note that when setting show_cv_progress to TRUE, you will see many messages
# in the R output saying "CV XX/30" and "New min RMAE: X:XXX". The
show_ac_progress <- TRUE
show_cv_progress <- TRUE
tune_each_cv <- TRUE

snow <- read_csv("../autocart_evaluation/data/ut2017_snow.csv") %>%
  distinct(LONGITUDE, LATITUDE, .keep_all = TRUE)

# Pick a few relevant features and take out rows with NA observations.
# Then add on cross-validation information as selected with the cluster_stations
# script.
snow <- snow %>%
  select(yr50, LONGITUDE, LATITUDE, ELEVATION, PPTWT, MCMT, MWMT) %>%
  mutate(log_yr50 = log(yr50)) %>%
  select(-yr50) %>%
  na.omit()

snow <- snow %>%
  mutate(fold = cluster_stations(snow$LONGITUDE, snow$LATITUDE, rep(0, nrow(snow)),
                                 dist_adj = 100, elev_adj = 0, h = 3))
# TEST BASIC
# --------------
train_ind <- sample(1:nrow(snow), 150)
snow_train <- SpatialPointsDataFrame(coords = snow[train_ind, ] %>% select(LONGITUDE, LATITUDE),
                                     snow[train_ind, ])
snow_test <- SpatialPointsDataFrame(coords = snow[-train_ind, ] %>% select(LONGITUDE, LATITUDE),
                                    snow[-train_ind, ])

# Tune something?
tune <- autotune(log_yr50 ~ . -fold, data = snow_train, k = 5,
                 alphaVals = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                 betaVals = c(0.0, 0.2), bandwidthVals = c(1),
                 powerVals = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.2, 1.4, 1.6),
                 rangedPowerOffset = 0.25,
                 outputProgress = TRUE)

model <- autocart(log_yr50 ~ . -fold, data = snow_train, tune$alpha, tune$beta)
modelf <- autoforest(log_yr50 ~ . -fold, data = snow_train, tune$alpha, tune$beta)

y_guess1 <- predict(model, snow_test)
y_guess2 <- predict(modelf, snow_test)
y_guess3 <- predict(model, snow_test, spatialNodes = TRUE, p = tune$power)
y_guess4 <- predict(modelf, snow_test, spatialNodes = TRUE, p = tune$power)

res_guess1 <- rmse(y_guess1, snow_test$log_yr50)
res_guess2 <- rmse(y_guess2, snow_test$log_yr50)
res_guess3 <- rmse(y_guess3, snow_test$log_yr50)
res_guess4 <- rmse(y_guess4, snow_test$log_yr50)

# Evaluate accuracy in each of the random forest elements
ind_tree_pred_nsp <- lapply(modelf, predict.autocart, snow_test)
ind_tree_pred <- lapply(modelf, predict.autocart, snow_test, spatialNodes = TRUE, p = tune$power)
ind_tree_acc_nsp <- unlist(lapply(ind_tree_pred_nsp, rmse, snow_test$log_yr50))
ind_tree_acc <- unlist(lapply(ind_tree_pred, rmse, snow_test$log_yr50))

res_guess1
res_guess2
res_guess3
res_guess4

# Inspect a tree
modelf[[1]]

# TEST AUTOTUNE
# --------------
y <- snow$log_yr50
X <- snow %>% dplyr::select(LONGITUDE, LATITUDE, ELEVATION, PPTWT, MCMT, MWMT)
loc <- snow %>% dplyr::select(LONGITUDE, LATITUDE) %>% as.matrix()
autotune(y, X, loc, control = autocartControl(),
         k = 5,
         alphaVals = c(0, 0.33, 0.66, 1.0),
         betaVals = c(0, 0.33),
         bandwidthVals = c(0.5, 1.0),
         powerVals = c(1.5, 2, 2.5, 3),
         outputProgress = TRUE)

# --------------

# These are the predictor variables we use (excluding longitude/latitude
# as those are handled separately)
covariates <- c("ELEVATION", "PPTWT", "MCMT", "MWMT")

# This stores the RMSE for each of the folds
# The abbreviations used below are:
# * rf: random forest
# * rp: rpart (regression trees)
# * ac: autocart
# * af: autoforest
# * w_loc: with location
# * wo_loc: without location
# * nointerp: No interpolation
# * fixp: Fixed value of p
# * rangep: Ranged value of p
fold_rmse_rf_wo_loc <- rep(NA, length(levels(snow$fold)))
fold_rmse_rf_w_loc <- rep(NA, length(levels(snow$fold)))
fold_rmse_rp_wo_loc <- rep(NA, length(levels(snow$fold)))
fold_rmse_rp_w_loc <- rep(NA, length(levels(snow$fold)))
fold_rmse_ac_wo_loc_nointerp <- rep(NA, length(levels(snow$fold)))
fold_rmse_ac_w_loc_nointerp <- rep(NA, length(levels(snow$fold)))
fold_rmse_ac_wo_loc_fixp <- rep(NA, length(levels(snow$fold)))
fold_rmse_ac_w_loc_fixp <- rep(NA, length(levels(snow$fold)))
fold_rmse_ac_wo_loc_rangep <- rep(NA, length(levels(snow$fold)))
fold_rmse_ac_w_loc_rangep <- rep(NA, length(levels(snow$fold)))
fold_rmse_af_wo_loc_nointerp <- rep(NA, length(levels(snow$fold)))
fold_rmse_af_w_loc_nointerp <- rep(NA, length(levels(snow$fold)))
fold_rmse_af_wo_loc_fixp <- rep(NA, length(levels(snow$fold)))
fold_rmse_af_w_loc_fixp <- rep(NA, length(levels(snow$fold)))
fold_rmse_af_wo_loc_rangep <- rep(NA, length(levels(snow$fold)))
fold_rmse_af_w_loc_rangep <- rep(NA, length(levels(snow$fold)))

fold_size <- rep(NA, length(levels(snow$fold)))

# Use these to record progress of the various optimal values of alpha, beta, etc.
opt_alpha_w_loc <- rep(NA, length(levels(snow$fold)))
opt_alpha_wo_loc <- rep(NA, length(levels(snow$fold)))
opt_beta_w_loc <- rep(NA, length(levels(snow$fold)))
opt_beta_wo_loc <- rep(NA, length(levels(snow$fold)))
opt_bandwidth_w_loc <- rep(NA, length(levels(snow$fold)))
opt_bandwidth_wo_loc <- rep(NA, length(levels(snow$fold)))
opt_power_w_loc <- rep(NA, length(levels(snow$fold)))
opt_power_wo_loc <- rep(NA, length(levels(snow$fold)))

num_xval_folds <- length(levels(snow$fold))
for (k in 1:num_xval_folds) {
  n_test <- nrow(snow[snow$fold == k, ])
  fold_size[k] <- n_test

  rf_wo_loc_pred <- rep(NA, n_test)
  rf_w_loc_pred <- rep(NA, n_test)
  rp_wo_loc_pred <- rep(NA, n_test)
  rp_w_loc_pred <- rep(NA, n_test)
  ac_wo_loc_pred_nointerp <- rep(NA, n_test)
  ac_w_loc_pred_nointerp <- rep(NA, n_test)
  ac_wo_loc_pred_fixp <- rep(NA, n_test)
  ac_w_loc_pred_fixp <- rep(NA, n_test)
  ac_wo_loc_pred_rangep <- rep(NA, n_test)
  ac_w_loc_pred_rangep <- rep(NA, n_test)
  af_wo_loc_pred_nointerp <- rep(NA, n_test)
  af_w_loc_pred_nointerp <- rep(NA, n_test)
  af_wo_loc_pred_fixp <- rep(NA, n_test)
  af_w_loc_pred_fixp <- rep(NA, n_test)
  af_wo_loc_pred_rangep <- rep(NA, n_test)
  af_w_loc_pred_rangep <- rep(NA, n_test)

  # Separate into training and testing sets
  train_X <- select(snow[snow$fold != k, ], all_of(covariates))
  train_y <- unlist(select(snow[snow$fold != k, ], log_yr50))
  train_loc <- as.matrix(select(snow[snow$fold != k, ], LONGITUDE, LATITUDE))
  test_X <- select(snow[snow$fold == k, ], all_of(covariates))
  test_y <- unlist(select(snow[snow$fold == k, ], log_yr50))
  test_loc <- as.matrix(select(snow[snow$fold == k, ], LONGITUDE, LATITUDE))

  # Create the data frames for the models that use formula specifications
  # (The above format works okay for autocart which uses matrix specifications)
  train_X_wo_loc <- cbind(log_yr50 = train_y, train_X)
  train_X_w_loc <- cbind(train_X_wo_loc, LONGITUDE = train_loc[, 1], LATITUDE = train_loc[, 2])
  train_X_w_loc_a <- select(train_X_w_loc, -log_yr50)
  test_X_wo_loc <- cbind(log_yr50 = test_y, test_X)
  test_X_w_loc <- cbind(test_X_wo_loc, LONGITUDE = test_loc[, 1], LATITUDE = test_loc[, 2])
  test_X_w_loc_a <- select(test_X_w_loc, -log_yr50)

  # Train our models
  rp_wo_loc <- rpart(log_yr50 ~ ., data = train_X_wo_loc)
  rp_w_loc <- rpart(log_yr50 ~ ., data = train_X_w_loc)
  rf_wo_loc <- randomForest(log_yr50 ~ ., data = train_X_wo_loc, ntree = 50)
  rf_w_loc <- randomForest(log_yr50 ~ ., data = train_X_w_loc, ntree = 50)

  # For the thorough testing procedure, set tune_each_cv to TRUE. If you want
  # quick results without tuning hyperparameters then set tune_each_cv = FALSE
  # above at the top of this script.
  if (tune_each_cv) {
    # Cross-validate good parameters for autocart. Note that this needs to be done
    # only with the training data (and not on the withheld fold) because we don't
    # want our hyperparameters to be tuned on something that it should be blind to.
    alpha_vals <- c(0, 0.33, 0.66, 1.0)
    beta_vals <- c(0, 0.33)
    bandwidth_vals <- c(1)
    power_vals <- c(1.5, 2.0, 2.5, 3.0)

    # Note that we store separate results for both
    ac_res <- autotune(train_y, train_X, train_loc, k = 5, alphaVals = alpha_vals,
                       betaVals = beta_vals, bandwidthVals = bandwidth_vals,
                       powerVals = power_vals, rangedPowerOffset = 0.25,
                       outputProgress = show_cv_progress)
    ac_res_w_loc <- autotune(train_y, train_X_w_loc_a, train_loc, k = 5, alphaVals = alpha_vals,
                             betaVals = beta_vals, bandwidthVals = bandwidth_vals,
                             powerVals = power_vals, rangedPowerOffset = 0.25,
                             outputProgress = show_cv_progress)

    # Store the optimal values of alpha, beta, etc.
    opt_alpha_w_loc[k] <- ac_res_w_loc$alpha
    opt_alpha_wo_loc[k] <- ac_res$alpha
    opt_beta_w_loc[k] <- ac_res_w_loc$beta
    opt_beta_wo_loc[k] <- ac_res$beta
    opt_bandwidth_w_loc[k] <- ac_res_w_loc$bandwidth
    opt_bandwidth_wo_loc[k] <- ac_res$bandwidth

    # Train the autocart + autoforest models
    ac <- autocart(train_y, train_X, train_loc, ac_res$alpha, ac_res$beta,
                   control = autocartControl(spatialBandwidthProportion = ac_res$bandwidth))
    ac_w_loc <- autocart(train_y, train_X_w_loc_a, train_loc, ac_res_w_loc$alpha, ac_res_w_loc$beta,
                         control = autocartControl(spatialBandwidthProportion = ac_res_w_loc$bandwidth))
    af <- autoforest(train_y, train_X, train_loc, ac_res$alpha, ac_res$beta,
                     control = autocartControl(spatialBandwidthProportion = ac_res$bandwidth),
                     numtrees = 50)
    af_w_loc <- autoforest(train_y, train_X_w_loc_a, train_loc, ac_res_w_loc$alpha, ac_res_w_loc$beta,
                           control = autocartControl(spatialBandwidthProportion = ac_res_w_loc$bandwidth),
                           numtrees = 50)

    # Create predictions for the test data
    rf_wo_loc_pred <- predict(rf_wo_loc, test_X_wo_loc)
    rf_w_loc_pred <- predict(rf_w_loc, test_X_w_loc)
    rp_wo_loc_pred <- predict(rp_wo_loc, test_X_wo_loc)
    rp_w_loc_pred <- predict(rp_w_loc, test_X_w_loc)

    ac_wo_loc_pred_nointerp <- predictAutocart(ac, test_X)
    ac_w_loc_pred_nointerp <- predictAutocart(ac_w_loc, test_X_w_loc_a)
    ac_wo_loc_pred_fixp <- spatialNodes(ac, test_X, test_loc, distpower = ac_res$power)
    ac_w_loc_pred_fixp <- spatialNodes(ac_w_loc, test_X_w_loc_a, test_loc, distpower = ac_res$power)
    ac_wo_loc_pred_rangep <- spatialNodes(ac, test_X, test_loc, distpowerRange = c(ac_res$power - 0.25, ac_res$power + 0.25))
    ac_w_loc_pred_rangep <- spatialNodes(ac_w_loc, test_X_w_loc_a, test_loc, distpowerRange = c(ac_res$power - 0.25, ac_res$power + 0.25))

    af_wo_loc_pred_nointerp <- predictAutoforest(af, test_X, test_loc)
    af_w_loc_pred_nointerp <- predictAutoforest(af_w_loc, test_X_w_loc_a, test_loc)
    af_wo_loc_pred_fixp <- predictAutoforest(af, test_X, test_loc,
                                             useSpatialNodes = TRUE,
                                             distpower = ac_res$power)
    af_w_loc_pred_fixp <- predictAutoforest(af_w_loc, test_X_w_loc_a,
                                            test_loc, useSpatialNodes = TRUE,
                                            distpower = ac_res$power)
    af_wo_loc_pred_rangep <- predictAutoforest(af, test_X, test_loc,
                                               useSpatialNodes = TRUE,
                                               distpowerRange = c(ac_res$power - 0.25, ac_res$power + 0.25))
    af_w_loc_pred_rangep <- predictAutoforest(af_w_loc, test_X_w_loc_a,
                                              test_loc, useSpatialNodes = TRUE,
                                              distpowerRange = c(ac_res$power - 0.25, ac_res$power + 0.25))
  } else {
    # Here we just pick some hyperparameters to cut down on computation time
    # for quicker testing results. This is not the procedure we use for the
    # formal testing as used in the final paper.
    use_alpha <- 0.4
    use_beta <- 0.1
    use_bandwidth <- 1.0

    ac <- autocart(train_y, train_X, train_loc, use_alpha, use_beta,
                   control = autocartControl(spatialBandwidthProportion = use_bandwidth))
    ac_w_loc <- autocart(train_y, train_X_w_loc_a, train_loc, use_alpha, use_beta,
                         control = autocartControl(spatialBandwidthProportion = use_bandwidth))
    af <- autoforest(train_y, train_X, train_loc, use_alpha, use_beta,
                     control = autocartControl(spatialBandwidthProportion = use_bandwidth),
                     numtrees = 50)
    af_w_loc <- autoforest(train_y, train_X_w_loc_a, train_loc, use_alpha, use_beta,
                           control = autocartControl(spatialBandwidthProportion = use_bandwidth),
                           numtrees = 50)
    stop("Only doing CV for now.")
  }

  # Calculate RMSE for each fold
  fold_rmse_rf_wo_loc[k] <- rmse(rf_wo_loc_pred, test_y)
  fold_rmse_rf_w_loc[k] <- rmse(rf_w_loc_pred, test_y)
  fold_rmse_rp_wo_loc[k] <- rmse(rp_wo_loc_pred, test_y)
  fold_rmse_rp_w_loc[k] <- rmse(rp_w_loc_pred, test_y)
  fold_rmse_ac_w_loc_nointerp[k] <- rmse(ac_w_loc_pred_nointerp, test_y)
  fold_rmse_ac_wo_loc_nointerp[k] <- rmse(ac_wo_loc_pred_nointerp, test_y)
  fold_rmse_ac_w_loc_fixp[k] <- rmse(ac_w_loc_pred_fixp, test_y)
  fold_rmse_ac_wo_loc_fixp[k] <- rmse(ac_wo_loc_pred_fixp, test_y)
  fold_rmse_ac_w_loc_rangep[k] <- rmse(ac_w_loc_pred_rangep, test_y)
  fold_rmse_ac_wo_loc_rangep[k] <- rmse(ac_wo_loc_pred_rangep, test_y)
  fold_rmse_af_w_loc_nointerp[k] <- rmse(af_w_loc_pred_nointerp, test_y)
  fold_rmse_af_wo_loc_nointerp[k] <- rmse(af_wo_loc_pred_nointerp, test_y)
  fold_rmse_af_w_loc_fixp[k] <- rmse(af_w_loc_pred_fixp, test_y)
  fold_rmse_af_wo_loc_fixp[k] <- rmse(af_wo_loc_pred_fixp, test_y)
  fold_rmse_af_w_loc_rangep[k] <- rmse(af_w_loc_pred_rangep, test_y)
  fold_rmse_af_wo_loc_rangep[k] <- rmse(af_wo_loc_pred_rangep, test_y)

  if (show_ac_progress) {
    print(paste0("Finished fold ", k, "/", num_xval_folds))
  }
}

# Look at results
total_obs <- sum(fold_size)
sum(fold_rmse_rf_wo_loc * fold_size) / total_obs
sum(fold_rmse_rf_w_loc * fold_size) / total_obs
sum(fold_rmse_rp_wo_loc * fold_size) / total_obs
sum(fold_rmse_rp_w_loc * fold_size) / total_obs

sum(fold_rmse_ac_wo_loc_nointerp * fold_size) / total_obs
sum(fold_rmse_ac_w_loc_nointerp * fold_size) / total_obs
sum(fold_rmse_ac_wo_loc_fixp * fold_size) / total_obs
sum(fold_rmse_ac_w_loc_fixp * fold_size) / total_obs
sum(fold_rmse_ac_wo_loc_rangep * fold_size) / total_obs
sum(fold_rmse_ac_w_loc_rangep * fold_size) / total_obs

sum(fold_rmse_af_wo_loc_nointerp * fold_size) / total_obs
sum(fold_rmse_af_w_loc_nointerp * fold_size) / total_obs
sum(fold_rmse_af_wo_loc_fixp * fold_size) / total_obs
sum(fold_rmse_af_w_loc_fixp * fold_size) / total_obs
sum(fold_rmse_af_wo_loc_rangep * fold_size) / total_obs
sum(fold_rmse_af_w_loc_rangep * fold_size) / total_obs

# one run saved in 850obs.RData
# save.image(file = "C:/Users/ethan/code/autocart_evaluation/600obs.RData")
