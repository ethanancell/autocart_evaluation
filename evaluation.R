library(autocart)
library(CAST)
library(rpart)
library(randomForest)
library(lubridate)
library(tidyverse)

# Set a seed for reproducibility
set.seed(1234)

# Function to evaluate RMSE
rmse <- function(y_pred, y_true) {
  sqrt(mean((y_pred - y_true)^2, na.rm = TRUE))
}

# ============================
# == KING COUNTY HOUSE DATA ==
# ============================

# Set this variable to TRUE if you want to see updates through the long amounts
# of computing that's required for this script (e.g. when various CVs are done)
# Note that when setting show_cv_progress to TRUE, you will see many messages
# in the R output saying "CV XX/30" and "New min RMAE: X:XXX". The
show_ac_progress <- TRUE
show_cv_progress <- TRUE
tune_each_cv <- FALSE

# Only keep the April - August months to reduce the effect of temporal trend
kc <- read_csv("../autocart_evaluation/data/kc_house_data.csv") %>%
  distinct(long, lat, .keep_all = TRUE)

months <- month(kc$date)
valid_months <- months >= 4 & months <= 8
kc <- kc[valid_months, ]

# Take out any repeat locations (there should be very few of these though,
# they would exist from repeat listings.)

# Subset the data to make autocart run in a feasible amount of time
# kc <- kc[sample(1:nrow(kc), 2000), ]

# Pick a few relevant features and take out rows with NA observations
kc <- kc %>%
  select(price, bedrooms, bathrooms, sqft_living, sqft_lot, floors,
         yr_built, long, lat) %>%
  #mutate(price_p_sqft = price / sqft_living)
  mutate(log_price = log(price))
kc <- na.omit(kc)

# Create folds (with cluster_stations()) and then evaluate
source("../autocart_evaluation/cluster_stations.R")
kc <- kc %>% mutate(fold = cluster_stations(kc$long, kc$lat, rep(0, nrow(kc)),
                                            dist_adj = 7, elev_adj = 0, h = 3))
kc <- SpatialPointsDataFrame(kc %>% dplyr::select(long, lat), kc)
# These are the predictor variables we use (excluding longitude/latitude
# as those are handled separately)
covariates <- c("bedrooms", "bathrooms", "sqft_living", "sqft_lot", "floors",
                "yr_built")

# TEST STUFF
# ------------
train_index <- sample(1:nrow(kc), 500)
kc_train <- kc[train_index, ]
kc_test <- kc[-train_index, ]

# Without coordinates
tune <- autotune(log_price ~ . -long -lat -price -fold, data = kc_train,
                 k=4, alphaVals = c(0, 0.33, 0.66, 1.0),
                 betaVals = c(0, 0.33), bandwidthVals = c(1.0),
                 powerVals = c(0, 1.0, 1.5, 2.0, 2.5),
                 outputProgress = TRUE)

# w/ coordinates
tune <- autotune(log_price ~ . -price -fold, data = kc_train,
                 k=4, alphaVals = c(0, 0.33, 0.66, 1.0),
                 betaVals = c(0, 0.33), bandwidthVals = c(1.0),
                 powerVals = c(0, 1.0, 1.5, 2.0, 2.5),
                 outputProgress = TRUE)

model1 <- autocart(price_p_sqft ~ . -price -sqft_living -fold, data = kc_train,
                   0, 0)
t1 <- predict(model1, kc_test)
t2 <- predict(model1, kc_test, spatialNodes = TRUE, p = 2)

rmse(t1, kc_test$price_p_sqft)
rmse(t2, kc_test$price_p_sqft)

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
fold_rmse_rf_wo_loc <- rep(NA, length(levels(kc$fold)))
fold_rmse_rf_w_loc <- rep(NA, length(levels(kc$fold)))
fold_rmse_rp_wo_loc <- rep(NA, length(levels(kc$fold)))
fold_rmse_rp_w_loc <- rep(NA, length(levels(kc$fold)))
fold_rmse_ac_wo_loc_nointerp <- rep(NA, length(levels(kc$fold)))
fold_rmse_ac_w_loc_nointerp <- rep(NA, length(levels(kc$fold)))
fold_rmse_ac_wo_loc_fixp <- rep(NA, length(levels(kc$fold)))
fold_rmse_ac_w_loc_fixp <- rep(NA, length(levels(kc$fold)))
fold_rmse_ac_wo_loc_rangep <- rep(NA, length(levels(kc$fold)))
fold_rmse_ac_w_loc_rangep <- rep(NA, length(levels(kc$fold)))
fold_rmse_af_wo_loc_nointerp <- rep(NA, length(levels(kc$fold)))
fold_rmse_af_w_loc_nointerp <- rep(NA, length(levels(kc$fold)))
fold_rmse_af_wo_loc_fixp <- rep(NA, length(levels(kc$fold)))
fold_rmse_af_w_loc_fixp <- rep(NA, length(levels(kc$fold)))
fold_rmse_af_wo_loc_rangep <- rep(NA, length(levels(kc$fold)))
fold_rmse_af_w_loc_rangep <- rep(NA, length(levels(kc$fold)))

fold_size <- rep(NA, length(levels(kc$fold)))

# Use these to record progress of the various optimal values of alpha, beta, etc.
opt_alpha_w_loc <- rep(NA, length(levels(kc$fold)))
opt_alpha_wo_loc <- rep(NA, length(levels(kc$fold)))
opt_beta_w_loc <- rep(NA, length(levels(kc$fold)))
opt_beta_wo_loc <- rep(NA, length(levels(kc$fold)))
opt_bandwidth_w_loc <- rep(NA, length(levels(kc$fold)))
opt_bandwidth_wo_loc <- rep(NA, length(levels(kc$fold)))

num_xval_folds <- length(levels(kc$fold))
for (k in 1:num_xval_folds) {
  n_test <- nrow(kc[kc$fold == k, ])
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
  train_X <- select(kc[kc$fold != k, ], all_of(covariates))
  train_y <- unlist(select(kc[kc$fold != k, ], log_price))
  train_loc <- as.matrix(select(kc[kc$fold != k, ], long, lat))
  test_X <- select(kc[kc$fold == k, ], all_of(covariates))
  test_y <- unlist(select(kc[kc$fold == k, ], log_price))
  test_loc <- as.matrix(select(kc[kc$fold == k, ], long, lat))

  # Create the data frames for the models that use formula specifications
  # (The above format works okay for autocart which uses matrix specifications)
  train_X_wo_loc <- cbind(log_price = train_y, train_X)
  train_X_w_loc <- cbind(train_X_wo_loc, long = train_loc[, 1], lat = train_loc[, 2])
  train_X_w_loc_a <- select(train_X_w_loc, -log_price)
  test_X_wo_loc <- cbind(log_price = test_y, test_X)
  test_X_w_loc <- cbind(test_X_wo_loc, long = test_loc[, 1], lat = test_loc[, 2])
  test_X_w_loc_a <- select(test_X_w_loc, -log_price)

  # Train our models
  rp_wo_loc <- rpart(log_price ~ ., data = train_X_wo_loc)
  rp_w_loc <- rpart(log_price ~ ., data = train_X_w_loc)
  rf_wo_loc <- randomForest(log_price ~ ., data = train_X_wo_loc, ntree = 10)
  rf_w_loc <- randomForest(log_price ~ ., data = train_X_w_loc, ntree = 10)

  # For the thorough testing procedure, set tune_each_cv to TRUE. If you want
  # quick results without tuning hyperparameters then set tune_each_cv = FALSE
  # above at the top of this script.
  if (tune_each_cv) {
    # Cross-validate good parameters for autocart. Note that this needs to be done
    # only with the training data (and not on the withheld fold) because we don't
    # want our hyperparameters to be tuned on something that it should be blind to.
    alpha_vals <- c(0, 0.2, 0.4, 0.6, 0.8)
    beta_vals <- c(0, 0.05)
    bandwidth_vals <- c(0.5, 1)

    # Note that we store separate results for both
    ac_res <- autotune(train_y, train_X, train_loc, k = 3, alphaVals = alpha_vals,
                       betaVals = beta_vals, bandwidthVals = bandwidth_vals,
                       outputProgress = show_cv_progress, useSpatialNodes = TRUE,
                       spatialNodesDistPowerRange = c(4, 5))
    ac_res_w_loc <- autotune(train_y, train_X, train_loc, k = 3, alphaVals = alpha_vals,
                             betaVals = beta_vals, bandwidthVals = bandwidth_vals,
                             outputProgress = show_cv_progress, useSpatialNodes = TRUE,
                             spatialNodesDistPowerRange = c(4, 5))

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
                     numtrees = 50, mtry = 5)
    af_w_loc <- autoforest(train_y, train_X_w_loc_a, train_loc, ac_res_w_loc$alpha, ac_res_w_loc$beta,
                           control = autocartControl(spatialBandwidthProportion = ac_res_w_loc$bandwidth),
                           numtrees = 50, mtry = 5)
  } else {
    # Here we just pick some hyperparameters to cut down on computation time
    # for quicker testing results. This is not the procedure we use for the
    # formal testing as used in the final paper.
    use_alpha <- 0.65
    use_beta <- 0
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
  }

  # Create predictions for the test data
  rf_wo_loc_pred <- predict(rf_wo_loc, test_X_wo_loc)
  rf_w_loc_pred <- predict(rf_w_loc, test_X_w_loc)
  rp_wo_loc_pred <- predict(rp_wo_loc, test_X_wo_loc)
  rp_w_loc_pred <- predict(rp_w_loc, test_X_w_loc)

  ac_wo_loc_pred_nointerp <- predictAutocart(ac, test_X)
  ac_w_loc_pred_nointerp <- predictAutocart(ac_w_loc, test_X_w_loc_a)
  ac_wo_loc_pred_fixp <- spatialNodes(ac, test_X, test_loc, distpower = 4)
  ac_w_loc_pred_fixp <- spatialNodes(ac_w_loc, test_X_w_loc_a, test_loc, distpower = 4)
  ac_wo_loc_pred_rangep <- spatialNodes(ac, test_X, test_loc, distpowerRange = c(4, 4.5))
  ac_w_loc_pred_rangep <- spatialNodes(ac_w_loc, test_X_w_loc_a, test_loc, distpowerRange = c(4, 4.5))

  af_wo_loc_pred_nointerp <- predictAutoforest(af, test_X, test_loc)
  af_w_loc_pred_nointerp <- predictAutoforest(af_w_loc, test_X_w_loc_a, test_loc)
  af_wo_loc_pred_fixp <- predictAutoforest(af, test_X, test_loc,
                                           useSpatialNodes = TRUE,
                                           distpower = 4)
  af_w_loc_pred_fixp <- predictAutoforest(af_w_loc, test_X_w_loc_a,
                                          test_loc, useSpatialNodes = TRUE,
                                          distpower = 4)
  af_wo_loc_pred_rangep <- predictAutoforest(af, test_X, test_loc,
                                               useSpatialNodes = TRUE,
                                               distpowerRange = c(4, 4.5))
  af_w_loc_pred_rangep <- predictAutoforest(af_w_loc, test_X_w_loc_a,
                                              test_loc, useSpatialNodes = TRUE,
                                              distpowerRange = c(4, 4.5))

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
