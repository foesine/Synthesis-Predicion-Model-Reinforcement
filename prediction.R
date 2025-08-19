## pred functions

## Calculate evaluation metrics for continuous targets
evaluation_metrics_cont <- function(predictions, true_values) {
  residuals <- predictions - true_values
  MAE       <- mean(abs(residuals))
  RMSE      <- sqrt(mean(residuals^2))
  ss_res    <- sum(residuals^2)
  ss_tot    <- sum((true_values - mean(true_values))^2)
  R2        <- if (ss_tot == 0) NA else 1 - ss_res/ss_tot
  data.frame(MAE = MAE, RMSE = RMSE, R_squared = R2)
}


## Calculate evaluation metrics for factored targets with error handling
evaluation_metrics_factor <- function(predictions, true_values) {
  # sanity check lengths
  if (length(predictions) != length(true_values)) {
    stop("Length mismatch: ", length(predictions),
         " predictions vs. ", length(true_values), " true values.")
  }
  true_values <- factor(true_values)
  predictions <- factor(predictions, levels = levels(true_values))
  if (nlevels(true_values) < 2) {
    stop("The target must have at least two levels.")
  }
  cm <- confusionMatrix(predictions, true_values, mode = "everything")
  acc     <- cm$overall["Accuracy"]
  support <- table(true_values)
  total   <- sum(support)
  if (nlevels(true_values) == 2) {
    sens <- cm$byClass["Sensitivity"]
    spec <- cm$byClass["Specificity"]
  } else {
    sens <- sum(cm$byClass[,"Sensitivity"] * support) / total
    spec <- sum(cm$byClass[,"Specificity"] * support) / total
  }
  data.frame(Accuracy = acc, Sensitivity = sens, Specificity = spec)
}
# evaluation_metrics_factor <- function(predictions, true_values) {
#   if (length(predictions) != length(true_values)) {
#     stop("Length mismatch: ", length(predictions),
#          " predictions vs. ", length(true_values), " true values.")
#   }
#   
#   true_values <- factor(true_values)
#   predictions <- factor(predictions, levels = levels(true_values))
#   if (nlevels(true_values) < 2) {
#     stop("The target must have at least two levels.")
#   }
#   
#   cm  <- caret::confusionMatrix(predictions, true_values, mode = "everything")
#   acc <- unname(cm$overall["Accuracy"])
#   byc <- cm$byClass
#   support <- table(true_values)
#   total   <- sum(support)
#   
#   if (is.null(dim(byc))) {
#     ## Binary: byClass is a named vector
#     sens <- unname(byc["Sensitivity"])
#     spec <- unname(byc["Specificity"])
#     f1   <- unname(byc["F1"])
#   } else {
#     # Align rows to the levels order in true_values
#     row_lvls <- sub("^Class: ", "", rownames(byc))
#     ord      <- match(levels(true_values), row_lvls)
#     
#     sens_vec <- byc[ord, "Sensitivity", drop = TRUE]
#     spec_vec <- byc[ord, "Specificity", drop = TRUE]
#     f1_vec   <- byc[ord, "F1",          drop = TRUE]
#     
#     # support aligned to levels(true_values)
#     supp_vec <- as.numeric(support[levels(true_values)])
#     
#     # support-weighted averages
#     sens <- sum(sens_vec * supp_vec, na.rm = TRUE) / total
#     spec <- sum(spec_vec * supp_vec, na.rm = TRUE) / total
#     f1   <- sum(f1_vec   * supp_vec, na.rm = TRUE) / total
#   }
#   
#   data.frame(
#     Accuracy    = as.numeric(acc),
#     Sensitivity = as.numeric(sens),
#     Specificity = as.numeric(spec),
#     F1          = as.numeric(f1),
#     row.names   = NULL
#   )
# }

discretize_df = function(df, breaks = 5) {
  for (var in colnames(df)) {
    # Check if the variable is not a factor
    if (is.numeric(df[[var]])) {
      
      # Count the frequency of each unique value
      freq_table <- table(df[[var]])
      
      # Calculate the proportion of zeros, ensuring NA is handled
      zero_proportion <- ifelse(!is.na(freq_table[as.character(0)]),
                                freq_table[as.character(0)] / sum(freq_table),
                                0)
      
      # Determine the number of breaks based on zero proportion
      if (zero_proportion > 4/5) {
        new_breaks = 1
      } else if (zero_proportion > 1/4) {
        new_breaks = breaks - 2
      } else if (zero_proportion > 1/5) {
        new_breaks = breaks - 1
      } else {
        new_breaks = breaks
      }
      
      # Separate zeros and non-zeros
      zero_portion = (df[[var]] == 0)
      non_zero_values = df[[var]][!zero_portion]
      
      # Discretize non-zero values
      if (length(non_zero_values) > 0) {
        # Calculate breaks for non-zero values
        range_values = range(non_zero_values, na.rm = TRUE)
        breaks_values = seq(range_values[1], range_values[2], length.out = new_breaks + 1)
        
        # Ensure correct number of labels are created
        labels = sapply(1:(length(breaks_values)-1), function(i)
          paste("(", breaks_values[i], "-", breaks_values[i+1], "]", sep=""))
        
        
        # Use cut to apply these breaks and labels
        discretized_non_zeros = cut(non_zero_values, breaks = breaks_values, labels = labels, include.lowest = TRUE)
        # Combine zero and discretized non-zeros into the original dataframe
        df[[var]] <- factor(ifelse(zero_portion, "0", as.character(discretized_non_zeros)))
      } else {
        # If all values are zero or the number of breaks is zero or negative
        df[[var]] <- factor("0")
      }
    }
  }
  return(df)
}

parm_pred <- function(data, target_var,
                      outer_folds = 5,
                      alpha        = 0.5,
                      inner_folds  = 3,
                      seed         = 123) {
  
  set.seed(seed)
  if (!target_var %in% names(data)) {
    stop("Target variable '", target_var, "' not found.")
  }
  df <- as.data.frame(data)
  
  # extract and prepare y
  y0 <- df[[target_var]]
  if (is.character(y0))    y0 <- factor(y0)
  if (is.factor(y0)) {
    task <- "classification"
    levels(y0) <- make.names(levels(y0))
    df[[target_var]] <- y0
  } else if (is.numeric(y0) && length(unique(y0)) == 2) {
    task <- "classification"
    df[[target_var]] <- y0 <- factor(y0)
  } else if (is.numeric(y0)) {
    task <- "regression"
  } else {
    stop("Target must be numeric or factor.")
  }
  
  family <- if (task == "regression") "gaussian" else {
    if (nlevels(y0) == 2) "binomial" else "multinomial"
  }
  
  # stratify on training folds to preserve all classes
  train_idx_list <- createFolds(df[[target_var]], 
                                k           = outer_folds,
                                returnTrain = TRUE)
  
  outer_results <- lapply(train_idx_list, function(train_idx) {
    # train/test
    test_idx  <- setdiff(seq_len(nrow(df)), train_idx)
    train_set <- df[ train_idx, ]
    test_set  <- df[ test_idx, ]
    
    # drop unused levels
    if (task == "classification") {
      train_set[[target_var]] <- droplevels(train_set[[target_var]])
      test_set[[target_var]]  <- factor(test_set[[target_var]],
                                        levels = levels(train_set[[target_var]]))
    }
    
    # design matrices
    X_train <- model.matrix(~ . - 1,
                            data = train_set[, setdiff(names(train_set), target_var), drop = FALSE])
    y_train <- train_set[[target_var]]
    X_test  <- model.matrix(~ . - 1,
                            data = test_set[,  setdiff(names(test_set),  target_var), drop = FALSE])
    y_test  <- test_set[[target_var]]
    
    # inner CV for lambda
    cvfit <- cv.glmnet(
      x                = X_train,
      y                = y_train,
      family           = family,
      alpha            = alpha,
      nfolds           = inner_folds,
      standardize      = TRUE,
      lambda.min.ratio = 1e-2,
      nlambda          = 10,
      type.multinomial = if (family == "multinomial") "grouped" else NULL
    )
    best_lambda <- cvfit$lambda.min
    
    # final fit
    final_mod <- glmnet(
      x                 = X_train,
      y                 = y_train,
      family            = family,
      alpha             = alpha,
      lambda            = best_lambda,
      standardize       = TRUE,
      type.multinomial  = if (family == "multinomial") "grouped" else NULL
    )
    
    # predict & evaluate
    if (task == "classification") {
      preds <- predict(final_mod,
                       newx = X_test,
                       s    = best_lambda,
                       type = "class")
      evaluation_metrics_factor(preds, y_test)
    } else {
      preds <- as.numeric(predict(final_mod,
                                  newx = X_test,
                                  s    = best_lambda))
      evaluation_metrics_cont(preds, y_test)
    }
  })
  
  # aggregate
  df_out <- bind_rows(outer_results)
  df_out %>% summarise(across(everything(), mean, na.rm = TRUE))
}


cart_pred <- function(data, target_var,
                      outer_folds = 5, cp_steps = 10,
                      inner_folds = 3, seed = 123) {
  set.seed(seed)
  if (!target_var %in% names(data)) stop("Target not in data.")
  # ensure factor/character targets are factors
  if (is.character(data[[target_var]])) {
    data[[target_var]] <- factor(data[[target_var]])
  }
  # choose task
  if (is.factor(data[[target_var]])) {
    data[[target_var]] <- make.names(data[[target_var]])
    data <- data %>% mutate(across(where(is.character), as.factor))
    task <- "class"
    summaryFunc <- defaultSummary
  } else if (is.numeric(data[[target_var]])) {
    task <- "reg"
    summaryFunc <- defaultSummary
  } else {
    stop("Target must be numeric or factor.")
  }
  
  outer_ctrl <- trainControl(method = "cv", number = outer_folds,
                             summaryFunction = summaryFunc,
                             verboseIter = FALSE)
  inner_ctrl <- trainControl(method = "cv", number = inner_folds,
                             summaryFunction = summaryFunc,
                             verboseIter = FALSE)
  
  cp_vals  <- 10^seq(log10(1e-4), log10(1e-2), length.out = cp_steps)
  tunegrid <- expand.grid(cp = cp_vals)
  
  if (task == "class") {
    outer_idx <- createFolds(data[[target_var]], k = outer_folds)
  } else {
    outer_idx <- createFolds(seq_len(nrow(data)), k = outer_folds)
  }
  
  results <- lapply(outer_idx, function(test_idx) {
    train_data <- data[-test_idx, ]; test_data <- data[test_idx, ]
    
    # inner CV for best cp
    if (task == "class") {
      inner_idx <- createFolds(train_data[[target_var]], k = inner_folds)
    } else {
      inner_idx <- createFolds(seq_len(nrow(train_data)), k = inner_folds)
    }
    best_cps <- sapply(inner_idx, function(idx_i) {
      dti <- train_data[-idx_i, ]; vti <- train_data[idx_i, ]
      m <- train(as.formula(paste(target_var, "~ .")), data = dti,
                 method    = "rpart", tuneGrid = tunegrid,
                 trControl = inner_ctrl,
                 control   = rpart.control(maxsurrogate = 0, maxcompete = 1))
      m$bestTune$cp
    })
    best_cp <- as.numeric(names(sort(table(best_cps), decreasing = TRUE))[1])
    
    final_m <- train(as.formula(paste(target_var, "~ .")), data = train_data,
                     method    = "rpart",
                     tuneGrid  = data.frame(cp = best_cp),
                     trControl = outer_ctrl,
                     control   = rpart.control(maxsurrogate = 0, maxcompete = 1))
    preds <- predict(final_m, newdata = test_data)
    
    if (task == "class") {
      evaluation_metrics_factor(preds, test_data[[target_var]])
    } else {
      evaluation_metrics_cont(preds, test_data[[target_var]])
    }
  })
  
  df <- do.call(rbind, results)
  df %>% summarise(across(everything(), mean, na.rm = TRUE))
}


bn_pred <- function(data, target_var, outer_folds = 5, inner_folds = 3, seed = 123) {
  set.seed(seed)
  if (!(target_var %in% colnames(data))) {
    stop("Target variable '", target_var, "' not found in the dataset.")
  }
  
  data <- discretize_df(data)
  
  # Convert the target variable to a factor
  # Convert all character columns to factors
  data <- data %>%
    dplyr::mutate(across(where(is.character), as.factor))
  
  # Make sure target is a factor with proper levels
  data[[target_var]] <- factor(data[[target_var]], levels = unique(data[[target_var]]))
  
  algorithms <- c("tabu")
  
  outer_cv_folds <- createFolds(data[[target_var]], k = outer_folds)
  outer_results <- list()
  
  for (i in seq_along(outer_cv_folds)) {
    outer_test_index <- outer_cv_folds[[i]]
    outer_testData <- data[outer_test_index, ]
    outer_trainData <- data[-outer_test_index, ]
    
    if (length(unique(outer_trainData[[target_var]])) < 2) {
      warning("Outer Fold ", i, ": The target variable has less than two levels in the training set.")
      next
    }
    
    # inner CV
    inner_folds_indices <- createFolds(outer_trainData[[target_var]], k = inner_folds)
    best_model <- NULL
    best_performance <- -Inf
    best_algorithm <- NULL
    
    for (algorithm in algorithms) {
      #cat("Trying algorithm:", algorithm, "\n")
      fold_results <- c()
      
      for (j in seq_along(inner_folds_indices)) {
        
        inner_test_index <- inner_folds_indices[[j]]
        inner_trainData <- outer_trainData[-inner_test_index, ]
        inner_testData <- outer_trainData[inner_test_index, ]
        
        if (length(unique(inner_trainData[[target_var]])) < 2) {
          warning("Inner Fold ", j, ": The target variable has less than two levels in the training set.")
          fold_results[j] <- NA
          next
        }
        
        # Build Bayesian Network model
        bn_model <- do.call(get(algorithm, envir = asNamespace("bnlearn")), list(inner_trainData))
        fitted_bn_model <- bnlearn::bn.fit(bn_model, inner_trainData)
        
        # Predictions
        predictions <- predict(fitted_bn_model, node = target_var, data = inner_testData, method = "bayes-lw")
        predictions <- factor(predictions, levels = levels(inner_trainData[[target_var]]))
        
        # Calculate accuracy
        accuracy <- mean(predictions == inner_testData[[target_var]], na.rm = TRUE)
        fold_results[j] <- accuracy
      }
      
      avg_performance <- mean(fold_results, na.rm = TRUE)
      
      if (!is.na(avg_performance) && avg_performance > best_performance) {
        best_performance <- avg_performance
        best_model <- fitted_bn_model
        best_algorithm <- algorithm
      }
    }
    #cat("Best algorithm selected in outer fold", i, ":", best_algorithm, "with accuracy:", best_performance, "\n")
    
    # Apply best model to outer test set
    if (!is.null(best_model) && length(unique(outer_testData[[target_var]])) >= 2) {
      predictions <- predict(best_model, node = target_var, data = outer_testData, method = "bayes-lw")
      predictions <- factor(predictions, levels = levels(outer_testData[[target_var]]))
      
      eval <- evaluation_metrics_factor(predictions = predictions, 
                                        true_values = outer_testData[[target_var]])
      outer_results[[i]] <- eval
      
    } else {
      cat("Skipping evaluation in outer fold", i, "due to insufficient class diversity.\n")
    }
  }
  
  # Average evaluation metrics over outer folds
  if (length(outer_results) > 0) {
    eval_avg_outer_folds <- do.call(rbind, outer_results) %>%
      dplyr::summarise(across(everything(), mean, na.rm = TRUE))
    
    return(eval_avg_outer_folds)
    
  } else {
    
    warning("No valid folds to evaluate.")
    return(NULL)
  }
}

svm_pred <- function(data, target_var,
                     outer_folds = 5,
                     cost_steps   = 10,
                     inner_folds  = 3,
                     seed         = 123) {
  
  set.seed(seed)
  if (!target_var %in% names(data)) {
    stop("Target variable ", target_var, " not found in the dataset")
  }
  if (length(data[[target_var]]) == 0) {
    stop("Target variable is empty. Please provide a valid dataset.")
  }
  # Extract y
  y <- data[[target_var]]
  # Determine task type
  if (is.character(y)) {
    y <- factor(y)
  }
  if (is.factor(y)) {
    task <- "classification"
  } else if (is.numeric(y)) {
    task <- "regression"
  } else {
    stop("Target must be numeric or factor.")
  }
  # Choose summary function
  summaryFunctionType <- if (task == "regression") defaultSummary else multiClassSummary
  outer_control <- trainControl(
    method          = "cv",
    number          = outer_folds,
    summaryFunction = summaryFunctionType,
    verboseIter     = FALSE,
    allowParallel   = FALSE
  )
  inner_control <- trainControl(
    method          = "cv",
    number          = inner_folds,
    summaryFunction = summaryFunctionType,
    verboseIter     = FALSE,
    allowParallel   = FALSE
  )
  
  cost_values <- 10^seq(log10(0.001), log10(100), length.out = cost_steps)
  tunegrid    <- expand.grid(C = cost_values) #, sigma = 0.1)
  
  # Outer folds
  outer_cv_folds <- createFolds(data[[target_var]], k = outer_folds)
  all_results <- lapply(outer_cv_folds, function(test_idx) {
    train_set <- data[-test_idx, ]
    test_set  <- data[ test_idx, ]
    
    # Inner tuning
    ctrl <- inner_control
    model <- train(
      formula(paste(target_var, "~ .")),
      data     = train_set,
      method   = "svmRadialCost",
      tuneGrid = tunegrid,
      trControl= ctrl
    )
    best_params <- model$bestTune
    
    # Final fit on outer-train
    final_model <- train(
      formula(paste(target_var, "~ .")),
      data     = train_set,
      method   = "svmRadialCost",
      tuneGrid = best_params,
      trControl= outer_control
    )
    
    # Predictions
    preds <- predict(final_model, newdata = test_set)
    
    # Evaluate
    if (task == "regression") {
      evaluation_metrics_cont(preds, test_set[[target_var]])
    } else {
      # classification
      evaluation_metrics_factor(preds, test_set[[target_var]])
    }
  })
  
  # Average across outer folds
  df_out <- bind_rows(all_results)
  df_out %>% summarise(across(everything(), mean, na.rm = TRUE))
}


# Run prediction
run_predictions <- function(data_list,
                            prediction_function,
                            seed = 123,
                            # default target sets for each p
                            targets_p6  = c("X1","X2","X3","C1","C2","C3"),
                            targets_p12 = c("X2","X4","X6","C2","C4","C6"),
                            targets_p18 = c("X3","X6","X9","C3","C6","C9")) {
  set.seed(seed)
  
  predictions <- imap(data_list, function(dataset, dataset_name) {
    # Parse p from names like "p6_mc1", "p12_mc3", etc.
    p_val <- as.integer(sub("^p(\\d+)_.*$", "\\1", dataset_name))
    # Choose the appropriate default targets
    if (p_val == 6) {
      selected_targets <- targets_p6
    } else if (p_val == 12) {
      selected_targets <- targets_p12
    } else if (p_val == 18) {
      selected_targets <- targets_p18
    } else {
      stop("Unexpected p = ", p_val, " in dataset name ", dataset_name)
    }
    # Split into continuous vs categorical for messaging
    cont_targets <- grep("^X", selected_targets, value = TRUE)
    cat_targets  <- grep("^C", selected_targets, value = TRUE)
    cat("\nDataset:", dataset_name, "\n",
        "  Continuous targets:", paste(cont_targets, collapse = ", "), "\n",
        "  Categorical targets:", paste(cat_targets, collapse = ", "), "\n")
    # Run prediction_function for each target
    res <- map(selected_targets, function(tv) {
      cat("    -> Predicting with target:", tv, "\n")
      if (!tv %in% names(dataset)) {
        warning("      Skipping missing target: ", tv)
        return(NULL)
      }
      prediction_function(
        data        = dataset,
        target_var  = tv,
        outer_folds = 5,
        inner_folds = 3,
        seed        = seed
      )
    })
    set_names(res, selected_targets)
  })
  
  save(predictions, file = "predictions_sim_data.RData")
  predictions
}