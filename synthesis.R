## synthesis functions

# Function to apply 3rd root transformation to continuous variables
transform_continuous <- function(col) {
  return(sign(col) * abs(col)^(1/3))
}

# Function to retransform continuous variables back to original scale
retransform_continuous <- function(col) {
  return(col^3)
}

discretize_df = function(df, breaks = 5) {
  for (var in colnames(df)) {
    # Check if the variable is not a factor
    if (!is.factor(df[[var]])) {
      
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

synthesize_data_parametric <- function(data, first_var,
                                       seed = 123) {
  set.seed(seed)
  
  # Check that first_var exists
  if (! first_var %in% names(data)) {
    stop("Column '", first_var, "' not found in the data.")
  }
  
  n   <- nrow(data)
  syn <- data.frame(matrix(NA, nrow = n, ncol = ncol(data)),
                    stringsAsFactors = FALSE)
  names(syn) <- names(data)
  
  # Marginally sample first var
  syn[[first_var]] <- sample(data[[first_var]], n, replace = TRUE)
  
  # (OPTIONAL) 3rd‐root transform on continuous cols
  # cont_cols <- names(data)[vapply(data, is.numeric, logical(1))]
  # for (v in cont_cols) data[[v]] <- transform_continuous(data[[v]])
  
  # Seq syn of remaining columns
  for (var_j in setdiff(names(data), first_var)) {
    message("Synthesizing ", var_j)
    filled_preds <- names(syn)[!is.na(syn[1, ])]
    preds        <- setdiff(filled_preds, var_j)
    
    # Prepare training & synthetic pred sets
    train_df     <- data[ , c(preds, var_j), drop = FALSE]
    syn_preds_df <- syn [ , preds,         drop = FALSE]
    
    if (is.factor(data[[var_j]]) || is.character(data[[var_j]])) {
      # categorical with glmnet multinomial
      y_fact <- factor(train_df[[var_j]])
      levs   <- levels(y_fact)
      xmat   <- as.matrix(train_df[ , preds, drop = FALSE])
      
      cvfit  <- cv.glmnet(
        xmat, y_fact,
        family = "multinomial",
        alpha  = 0.5
      )
      
      prob_array <- predict(
        cvfit,
        newx = as.matrix(syn_preds_df),
        type = "response",
        s    = "lambda.min"
      )
      
      if (length(dim(prob_array)) == 3) {
        probs <- prob_array[,,1]
      } else {
        probs <- prob_array
      }
      colnames(probs) <- levs
      
      syn_cat <- apply(probs, 1, function(pr) {
        sample(levs, 1, prob = pr)
      })
      syn[[var_j]] <- factor(syn_cat, levels = levs)
      
    } else {
      # continuous with OLS + Normal residual sampling
      df_fit <- data.frame(
        train_df[ , preds, drop = FALSE],
        y = train_df[[var_j]],
        check.names = FALSE
      )
      lmfit <- lm(y ~ ., data = df_fit)
      
      # get residual standard deviation
      sigma_hat <- summary(lmfit)$sigma
      
      # synthetic predictions
      syn_df   <- data.frame(syn_preds_df, check.names = FALSE)
      syn_pred <- predict(lmfit, newdata = syn_df)
      
      # add i.i.d. normal noise
      syn_cont <- syn_pred + rnorm(n, mean = 0, sd = sigma_hat)
      
      # (OPTIONAL) back‐transform if you used transform_continuous()
      # syn_cont <- retransform_continuous(syn_cont)
      
      syn[[var_j]] <- syn_cont
    }
  }
  
  # Restore factors
  for (v in names(data)) {
    if (is.factor(data[[v]])) {
      syn[[v]] <- factor(syn[[v]], levels = levels(data[[v]]))
    }
  }
  
  return(syn)
}

synthesize_data_cart <- function(data, first_var, seed = 123) {    
  syn_obj <- synthpop::syn(data, visit.sequence = c(first_var, setdiff(colnames(data), first_var)), method = "cart", seed = seed)
  return(syn_obj$syn)
}

synthesize_data_bn <- function(data, seed = 123) {
  # Discretize data
  data <- discretize_df(data)
  
  set.seed(seed)
  
  # learn structure
  bn_structure <- tabu(data)  
  
  # fit with parameters
  bn_fitted <- bn.fit(bn_structure, data, method = "bayes", iss = 1) # to ensure no zero probabilities
  
  # rbn() function generates synthetic data based on the fitted Bayesian network
  syn_data <- rbn(bn_fitted, n = nrow(data))  # gen same number of obs as original dataframe
  
  # ensure factor levels match original dataset
  for (var in colnames(data)) {
    if (is.factor(data[[var]])) {
      syn_data[[var]] <- factor(syn_data[[var]], levels = levels(data[[var]]))
    }
  }
  
  return(syn_data)
}

synthesize_data_svm <- function(data, first_var,
                                C_class      = 1,
                                C_reg        = 1,
                                epsilon      = 0.1,
                                k_neighbors  = 10,
                                seed         = 123) {
  set.seed(seed)
  
  # Copy original data
  orig <- data.frame(data, stringsAsFactors = FALSE)
  
  # Initialize synthetic data frame
  syn_data <- data.frame(matrix(NA, nrow = nrow(orig), ncol = ncol(orig)))
  names(syn_data) <- names(orig)
  
  # Sample first_var 
  syn_data[[first_var]] <- sample(orig[[first_var]], nrow(orig), replace = TRUE)
  
  # Seq syn
  for (var_j in setdiff(names(orig), first_var)) {
    
    # predictors
    preds_filled <- names(syn_data)[!is.na(syn_data[1, ])]
    predictors   <- setdiff(preds_filled, var_j)
    
    # training data
    train_df <- orig[, c(predictors, var_j), drop = FALSE]
    
    # synpreds for var_j
    syn_preds_df <- syn_data[, predictors, drop = FALSE]
    if (is.factor(orig[[var_j]]) || is.character(orig[[var_j]])) {
      
      # Classification
      train_df[[var_j]] <- factor(train_df[[var_j]])
      levs <- levels(train_df[[var_j]])
      model_cl <- svm(
        as.formula(paste(var_j, "~ .")),
        data        = train_df,
        probability = TRUE,
        kernel      = "radial",
        cost        = C_class
      )
      pred_obj <- predict(model_cl, newdata = syn_preds_df, probability = TRUE)
      prob_mat <- attr(pred_obj, "probabilities")
      
      # sample each row according to its probability vector
      syn_data[[var_j]] <- apply(prob_mat, 1, function(p_row) {
        sample(levs, size = 1, prob = p_row)
      })
      syn_data[[var_j]] <- factor(syn_data[[var_j]], levels = levs)
    } else if (is.numeric(orig[[var_j]])) {
      
      # Regression
      model_rg <- svm(
        as.formula(paste(var_j, "~ .")),
        data    = train_df,
        type    = "eps-regression",
        kernel  = "radial",
        cost    = C_reg,
        epsilon = epsilon
      )
      
      # train predictions & residuals
      preds_train <- predict(model_rg, newdata = train_df)
      resid_train <- train_df[[var_j]] - preds_train
      
      # synthetic predictions
      syn_preds <- predict(model_rg, newdata = syn_preds_df)
      
      # Find kNN among train preds for each syn_pred
      nn <- get.knnx(
        data  = matrix(preds_train, ncol = 1),
        query = matrix(syn_preds,   ncol = 1),
        k     = k_neighbors
      )$nn.index  
      
      # sample one residual each
      sampled_resid <- apply(nn, 1, function(ix) {
        sample(resid_train[ix], 1)
      })
      syn_data[[var_j]] <- syn_preds + sampled_resid
    } else {
      stop("Variable ", var_j, " is neither factor nor numeric.")
    }
  }
  return(syn_data)
}

# run synthesis
synthesize_data_model <- function(simulated_data, synth_method, mc_reps = 5) {
  
  purrr::imap(simulated_data, ~ {
    
    cat("Synthesizing dataset:", .y, "with", synth_method$name, "\n")
    df <- as.data.frame(.x)
    
    # Identify categorical and continuous variables
    categorical_vars <- grep("^C", names(df), value = TRUE)
    continuous_vars <- grep("^X", names(df), value = TRUE)
    
    if (length(categorical_vars) > 0) {
      
      # Compute balance and cardinality for categorical variables
      cat_scores <- sapply(categorical_vars, function(var) {
        tab <- table(df[[var]])
        props <- tab / sum(tab)
        k <- length(props)
        gini_raw <- 1 - sum(props^2)
        gini_norm <- gini_raw / (1 - 1 / k)
        balance <- gini_norm  # Gini-based balance score in [0, 1]
        
        cardinality <- length(tab)
        max_cardinality <- max(sapply(categorical_vars, function(v) length(table(df[[v]]))), na.rm = TRUE)
        
        score <- 0.8 * balance + 0.2 * (1 - (cardinality / max_cardinality))
        return(score)
      })
      
      first_var <- names(which.max(cat_scores))
      
      # Order categorical variables by increasing cardinality, excluding first_var
      remaining_categorical_vars <- setdiff(categorical_vars, first_var)
      remaining_categorical_vars <- remaining_categorical_vars[
        order(sapply(remaining_categorical_vars, function(var) length(unique(df[[var]]))))
      ]
      
      if (length(continuous_vars) > 0) {
        ordered_vars <- c(first_var, continuous_vars, remaining_categorical_vars)
      } else {
        ordered_vars <- c(first_var, remaining_categorical_vars)  # No continuous variables
      }
      
    } else if (length(continuous_vars) > 0) {
      first_var <- continuous_vars[1]
      ordered_vars <- continuous_vars  # Keep order as given, no categorical variables
      
    } else {
      warning("Dataset", .y, "contains no variables. Skipping synthesis.")
      return(NULL)
    }
    
    cat("  Using first_var:", first_var, "\n")
    
    # Reorder dataset
    df <- df[ordered_vars]
    
    # Perform synthesis with the given method
    synthetic_results <- purrr::map(1:mc_reps, function(mc_iter) {
      set.seed(mc_iter)
      cat("  Applying", synth_method$name, "with seed", mc_iter, "on", .y, "\n")
      
      # Function arguments
      args <- list(data = df)
      if ("first_var" %in% names(formals(synth_method$func))) {
        args$first_var <- first_var
      }
      if ("seed" %in% names(formals(synth_method$func))) {
        args$seed <- mc_iter
      }
      
      # Synthesize data
      do.call(synth_method$func, args)
    })
    
    # Flatten list, name elements
    names(synthetic_results) <- paste0(.y, "_synmc", 1:mc_reps)
    return(synthetic_results)
  }) %>% purrr::flatten()
}