## simulation functions

generate_data_linear <- function(n = 1000,
                                 p = 10,
                                 #rho = 0.8,
                                 min_obs = 20,
                                 max_levels = 7,
                                 seed = 123) {
  set.seed(seed)
  
  # continuous vs. categorical
  p_cont <- floor(p/2)
  p_cat  <- p - p_cont
  
  # latent dimension
  d <- max(max_levels, p_cont)
  
  # build a positive‐definite covariance
  cov_mat <- matrix(runif(d * d, min = 0.05, max = 0.95), nrow = d)
  Sigma   <- t(cov_mat) %*% cov_mat
  mu      <- rep(0, d)
  
  # draw latent MVN once
  base_mvn <- MASS::mvrnorm(n, mu, Sigma)
  
  # continuous 
  if (p_cont > 0) {
    cont_dims <- sample.int(d, p_cont)
    cont_df   <- as.data.frame(base_mvn[, cont_dims, drop = FALSE])
    names(cont_df) <- paste0("X", seq_len(p_cont))
  } else {
    cont_df <- NULL
  }
  
  # categorical
  if (p_cat > 0) {
    cat_levels <- sample(2:max_levels, p_cat, replace = TRUE)
    cat_df <- data.frame(matrix(nrow = n, ncol = p_cat))
    for (i in seq_len(p_cat)) {
      num_levels    <- cat_levels[i]
      selected_dims <- sample.int(d, num_levels - 1)
      X_pred        <- base_mvn[, selected_dims, drop = FALSE]
      CovRes        <- diag(num_levels)
      resid         <- MASS::mvrnorm(n, rep(0, num_levels), CovRes)
      betas         <- matrix(
        runif(length(selected_dims) * (num_levels - 1),
              min = -3, max = 3),
        nrow = length(selected_dims),
        ncol = num_levels - 1
      )
      betas         <- cbind(0, betas)  # intercept = 0 for first level
      ylat          <- X_pred %*% betas + resid
      cat_df[[i]]   <- max.col(ylat)
      
      # enforce min_obs per level
      counts <- table(factor(cat_df[[i]], levels = 1:num_levels))
      bad    <- as.integer(names(counts[counts < min_obs]))
      good   <- setdiff(seq_len(num_levels), bad)
      
      if (length(good) < 2) {
        stop("Too few valid categories in variable ", i)
      }
      
      idx_bad <- which(cat_df[[i]] %in% bad)
      if (length(idx_bad) > 0) {
        ylat_sub <- ylat[idx_bad, good, drop = FALSE]
        cat_df[[i]][idx_bad] <- good[max.col(ylat_sub)]
      }
      
      # re‐label
      obs_lvls   <- sort(unique(cat_df[[i]]))
      relabeled  <- match(cat_df[[i]], obs_lvls)
      cat_df[[i]] <- factor(relabeled, levels = seq_along(obs_lvls))
    }
    
    names(cat_df) <- paste0("C", seq_len(p_cat))
    
  } else {
    cat_df <- NULL
  }
  
  # combine and return
  data.frame(cont_df, cat_df, check.names = FALSE)
}

# Helper: sample a truncated normal quantile
sample_quantile <- function(mu = 0.5, sigma = 0.2, lower = 0.1, upper = 0.9) {
  q <- rnorm(1, mean = mu, sd = sigma)
  pmin(pmax(q, lower), upper)
}

# Recursive tree builder: returns a nested list of nodes
build_tree <- function(indices, df_prev, depth = 1,
                       max_depth = 5, min_split = 200, min_bucket = 50) {
  recurse <- function(idx, current_depth) {
    node <- list(indices = idx)
    n <- length(idx)
    # stop if too small or depth limit
    if (n < min_split || current_depth > max_depth) {
      node$is_leaf <- TRUE
      return(node)
    }
    # try up to 5 splits
    for (try in 1:5) {
      split_var <- sample(names(df_prev), 1)
      values    <- df_prev[idx, split_var]
      # if numeric but no variation, make leaf
      if (is.numeric(values) && diff(range(values)) == 0) {
        node$is_leaf <- TRUE
        return(node)
      }
      # if categorical but only one level present, make leaf
      if (!is.numeric(values) && length(unique(values)) < 2) {
        node$is_leaf <- TRUE
        return(node)
      }
      if (is.numeric(values)) {
        q      <- sample_quantile()
        thresh <- quantile(values, q, names = FALSE)
        left_idx  <- idx[values <= thresh]
        right_idx <- idx[values >  thresh]
        split_value <- thresh
      } else {
        levs           <- unique(values)
        k              <- length(levs)
        pick_n         <- sample(1:(k-1), 1)
        subset_levels  <- sample(levs, pick_n)
        left_idx  <- idx[values %in% subset_levels]
        right_idx <- idx[!(values %in% subset_levels)]
        split_value <- subset_levels
      }
      if (length(left_idx) >= min_bucket && length(right_idx) >= min_bucket) {
        node$is_leaf     <- FALSE
        node$split_var   <- split_var
        node$split_value <- split_value
        node$left        <- recurse(left_idx,  current_depth + 1)
        node$right       <- recurse(right_idx, current_depth + 1)
        return(node)
      }
      message(sprintf(
        "  Split attempt %d at depth %d on %s failed (sizes: %d, %d)",
        try, current_depth, split_var, length(left_idx), length(right_idx)
      ))
    }
    warning(sprintf("Could not split node at depth %d; making leaf", current_depth))
    node$is_leaf <- TRUE
    return(node)
  }
  recurse(indices, depth)
}

# Generate continuous variable based on tree leaves
generate_continuous <- function(tree, df_prev, a, b) {
  N <- nrow(df_prev)
  X <- numeric(N)
  assign_leaf <- function(node, interval) {
    if (node$is_leaf) {
      idx   <- node$indices
      l     <- interval[1]; u <- interval[2]
      mu    <- (l + u) / 2
      sigma <- (u - l) / 8
      X[idx] <<- rnorm(length(idx), mean = mu, sd = sigma)
    } else {
      assign_leaf(node$left,  interval)
      assign_leaf(node$right, interval)
    }
  }
  assign_leaf(tree, c(a, b))
  X
}

# Generate categorical variable based on tree leaves
generate_categorical <- function(tree, df_prev) {
  N <- nrow(df_prev)
  # collect leaves
  leaves <- list()
  collect <- function(node) {
    if (node$is_leaf) leaves[[length(leaves) + 1]] <<- node
    else {
      collect(node$left)
      collect(node$right)
    }
  }
  collect(tree)
  L <- length(leaves)
  K <- sample(2:7, 1)
  
  # ensure each category assigned to at least one leaf
  leaf_cats <- integer(L)
  ord       <- sample(L)
  leaf_cats[ord[1:K]]       <- seq_len(K)
  leaf_cats[ord[(K+1):L]]   <- sample(seq_len(K), L - K, replace = TRUE)
  
  # label each observation
  C <- integer(N)
  for (i in seq_along(leaves)) {
    idx        <- leaves[[i]]$indices
    assigned_c <- leaf_cats[i]
    for (j in idx) {
      C[j] <- if (runif(1) < 0.8) assigned_c else sample((1:K)[-assigned_c], 1)
    }
  }
  if (any(table(C) < 20)) {
    stop("Could not assign categories with ≥20 observations each after 5 tries")
  }
  factor(C)  # return as factor
}

# Main simulation function producing numeric Xs and factor Cs
generate_data_hierarchical <- function(
    N = 10000,
    a = 0, b = 10,
    p = 5,
    seed = 123,
    return_trees = FALSE
) {
  set.seed(seed)
  df    <- data.frame(id = seq_len(N))
  trees <- list()
  
  for (j in seq_len(p)) {
    var_name <- if (j %% 2 == 1) paste0("X", (j + 1) %/% 2) else paste0("C", j %/% 2)
    is_cont  <- (j %% 2 == 1)
    message(sprintf("Generating variable %s (%s)…", 
                    var_name, ifelse(is_cont, "continuous", "categorical")))
    if (j == 1 && is_cont) {
      df[[var_name]] <- runif(N, min = a, max = b)
      if (return_trees) trees[[var_name]] <- NULL
    } else {
      df_prev     <- df[, setdiff(names(df), "id"), drop = FALSE]
      tree        <- build_tree(seq_len(N), df_prev)
      if (is_cont) {
        df[[var_name]] <- generate_continuous(tree, df_prev, a, b)
      } else {
        df[[var_name]] <- generate_categorical(tree, df_prev)
      }
      if (return_trees) trees[[var_name]] <- tree
    }
  }
  
  result <- list(data = df[, -1])
  if (return_trees) result$trees <- trees
  invisible(result)
}

