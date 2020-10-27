prepare_geno <- function(Geno, with_markers) {
  if (with_markers) {
    Geno <- Geno[order(Geno$Line), ]
    rownames(Geno) <- Geno$Line
    Geno <- select(Geno, -Line)
    Geno <- as.matrix(Geno)
  } else {
    lines_names <- Geno$Line
    Geno$Line <- NULL

    colnames(Geno) <- lines_names
    rownames(Geno) <- lines_names

    Geno <- as.matrix(Geno)

    GIDs_Geno <- colnames(Geno)
    GIDs_Geno_Ord <- sort(GIDs_Geno)
    Geno <- Geno[GIDs_Geno_Ord, GIDs_Geno_Ord]
  }
}

prepare_pheno <- function(Geno) {
  Pheno <- Pheno[order(Pheno$Env, Pheno$Line), ]
  rownames(Pheno) <- 1:nrow(Pheno)

  return(Pheno)
}

prepare_mtry <- function(mtry, x_n_cols) {
  if (is.null(mtry)) {
    return(x_n_cols / 3)
  } else {
    return(ceiling(mtry * x_n_cols))
  }
}

prepare_geno_pheno <- function(Pheno, Geno, y, with_markers, is_uni_env, env) {
  Pheno <- prepare_pheno(Pheno)
  Geno <- prepare_geno(Geno, with_markers)

  if (is_uni_env) {
    env_indices <- which(Pheno$Env == env)
    Pheno <- droplevels(Pheno[env_indices, ])
    pheno_lines <- Pheno$Line
  }

  if (with_markers) {
    geno_lines <- rownames(Geno)
  } else {
    geno_lines <- colnames(Geno)
  }

  pheno_lines <- Pheno$Line

  lines_in_both <- unique(geno_lines[which(geno_lines %in% pheno_lines)])
  pheno_indices <- which(Pheno$Line %in% lines_in_both)
  y <- y[pheno_indices]

  Pheno <- Pheno[pheno_indices, ]
  # Reseet factor levels
  Pheno <- droplevels(Pheno)

  geno_indices <- which(geno_lines %in% lines_in_both)
  if (with_markers) {
    Geno <- Geno[geno_indices, ]
  } else {
    Geno <- Geno[geno_indices, geno_indices]
  }
  Geno <- as.matrix(Geno)

  return(list(Pheno = Pheno, Geno = Geno, y = y))
}

rename_env <- function(name) {
  return(gsub( "(as\\.factor\\.Pheno\\.Env\\.)|as\\.factor\\(Pheno\\$Env\\)",
              "", name))
}

rename_interaction <- function(name) {
  new_name <- rename_env(name)

  return(gsub("ZG[0-9]?", "", new_name))
}

rename_line <- function(name) {
  return(name)
}

extract_importances <- function(Importances, forest_model, with_interaction,
                                is_uni_env, n_envs, n_lines, x_cols_names,
                                lines_names) {
  NewImportances <- as.data.frame(t(forest_model$importance))

  betas <- list()

  if (!is_uni_env) {
    envs_indices <- 1:(n_envs)
    betas$Env <- NewImportances[, envs_indices]
    names(betas$Env) <- sapply(names(betas$Env), rename_env)
    Importances$Env <- rbind(Importances$Env, betas$Env)

    lines_indices <- (n_envs + 1):(n_envs + n_lines)
    betas$Line <- NewImportances[, lines_indices]
    colnames(betas$Line) <- lines_names
    Importances$Line <- rbind(Importances$Line, betas$Line)

    if (with_interaction) {
      interaction_indices <- (n_envs + n_lines + 1):ncol(NewImportances)
      betas$GE <- NewImportances[, interaction_indices]
      interaction_names <- x_cols_names[interaction_indices]
      colnames(betas$GE) <- sapply(interaction_names, rename_interaction)
      Importances$GE <- rbind(Importances$GE, betas$GE)
    }
  } else {
    betas$Line <- NewImportances
    colnames(betas$Line) <- lines_names
    Importances$Line <- rbind(Importances$Line, betas$Line)
  }

  return(Importances)
}

train_gen_zap <- function(Pheno, y, Geno,
                          with_interaction,
                          mult_env_anal,

                          ntree_theta,
                          mtry_theta,
                          nodesize_theta,
                          ntree_lambda,
                          mtry_lambda,
                          nodesize_lambda,
                          importance,
                          type,
                          loss_function,

                          cross_validation,
                          number_of_folds,
                          proportion_of_testing,

                          type_of_tuning,
                          tuning_cross_validation,
                          tuning_number_of_folds,
                          tuning_proportion_of_testing,
                          sample_proportion,

                          digits,
                          results_dir,
                          verbose,
                          env,
                          with_markers) {
  is_uni_env <- !is.null(env)
  if (is_uni_env) {
    if (env == "PRN" && .Platform$OS.type == "windows") {
      env <- "PRN_ENV"
    }
    results_dir <- file.path(results_dir, env)
  }

  tmp <- prepare_geno_pheno(Pheno, Geno, y, with_markers, is_uni_env, env)
  Pheno <- tmp$Pheno
  Geno <- tmp$Geno
  y <- tmp$y

  X <- prepare_X(Pheno = Pheno, Geno = Geno, with_markers = with_markers,
                 with_interaction = with_interaction,
                 is_uni_env = is_uni_env)
  mtry_theta <- prepare_mtry(mtry_theta, ncol(X))
  mtry_lambda <- prepare_mtry(mtry_lambda, ncol(X))

  Results <- list(All=data.frame())
  if (importance) {
    Results$ThetaImportances <- list()
    Results$LambdaImportances <- list()
  }

  if (type_of_tuning == "global") {
    if (verbose) {
      cat("**** Tuning ****\n")
    }
    tuning_results <- tune.zap.rfsrc(
      X, y,
      ntree_theta = ntree_theta,
      mtry_theta = mtry_theta,
      nodesize_theta = nodesize_theta,
      ntree_lambda = ntree_lambda,
      mtry_lambda = mtry_lambda,
      nodesize_lambda = nodesize_lambda,
      type = type,
      importance = importance,
      loss_function = loss_function,

      number_of_folds = tuning_number_of_folds,
      cross_validation = tuning_cross_validation,
      proportion_of_testing = tuning_proportion_of_testing,
      sample_proportion = sample_proportion,
      verbose = verbose)

    best_params <- tuning_results$best_params
  }

  folds <- get_folds(cross_validation, number_of_folds, proportion_of_testing,
                     n_records = length(y))
  for (n_current_fold in 1:number_of_folds) {
    if (verbose) {
      cat("Fold", n_current_fold, "\n")
    }
    fold <- folds[[n_current_fold]]
    X_training <- X[fold$training, ]
    y_training <- y[fold$training]
    X_testing <- X[fold$testing, ]
    y_testing <- y[fold$testing]

    if (type_of_tuning == "local") {
      if (verbose) {
        cat("\t**** Tuning ****\n")
      }
      tuning_results <- tune.zap.rfsrc(
        X_training, y_training,
        ntree_theta = ntree_theta,
        mtry_theta = mtry_theta,
        nodesize_theta = nodesize_theta,
        ntree_lambda = ntree_lambda,
        mtry_lambda = mtry_lambda,
        nodesize_lambda = nodesize_lambda,
        type = type,
        importance = importance,
        loss_function = loss_function,

        number_of_folds = tuning_number_of_folds,
        cross_validation = tuning_cross_validation,
        proportion_of_testing = tuning_proportion_of_testing,
        sample_proportion = sample_proportion,
        verbose = verbose)

      best_params <- tuning_results$best_params
    }

    if (type_of_tuning == "global" || nrow(tuning_results$combinations) > 1) {
      fitted_model <- suppressWarnings(zap.rfsrc(
        X_training, y_training,
        importance = importance,
        ntree_theta = best_params$ntree_theta,
        mtry_theta = best_params$mtry_theta,
        nodesize_theta = best_params$nodesize_theta,
        ntree_lambda = best_params$ntree_lambda,
        mtry_lambda = best_params$mtry_lambda,
        nodesize_lambda = best_params$nodesize_lambda))
    } else {
      fitted_model <- tuning_results$zap_fitted_model
    }

    predictions <- predict(fitted_model, X_testing, best_params$type)
    FoldResults <- data.frame(Position = fold$testing,
                              Line = Pheno$Line[fold$testing],
                              Environment = Pheno$Env[fold$testing],
                              Partition = n_current_fold,
                              Observed = y_testing,
                              Predicted = predictions$predicted)
    Results$All <- rbind(Results$All, FoldResults)

    if (importance) {
      n_envs <- length(unique(Pheno$Env))
      lines_names <- rownames(Geno)
      n_lines <- length(lines_names)
      x_cols_names <- colnames(X)

      Results$ThetaImportances <- extract_importances(
        Results$ThetaImportances, fitted_model$theta_forest, with_interaction,
        is_uni_env, n_envs, n_lines, x_cols_names, lines_names)
      Results$LambdaImportances <- extract_importances(
        Results$LambdaImportances, fitted_model$theta_forest, with_interaction,
        is_uni_env, n_envs, n_lines, x_cols_names, lines_names)
    }
  }

  return(Results)
}

cv.zap.rf <- function(Pheno, Geno = NULL, Markers = NULL,
                      with_interaction = TRUE,
                      mult_env_anal = TRUE,

                      ntree_theta = c(500),
                      mtry_theta = NULL,
                      nodesize_theta = c(5),
                      ntree_lambda = c(500),
                      mtry_lambda = NULL,
                      nodesize_lambda = c(5),
                      importance = FALSE,
                      type = "original",
                      loss_function = mse,

                      cross_validation = "k_fold",
                      number_of_folds = 5,
                      proportion_of_testing = 0.2,

                      type_of_tuning = "global",
                      tuning_cross_validation = "k_fold",
                      tuning_number_of_folds = 5,
                      tuning_proportion_of_testing = 0.2,
                      sample_proportion = 1,

                      digits = 4,
                      seed = NULL,
                      results_dir = "gene_zap_random_forest_results",
                      verbose = TRUE) {
  old_state <- NULL
  if (!is.null(seed)) {
    old_state <- get_rand_state()

    set.seed(seed)
  }

  on.exit(set_rand_state(old_state))

  validate_gene_zap_params(
    Pheno, Geno, Markers, with_interaction, mult_env_anal,
    ntree_theta, mtry_theta, nodesize_theta,
    ntree_lambda, mtry_lambda, nodesize_lambda,
    importance, type, loss_function,
    cross_validation, number_of_folds, proportion_of_testing,
    type_of_tuning, tuning_cross_validation, tuning_number_of_folds,
    tuning_proportion_of_testing, sample_proportion,
    results_dir, verbose, digits, seed)

  y <- Pheno$Response

  with_markers <- !is.null(Markers)
  if (is.null(Geno)) {
    Geno <- Markers
    Markers <- NULL
  }

  if (mult_env_anal) {
    Results <- train_gen_zap(
      Pheno, y, Geno, with_interaction, mult_env_anal = TRUE,
      ntree_theta, mtry_theta, nodesize_theta,
      ntree_lambda, mtry_lambda, nodesize_lambda,
      importance, type, loss_function,
      cross_validation, number_of_folds, proportion_of_testing,
      type_of_tuning, tuning_cross_validation,
      tuning_number_of_folds, tuning_proportion_of_testing,
      sample_proportion, digits, results_dir, verbose, env = NULL,
      with_markers)
  } else {
    envs <- as.character(unique(Pheno$Env))
    Results <- list()
    for (env in envs) {
      if (verbose) {
        cat(env, ":\n", sep = "")
      }

      Results[[env]] <- train_gen_zap(
        Pheno, y, Geno, with_interaction, mult_env_anal = FALSE,
        ntree_theta, mtry_theta, nodesize_theta,
        ntree_lambda, mtry_lambda, nodesize_lambda,
        importance, type, loss_function,
        cross_validation, number_of_folds, proportion_of_testing,
        type_of_tuning, tuning_cross_validation,
        tuning_number_of_folds, tuning_proportion_of_testing,
        sample_proportion, digits, results_dir, verbose, env,
        with_markers)
    }
  }

  return(Results)
}