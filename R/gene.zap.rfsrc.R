prepare_geno <- function(Geno, with_markers) {
  if (with_markers) {
    Geno <<- Geno[order(Geno$Line), ]
    rownames(Geno) <<- Geno$Line
    Geno <<- select(Geno, -Line)
    Geno <<- as.matrix(Geno)
  } else {
    lines_names <- Geno$Line

    colnames(Geno) <<- lines_names
    rownames(Geno) <<- lines_names

    Geno <<- as.matrix(Geno)

    GIDs_Geno <- colnames(Geno)
    GIDs_Geno_Ord <- sort(GIDs_Geno)
    Geno <<- Geno[GIDs_Geno_Ord, GIDs_Geno_Ord]
  }
}

prepare_pheno <- function(Geno) {
  Pheno <<- Pheno[order(Pheno$Env, Pheno$Line), ]
  rownames(Pheno) <<- 1:nrow(Pheno)
}

prepare_mtry <- function(mtry, x_n_cols) {
  if (is.null(mtry)) {
    return(x_n_cols / 3)
  } else {
    return(ceiling(mtry * x_n_cols))
  }
}

prepare_geno_pheno <- function(Pheno, Geno, is_uni_env, env) {
  Pheno <- prepare_pheno(Pheno)
  Geno <- prepare_geno(Geno)

  if (is_uni_env) {
    env_indices <- which(Pheno$Env == env)
    Pheno <<- droplevels(Pheno[env_indices, ])
    pheno_lines <- Pheno$Line
  }

  if (with_markers) {
    geno_lines <- rownames(Geno)
  } else {
    geno_lines <- colnames(Geno)
  }
  pheno_lines <- Pheno$Line

  lines_in_both <- unique(geno_lines[which(geno_lines %in% pheno_lines)])

  Pheno <<- subset(Pheno, subset = Line %in% lines_in_both)
  # Reseet factor levels
  Pheno <<- droplevels(Pheno)

  geno_indices <- which(geno_lines %in% lines_in_both)
  if (with_markers) {
    Geno <<- Geno[geno_indices, ]
  } else {
    Geno <<- Geno[geno_indices, geno_indices]
  }
  Geno <<- as.matrix(Geno)

  return(list(Pheno=Pheno, Geno=Geno))
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
                          verbose) {
  is_uni_env <- !is.null(env)
  if (is_uni_env) {
    if (env == "PRN" && .Platform$OS.type == "windows") {
      env <- "PRN_ENV"
    }
    results_dir <- file.path(results_dir, env)
  }

  tmp <- prepare_geno_pheno(Pheno, Geno, is_uni_env, env)
  Pheno <- tmp$Pheno
  Geno <- tmp$Geno

  X <- prepare_X(Pheno = Pheno, Geno = Geno, with_markers = with_markers,
                 with_interaction = with_interaction,
                 is_uni_env = is_uni_env)
  mtry_theta <- prepare_mtry(mtry_theta, ncol(X))
  mtry_lambda <- prepare_mtry(mtry_lambda, ncol(X))

  Results <- list()

  if (type_of_tuning == "global") {
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
      cross_validation = tunig_cross_validation,
      proportion_of_testing = tuning_proportion_of_testing,
      sample_proportion = sample_proportion,
      verbose = verbose)

    best_params <- tuning_results$best_params
  }

  folds <- get_folds(cross_validation, number_of_folds, proportion_of_testing,
                     n_records = length(y))
  for (n_current_fold in 1:number_of_folds) {
    fold <- folds[[n_current_fold]]
    X_training <- X[fold$training, ]
    y_training <- y[fold$training]
    X_testing <- X[fold$testing, ]
    y_testing <- y[fold$testing]

    if (type_of_tuning == "local") {
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
        cross_validation = tunig_cross_validation,
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
  }
}

gene.zap.rfsrc <- function(Pheno, y, Geno = NULL, Markers = NULL,
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
    Pheno, y, Geno, Markers, with_interaction, mult_env_anal,
    ntree_theta, mtry_theta, nodesize_theta,
    ntree_lambda, mtry_lambda, nodesize_lambda,
    importance, type, loss_function,
    cross_validation, number_of_folds, proportion_of_testing,
    type_of_tuning, tuning_cross_validation, tuning_number_of_folds,
    tuning_proportion_of_testing, sample_proportion,
    results_dir, verbose, digits)

  with_markers <- !is.null(Markers)
  if (is.null(Geno)) {
    Geno <- Markers
    Markers <- NULL
  }

  if (mult_env_anal) {
    train_gen_zap(
      y = Pheno[[trait_name]], trait_name = trait_name,
      print_level = print_level + 1, env = NULL
    )
  } else {
    envs <- as.character(unique(Pheno$Env))
    for (env in envs) {
      if (verbose) {
        cat(env, ":\n", sep = "")
      }

      train_gen_zap(
        y = Pheno[[trait_name]], trait_name = trait_name, env = env,
        print_level = print_level + 2
      )
    }
  }
}