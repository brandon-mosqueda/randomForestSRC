is_number <- function(value) {
  if (is.factor(value)) {
    return(FALSE)
  }

  suppressWarnings(value <- as.numeric(value))

  return(!is.na(value))
}

is_int <- function(number) {
  return(is_number(number) && as.numeric(number) %% 1 == 0)
}

is_discrete <- function(number) {
  if (!all(is_int(number))) {
    return(FALSE)
  }

  return(all(number >= 0))
}

get_formula <- function(response, predictors) {
  return(formula(paste0(response, " ~ ",
                        paste0(predictors,  collapse = " + "))))
}

#' @title Cholesky
#'
#' @description Compute the Cholesky factorization of a non-real symmetric
#'              positive-definite square matrix.
#'
#' @param G (\code{numeric - matrix}) an object to apply this method, it could
#'        be non positive-definite matrices.
#' @param tolerance (\code{double}) Tolerance level, by default is 1e-10.
#'
#' @export
cholesky <- function(G, tolerance = 1e-10) {
  data_names <- colnames(G)

  tryCatch({
    result <- t(chol(G))
  }, error=function(error_condition) {
      G <- (G + t(G)) / 2
      EigenA <- eigen(G)
      d_A    <- EigenA$values
      V_A    <- EigenA$vectors
      d_A[which(d_A < tolerance)] <- tolerance*100L
      pos_A1 <- which(d_A > tolerance)
      if (identical(pos_A1, integer(0))) {
        pos_A <- 1L
      } else {
        pos_A <- pos_A1
      }
      d_A_Star <- d_A[pos_A]
      V_A_Star <- V_A[, pos_A]

      if (length(pos_A) == 1L) {
        d_A_Star <- 1L / d_A_Star
        LG <- d_A_Star * sqrt(V_A_Star)
      } else {
        d_A_Star <- diag(d_A_Star)
        LG <- (V_A_Star %*% sqrt(d_A_Star))
      }

      result <- LG
  })

  rownames(result) <- data_names
  colnames(result) <- data_names

  return(result)
}

prepare_X_with_markers <- function(Pheno, Geno, with_interaction, is_uni_env) {
  ZG <- model.matrix(~0 + as.factor(Pheno$Line))
  ZG1 <- ZG %*% Geno
  X <- ZG1

  # It's important to respect the order: Env, Lines and Interaction effect
  if (!is_uni_env) {
    #####of environment and GxE###################
    Z.E <- model.matrix(~0 + as.factor(Pheno$Env))
    X <- cbind(Z.E, X)

    if (with_interaction) {
      ####Linear kernel for the GE term
      Z.GE <- model.matrix(~0 + ZG1:as.factor(Pheno$Env))
      X <- cbind(X, Z.GE)
    }
  }

  return(X)
}

prepare_X_without_markers <- function(Pheno, Geno, with_interaction,
                                      is_uni_env) {
  # Of lines
  ZG <- model.matrix(~0 + as.factor(Pheno$Line))

  # Creating the linear kernel of lines for RKHS
  ZL <- cholesky(Geno)
  ZG1 <- ZG %*% ZL
  X <- ZG1

  # It's important to respect the order: Env, Lines and Interaction effect
  if (!is_uni_env) {
    # Of environment
    Z.E <- model.matrix(~0 + as.factor(Pheno$Env))
    X <- cbind(Z.E, X)

    if (with_interaction) {
      # GxE
      Z.GE <- model.matrix(~0 + ZG1:as.factor(Pheno$Env))
      X <- cbind(X, Z.GE)
    }
  }

  return(X)
}

prepare_X <- function(Pheno, Geno, with_markers, with_interaction,
                      is_uni_env) {
  if (with_markers) {
    X <- prepare_X_with_markers(Pheno, Geno, with_interaction, is_uni_env)
  } else {
    X <- prepare_X_without_markers(Pheno, Geno, with_interaction, is_uni_env)
  }

  X <- as.data.frame(X)
  colnames(X) <- make.names(colnames(X))
  rownames(X) <- make.names(rownames(X))

  return(X)
}

#' @title Mean square error of prediction
mse <- function(actual, predicted) {
  return(mean((actual - predicted)^2, na.rm=TRUE))
}

get_folds <- function(cross_validation, number_of_folds,
                      proportion_of_testing, n_records) {
  folds <- list()

  if (cross_validation == "k_fold") {
      folds_vector <- findInterval(cut(sample(1:n_records, n_records),
                                       breaks=number_of_folds), 1:n_records)

      for (fold_num in 1:number_of_folds) {
        current_fold <- list()
        current_fold$testing <- which(folds_vector == fold_num)
        current_fold$training <- setdiff(1:n_records, current_fold$testing)

        folds[[fold_num]] <- current_fold
      }

      return(folds)

  } else if (cross_validation == "random_partition") {
    temp_folds <- replicate(number_of_folds,
                            sample(n_records,
                                   ceiling(proportion_of_testing * n_records)))

    for (fold_num in 1:number_of_folds) {
      current_fold <- list(testing=temp_folds[, fold_num])
      current_fold$training <- setdiff(1:n_records, current_fold$testing)

      folds[[fold_num]] <- current_fold
    }
  }

  return(folds)
}

get_combinations <- function(ntree_theta, mtry_theta,
                             nodesize_theta, ntree_lambda,
                             mtry_lambda, nodesize_lambda,
                             sample_proportion, type) {
  FLAGS <- list(ntree_theta = ntree_theta, mtry_theta = mtry_theta,
                nodesize_theta = nodesize_theta, ntree_lambda = ntree_lambda,
                mtry_lambda = mtry_lambda, nodesize_lambda = nodesize_lambda,
                type = type)
  all_params_combinations <- FLAGS %>% cross()

  if (sample_proportion == 1) {
    return(all_params_combinations)
  }

  n_combinations <- length(all_params_combinations)
  final_combinations <- sample(n_combinations,
                               ceiling(n_combinations * sample_proportion))

  return(all_params_combinations[final_combinations])
}

get_function_name <- function(funct) {
  return(make.names(as.character(substitute(funct))))
}

get_rand_state <- function() {
  # Using `get0()` here to have `NULL` output in case object doesn't
  # exist. Also using `inherits = FALSE` to get value exactly from
  # global environment and not from one of its parent.
  get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
}

set_rand_state <- function(state) {
  # Assigning `NULL` state might lead to unwanted consequences
  if (!is.null(state)) {
    assign(".Random.seed", state, envir = .GlobalEnv, inherits = FALSE)
  }
}