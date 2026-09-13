#' Function to fit a weighted lassosum2 model over a grid of (lambda, delta)
#'
#' This function is the summary statistic analogue of \code{fit_weighted_model}. It fits the
#' lassosum2 penalized regression model of Prive et al. using only GWAS summary statistics and
#' an LD reference panel, but it penalizes each variant differently depending on whether it is
#' associated with only the phenotype, only the endophenotype, or both.
#'
#' It reproduces \code{bigsnpr::snp_lassosum2} exactly except that the per-variant penalty factor
#' is multiplied by the endoPRS weights. The compiled lassosum2 solver in bigsnpr already accepts
#' per-variant \code{lambda} and \code{delta_plus_one} vectors, so no modification of bigsnpr is
#' required -- this function calls that solver directly.
#'
#' @param corr An object of class SFBM containing the correlation (LD) matrix restricted and ordered to match the rows of df_beta. This is the output of build_corr_ss.
#' @param df_beta A data frame of summary statistics containing at minimum the columns beta, beta_se, and n_eff. Rows must be in the same order as the columns of corr and as penalty.
#' @param penalty A numeric vector of strictly positive per-variant penalty factors, one per row of df_beta. This is the penalty column of the output of create_penalty_table.
#' @param delta A vector of shrinkage parameters to try (L2 regularization). Default is c(0.001, 0.01, 0.1, 1).
#' @param nlambda Number of different lambdas to try (L1 regularization). Default is 30.
#' @param lambda.min.ratio Ratio between last and first lambdas to try. Default is 0.01.
#' @param dfmax Maximum number of non-zero effects in the model. Default is 200e3.
#' @param maxiter Maximum number of iterations before convergence. Default is 1000.
#' @param tol Tolerance parameter for assessing convergence. Default is 1e-5.
#' @param ind.corr Indices of the columns of corr to use. Defaults to all of them.
#' @param ncores Number of cores to use. Default is 1.
#'
#' @return A matrix of effect sizes with one column per row of \code{attr(<res>, "grid_param")}.
#'   Missing values are returned for grid points where strong divergence was detected.
#'
weighted_lassosum2 = function(corr, df_beta, penalty,
                              delta = c(0.001, 0.01, 0.1, 1),
                              nlambda = 30, lambda.min.ratio = 0.01,
                              dfmax = 200e3, maxiter = 1000, tol = 1e-5,
                              ind.corr = cols_along(corr),
                              ncores = 1){

  ## Same input checks as bigsnpr::snp_lassosum2
  assert_df_with_names(df_beta, c("beta", "beta_se", "n_eff"))
  assert_lengths(ind.corr, rows_along(df_beta))
  stopifnot(all(ind.corr %in% cols_along(corr)))
  assert_pos(df_beta$beta_se, strict = TRUE)
  assert_pos(delta, strict = TRUE)
  assert_cores(ncores)

  ## endoPRS specific check - one strictly positive penalty per variant
  assert_lengths(penalty, rows_along(df_beta))
  assert_pos(penalty, strict = TRUE)

  ## Put the effect sizes on the standardized scale used by lassosum2
  N <- df_beta$n_eff
  scale <- sqrt(N * df_beta$beta_se^2 + df_beta$beta^2)
  beta_hat <- df_beta$beta / scale

  ## This is the only line that differs from bigsnpr::snp_lassosum2. In lassosum2 the penalty
  ## factor corrects for differing per-variant sample sizes; here it additionally carries the
  ## endoPRS weights, so variants are shrunk differently by SNP group.
  pf <- sqrt(max(N) / N) * penalty

  ## Build the (lambda, delta) grid
  lambda0 <- max(abs(beta_hat / pf))
  seq_lam <- seq_log(lambda0, lambda.min.ratio * lambda0, nlambda + 1)[-1]
  grid_param <- expand.grid(lambda = seq_lam, delta = delta)

  ## Fit from least to most sparse so that warm starts are cheap
  ord <- with(grid_param, order(lambda * (1 + delta)))
  inv_ord <- match(seq_along(ord), ord)

  bigparallelr::register_parallel(ncores)

  res_grid <- foreach(ic = ord, .packages = "bigsnpr") %dopar% {

    ## lassosum2 is the compiled solver inside bigsnpr. It is not exported, and it is looked up
    ## inside the worker rather than passed in so that this also works on PSOCK clusters, where
    ## the external pointer to the compiled routine would not survive serialization.
    lassosum2_solver <- get("lassosum2", envir = asNamespace("bigsnpr"))

    time <- system.time(
      res <- lassosum2_solver(
        corr           = corr,
        beta_hat       = beta_hat,
        lambda         = pf * grid_param$lambda[ic],
        delta_plus_one = pf * grid_param$delta[ic] + 1,
        ind_sub        = ind.corr - 1L,
        dfmax          = dfmax,
        maxiter        = maxiter,
        tol            = tol
      )
    )

    res$time <- time[["elapsed"]]
    res
  }
  res_grid <- res_grid[inv_ord]

  grid_param$num_iter <- sapply(res_grid, function(.) .$num_iter)
  grid_param$time <- sapply(res_grid, function(.) .$time)
  beta_grid <- do.call("cbind", lapply(res_grid, function(.) .$beta_est))
  grid_param$sparsity <- colMeans(beta_grid == 0)

  ## Return effect sizes on the original scale
  structure(sweep(beta_grid, 1, scale, '*'), grid_param = grid_param)
}
