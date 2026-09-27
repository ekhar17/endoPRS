#' Function to fit endoPRS-SS for a single threshold pair and weight combination
#'
#' This function fits one cell of the endoPRS-SS tuning grid and writes its validation performance
#' and effect sizes to disk. It is the summary statistic analogue of \code{fit_endoPRS_single_iter},
#' and exists so that the grid can be fit in parallel as a job array on a high performance computing
#' cluster. Use it together with \code{select_endoPRS_ss}.
#'
#' Note that the full lassosum2 (lambda, delta) grid is still fit within a single call, because
#' those models share the same LD matrix and are cheap to fit once it has been built.
#'
#' @inheritParams fit_endoPRS_ss
#' @param thresh The p-value threshold to apply to the phenotype GWAS.
#' @param thresh_endo The p-value threshold to apply to the endophenotype GWAS.
#' @param iteration A numerical value corresponding to which row of the grid of weights to fit.
#' @param save_folder A path to a directory that files can be written to. The validation performance and effect sizes of this model are saved there.
#'
#' @return A list with two elements:
#' \itemize{
#' \item{beta: A data frame of the effect sizes of the best (lambda, delta) combination for this cell of the grid.}
#' \item{val_perf: A one row data frame of the validation performance of that model.}
#' }
#'
#' @export
fit_endoPRS_ss_single_iter = function(G, map, fam,
                                      val_pheno, val_covar = NULL,
                                      pheno_gwas, endo_gwas,
                                      ld_ref,
                                      filter_hapmap = FALSE, hapmap = NULL, type = NULL,
                                      thresh, thresh_endo = NULL,
                                      grid = NULL, iteration,
                                      delta = c(0.001, 0.01, 0.1, 1), nlambda = 30,
                                      NCORES = NULL,
                                      save_folder){

  if(is.null(thresh_endo)) thresh_endo = thresh

  ## Allowed, but it is outside the grid fit_endoPRS_ss tunes over, so it is most likely a swapped argument
  if(thresh_endo > thresh){
    warning(paste0("thresh_endo (", thresh_endo, ") is less stringent than thresh (", thresh,
                   "). fit_endoPRS_ss only tunes over thresh_endo <= thresh; check that the",
                   " arguments are not swapped."))
  }

  if(is.null(NCORES)) NCORES = nb_cores()
  if(is.null(grid)){
    grid = expand.grid(w2 = c(1e-1, 0.5, 1, 2, 10), w3 = c(1e-1, 0.5, 1, 2, 10))
  }

  ## Validation individuals
  geno_id = paste(fam[, 1], fam[, 2], sep = "_")
  val_id = paste(val_pheno[, 1], val_pheno[, 2], sep = "_")
  if(sum(val_id %in% geno_id) == 0){
    stop("There is no overlap between validation ID's and genotype ID's. Please reformat.")
  }
  val_index = match(val_id, geno_id)
  y.val = c(val_pheno[, 3])

  if(is.null(type)){
    type = if(sum(!(y.val %in% c(0, 1))) > 0) "linear" else "logistic"
  }
  covar.val = if(is.null(val_covar)) NULL else covar_from_df(val_covar[, -c(1, 2), drop = FALSE])

  geno_map = data.frame(chr = as.numeric(map$chromosome), pos = map$physical.pos,
                        a0 = map$allele2, a1 = map$allele1)

  print(paste("Fitting model corresponding to thresh:", thresh, "thresh_endo:", thresh_endo,
              "and weights:"))
  print(grid[iteration, ])

  ## Align summary statistics, select variants, and build the LD matrix
  pheno_gwas_matched = format_gwas_ss(pheno_gwas, ld_ref$map, filter_hapmap, hapmap)
  endo_gwas_matched  = format_gwas_ss(endo_gwas,  ld_ref$map, filter_hapmap, hapmap)

  snp_groups = extract_snp_groups_ss(pheno_gwas_matched, endo_gwas_matched, thresh, thresh_endo)
  df_beta = pheno_gwas_matched[pheno_gwas_matched$`_NUM_ID_` %in% snp_groups$`_NUM_ID_`, ]
  corr = build_corr_ss(df_beta, ld_ref)

  geno_align = df_beta[, c("chr", "pos", "a0", "a1")]
  geno_align$beta = 1
  geno_align$ss_row = seq_len(nrow(df_beta))
  geno_align = snp_match(geno_align, geno_map, join_by_pos = TRUE)

  ## Fit the weighted model and score it in the validation set
  penalty_table = create_penalty_table_ss(snp_groups, df_beta,
                                          w2 = grid$w2[iteration], w3 = grid$w3[iteration])

  beta_grid = weighted_lassosum2(corr, df_beta, penalty_table$penalty,
                                 delta = delta, nlambda = nlambda, ncores = NCORES)
  params = attr(beta_grid, "grid_param")

  pred = big_prodMat(G,
                     beta_grid[geno_align$ss_row, , drop = FALSE] * geno_align$beta,
                     ind.row = val_index, ind.col = geno_align$`_NUM_ID_`, ncores = NCORES)

  params$val_res = score_pred_ss(pred, y.val, covar.val, type)
  params$thresh = thresh
  params$thresh_endo = thresh_endo
  params$w2 = grid$w2[iteration]
  params$w3 = grid$w3[iteration]

  if(all(is.na(params$val_res))){
    stop("No (lambda, delta) combination produced a usable score in the validation set.")
  }

  best = which.max(params$val_res)
  res = params[best, ]
  beta_info = data.frame(df_beta[, c("chr", "pos", "a0", "a1")], beta = beta_grid[, best])
  beta_info = beta_info[!is.na(beta_info$beta) & beta_info$beta != 0, ]

  ## Save so that select_endoPRS_ss can pick the winner across jobs
  stem = paste0("endoPRS_ss_thresh", thresh, "_threshendo", thresh_endo,
                "_w2", grid$w2[iteration], "_w3", grid$w3[iteration])
  write.csv(res, file.path(save_folder, paste0(stem, "_validationresults.csv")), row.names = FALSE)
  write.csv(beta_info, file.path(save_folder, paste0(stem, "_beta.csv")), row.names = FALSE)

  return(list(beta = beta_info, val_perf = res))
}


#' Function to select the best endoPRS-SS model fit across a job array
#'
#' This function reads the validation results written by \code{fit_endoPRS_ss_single_iter}, picks
#' the best performing model, and returns its effect sizes. It is the summary statistic analogue of
#' \code{refit_endoPRS}.
#'
#' @param save_folder The directory that contains the results written by fit_endoPRS_ss_single_iter.
#'
#' @return A list with three elements: beta, best_params, and val_results, matching the output of fit_endoPRS_ss.
#'
#' @export
select_endoPRS_ss = function(save_folder){

  files = list.files(save_folder, pattern = "_validationresults\\.csv$", full.names = TRUE)
  if(length(files) == 0){
    stop(paste("No validation result files found in", save_folder,
               "- please check that fit_endoPRS_ss_single_iter has been run."))
  }

  val_results = do.call(rbind, lapply(files, read.csv))
  best_params = val_results[which.max(val_results$val_res), ]

  print(paste("Best performing model corresponds to thresh:", best_params$thresh,
              "thresh_endo:", best_params$thresh_endo,
              "w2:", best_params$w2, "w3:", best_params$w3))

  beta_file = file.path(save_folder,
                        paste0("endoPRS_ss_thresh", best_params$thresh,
                               "_threshendo", best_params$thresh_endo,
                               "_w2", best_params$w2, "_w3", best_params$w3, "_beta.csv"))
  if(!file.exists(beta_file)){
    stop(paste("Cannot find the effect sizes of the selected model at", beta_file))
  }

  return(list(beta = read.csv(beta_file), best_params = best_params, val_results = val_results))
}
