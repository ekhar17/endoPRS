#' Function to score a grid of candidate PRS in the validation set
#'
#' Internal helper used by the endoPRS-SS functions. Covariates, when supplied, are regressed out of
#' the predicted scores so that the reported performance is not driven by them.
#'
#' @param pred A matrix of predicted scores, one column per candidate model.
#' @param y The validation phenotype.
#' @param covar A matrix of validation covariates, or NULL.
#' @param type Either "linear" or "logistic".
#'
#' @return A numeric vector with one performance value per column of pred. R^2 is used for linear models and AUC for logistic models.
#'
score_pred_ss = function(pred, y, covar, type){

  apply(pred, 2, function(x){

    ## lassosum2 returns NA for grid points where the solver diverged
    if(all(is.na(x)) | stats::sd(x, na.rm = TRUE) == 0) return(NA_real_)

    if(!is.null(covar)) x = stats::residuals(stats::lm(x ~ covar))

    if(type == "linear")  return(stats::cor(x, y)^2)
    if(type == "logistic") return(AUC(x, y))
  })
}


#' Function to apply endoPRS-SS to generate a polygenic risk score model from summary statistics
#'
#' This function applies the summary statistic version of endoPRS. Like the individual level
#' version, it penalizes SNPs differently depending on whether they are associated with only the
#' phenotype, only the endophenotype, or both. Unlike the individual level version it never needs
#' individual level training genotypes: the model is fit from GWAS summary statistics and an LD
#' reference panel using a weighted lassosum2 model. Individual level data is only needed for the
#' validation set, which is used to tune the p-value thresholds, the endophenotype weights, and the
#' lassosum2 regularization parameters.
#'
#' endoPRS-SS also allows the phenotype and the endophenotype to use different p-value thresholds.
#' The individual level version of endoPRS applies a single threshold to both traits.
#'
#' @param G An object of class FBM from the bigsnpr package that contains the genotypes of the individuals used for validation.
#' @param map A data frame containing information about SNPs in G (genotype matrix). This can be the map object from the bigSNP class. It should contain columns labeled chromosome, physical.pos, allele1, and allele2.
#' @param fam A data frame containing information about the individuals in G (genotype matrix). Column one must correspond to FID and column 2 must correspond to IID.
#' @param val_pheno A data frame with three columns corresponding to the phenotypes of individuals used for validation. Column one is FID, column two is IID, and column three is y the phenotype. If the phenotype is binary, y must consist of 0's and 1's.
#' @param val_covar An optional data frame of covariates of the individuals used for validation (ie genetic PC's, sex, age, etc). The first column must be FID, and the second column must be IID. The order of individuals must be the same as val_pheno. If supplied, the covariates are regressed out of the predicted scores before the validation performance is computed.
#' @param pheno_gwas A data frame containing the results of the GWAS run on the phenotype. It must not include any individuals from the validation set. It must contain columns corresponding to chromosome, position, allele0, allele1, effect size, standard error, and sample size.
#' @param endo_gwas A data frame containing the results of the GWAS run on the endophenotype. It must not include any individuals from the validation set. It must contain the same columns as pheno_gwas.
#' @param ld_ref A list with elements map and corr_files describing the LD reference panel, as produced by load_ld_ref.
#' @param filter_hapmap An optional logical value. It corresponds to whether to only run the model on variants in the hapmap3 set.
#' @param hapmap An optional data frame with hapmap variants. This must be provided if filter_hapmap is set to TRUE. It should contain columns CHR and POS.
#' @param type An optional character vector of "linear" or "logistic." If not provided, it will be learned from the phenotype.
#' @param threshes An optional vector of p-value thresholds applied to the phenotype GWAS. If not provided, the values from the endoPRS manuscript will be used (0.01, 1e-4, and 1e-6).
#' @param threshes_endo An optional vector of p-value thresholds applied to the endophenotype GWAS. Every combination with threshes for which thresh_endo is at least as stringent as thresh is tuned over. If not provided, threshes is used.
#' @param grid An optional data frame of weights to use for the weighted penalty. The first column must be w2 and correspond to the weights applied to SNPs associated with only the endophenotype, and the second column must be w3 and correspond to the weights applied to SNPs associated with both traits. SNPs associated with only the phenotype are given a weight of 1. If not provided, the grid from the endoPRS manuscript (0.1, 0.5, 1, 2, 10) is used for both.
#' @param delta An optional vector of lassosum2 L2 regularization parameters. Default is c(0.001, 0.01, 0.1, 1).
#' @param nlambda An optional number of lassosum2 L1 regularization parameters to try. Default is 30.
#' @param NCORES An optional value corresponding to the number of cores to use. Otherwise, the number of cores will be learned using nb_cores().
#' @param pheno_gwas_refit An optional data frame containing the results of the GWAS run on the phenotype using both the training and validation set. If supplied together with endo_gwas_refit, the selected model is refit using these summary statistics.
#' @param endo_gwas_refit An optional data frame containing the results of the GWAS run on the endophenotype using both the training and validation set.
#' @param save_folder An optional path to a directory that files can be written to. If specified, the validation results are written there as they are produced.
#'
#' @return A list with three elements:
#' \itemize{
#' \item{beta: A data frame of the SNPs included in the final model and their coefficients, on the allele coding of the LD reference panel. This can be used to apply the PRS to an external data set using software such as PLINK.}
#' \item{best_params: A one row data frame giving the selected thresh, thresh_endo, w2, w3, lambda, and delta, and the validation performance.}
#' \item{val_results: A data frame of the validation performance of every model that was fit.}
#' }
#'
#' @export
fit_endoPRS_ss = function(G, map, fam,
                          val_pheno, val_covar = NULL,
                          pheno_gwas, endo_gwas,
                          ld_ref,
                          filter_hapmap = FALSE, hapmap = NULL, type = NULL,
                          threshes = c(1e-2, 1e-4, 1e-6),
                          threshes_endo = NULL,
                          grid = NULL,
                          delta = c(0.001, 0.01, 0.1, 1), nlambda = 30,
                          NCORES = NULL,
                          pheno_gwas_refit = NULL, endo_gwas_refit = NULL,
                          save_folder = NULL){

  #######################################################################################
  ###########             Part 1: Check and Process Input               #################
  #######################################################################################

  if(any(is.na(val_pheno))){
    stop("Currently endoPRS cannot handle any missing values in the validation phenotype. Please reformat.")
  }
  if(!is.null(val_covar)){
    if(any(is.na(val_covar))){
      stop("Currently endoPRS cannot handle any missing values in the validation covariates. Please reformat.")
    }
    if(!isTRUE(all.equal(val_pheno[, 1:2], val_covar[, 1:2]))){
      stop("First two columns of val pheno and val covar do not match. Please reformat.")
    }
  }

  if(nrow(G) != nrow(fam)){
    stop("Number of rows of G does not match the number of rows of fam. Please check your input.")
  }
  if(ncol(G) != nrow(map)){
    stop("Number of columns of G does not match the number of rows of map. Please check your input.")
  }

  if(is.null(NCORES)) NCORES = nb_cores()

  ## Grid of endophenotype weights
  if(is.null(grid)){
    grid = expand.grid(w2 = c(1e-1, 0.5, 1, 2, 10), w3 = c(1e-1, 0.5, 1, 2, 10))
  }

  ## Pairs of p-value thresholds. The endophenotype threshold is never looser than the phenotype
  ## one, which is the set of combinations used in the endoPRS-SS manuscript.
  if(is.null(threshes_endo)) threshes_endo = threshes
  thresh_pairs = expand.grid(thresh = threshes, thresh_endo = threshes_endo)
  thresh_pairs = thresh_pairs[thresh_pairs$thresh_endo <= thresh_pairs$thresh, ]

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

  covar.val = if(is.null(val_covar)) NULL else covar_from_df(val_covar[, -c(1, 2)])

  ## Genotype map in the allele naming used by snp_match. bigsnpr's convention is that allele1 is
  ## the effect allele (a1) and allele2 is the other allele (a0).
  geno_map = data.frame(chr = as.numeric(map$chromosome), pos = map$physical.pos,
                        a0 = map$allele2, a1 = map$allele1)

  #######################################################################################
  ###########      Part 2: Align Summary Statistics to LD Reference     #################
  #######################################################################################

  print("Aligning summary statistics to the LD reference panel.")
  pheno_gwas_matched = format_gwas_ss(pheno_gwas, ld_ref$map, filter_hapmap, hapmap)
  endo_gwas_matched  = format_gwas_ss(endo_gwas,  ld_ref$map, filter_hapmap, hapmap)

  #######################################################################################
  ###########     Part 3: Fit Grid of Models and Tune on Validation     #################
  #######################################################################################

  val_results = data.frame()
  beta_store = list()

  print(paste("Fitting", nrow(grid), "weight combinations over",
              nrow(thresh_pairs), "threshold pairs."))

  for(p in seq_len(nrow(thresh_pairs))){

    thresh = thresh_pairs$thresh[p]
    thresh_endo = thresh_pairs$thresh_endo[p]

    ## Assign each selected variant to the pheno only, endo only, or both group
    snp_groups = extract_snp_groups_ss(pheno_gwas_matched, endo_gwas_matched, thresh, thresh_endo)

    ## Effect sizes come from the phenotype GWAS for every selected variant
    df_beta = pheno_gwas_matched[pheno_gwas_matched$`_NUM_ID_` %in% snp_groups$`_NUM_ID_`, ]

    ## The LD matrix depends only on the selected variants, so it is built once per threshold pair
    ## rather than once per weight combination.
    corr = build_corr_ss(df_beta, ld_ref)

    ## Map the selected variants onto the columns of G. snp_match returns the sign needed to put
    ## each effect size onto the genotype allele coding.
    geno_align = df_beta[, c("chr", "pos", "a0", "a1")]
    geno_align$beta = 1
    geno_align$ss_row = seq_len(nrow(df_beta))
    geno_align = snp_match(geno_align, geno_map, join_by_pos = TRUE)

    print(paste0("thresh: ", thresh, ", thresh_endo: ", thresh_endo, " - ",
                 nrow(df_beta), " variants selected, ", nrow(geno_align),
                 " of them present in the validation genotypes."))
    print(table(snp_groups$group))

    pb = txtProgressBar(min = 0, max = nrow(grid), initial = 0, style = 3)

    for(iter in seq_len(nrow(grid))){

      ## Weighted penalty for this combination of weights
      penalty_table = create_penalty_table_ss(snp_groups, df_beta,
                                              w2 = grid$w2[iter], w3 = grid$w3[iter])

      ## Fit the weighted lassosum2 model over its own (lambda, delta) grid
      beta_grid = weighted_lassosum2(corr, df_beta, penalty_table$penalty,
                                     delta = delta, nlambda = nlambda, ncores = NCORES)
      params = attr(beta_grid, "grid_param")

      ## Score every candidate model in the validation set
      pred = big_prodMat(G,
                         beta_grid[geno_align$ss_row, , drop = FALSE] * geno_align$beta,
                         ind.row = val_index,
                         ind.col = geno_align$`_NUM_ID_`,
                         ncores = NCORES)

      params$val_res = score_pred_ss(pred, y.val, covar.val, type)
      params$thresh = thresh
      params$thresh_endo = thresh_endo
      params$w2 = grid$w2[iter]
      params$w3 = grid$w3[iter]
      params$iter = iter

      ## Keep only the best model from this fit so that the full beta grid is not held in memory
      if(all(is.na(params$val_res))){
        setTxtProgressBar(pb, iter)
        next
      }
      best = which.max(params$val_res)
      beta_store[[paste(thresh, thresh_endo, iter, sep = "_")]] =
        data.frame(df_beta[, c("chr", "pos", "a0", "a1")], beta = beta_grid[, best])

      val_results = rbind(val_results, params[best, ])

      setTxtProgressBar(pb, iter)
    }

    close(pb)

    ## The LD matrix for this threshold pair is no longer needed, and at biobank scale these
    ## backing files are large enough that leaving them until the session ends matters.
    unlink(corr$backingfile)
    rm(corr)

    if(!is.null(save_folder)){
      write.csv(val_results, file.path(save_folder, "endoPRS_ss_validation_results.csv"),
                row.names = FALSE)
    }
  }

  if(nrow(val_results) == 0){
    stop("No model produced a usable score in the validation set. Please check your input.")
  }

  #######################################################################################
  ###########           Part 4: Select and Optionally Refit             #################
  #######################################################################################

  best_params = val_results[which.max(val_results$val_res), ]
  print(paste("Best performing model corresponds to thresh:", best_params$thresh,
              "thresh_endo:", best_params$thresh_endo,
              "w2:", best_params$w2, "w3:", best_params$w3,
              "lambda:", signif(best_params$lambda, 3), "delta:", best_params$delta))

  beta_info = beta_store[[paste(best_params$thresh, best_params$thresh_endo,
                                best_params$iter, sep = "_")]]

  ## If GWAS run on the combined training and validation set are supplied, refit the selected model
  ## on them. This mirrors the refitting step of the individual level version of endoPRS.
  if(!is.null(pheno_gwas_refit) & !is.null(endo_gwas_refit)){

    print("Refitting the selected model using the combined training and validation summary statistics.")

    pheno_refit_matched = format_gwas_ss(pheno_gwas_refit, ld_ref$map, filter_hapmap, hapmap)
    endo_refit_matched  = format_gwas_ss(endo_gwas_refit,  ld_ref$map, filter_hapmap, hapmap)

    snp_groups = extract_snp_groups_ss(pheno_refit_matched, endo_refit_matched,
                                       best_params$thresh, best_params$thresh_endo)
    df_beta = pheno_refit_matched[pheno_refit_matched$`_NUM_ID_` %in% snp_groups$`_NUM_ID_`, ]
    corr = build_corr_ss(df_beta, ld_ref)

    penalty_table = create_penalty_table_ss(snp_groups, df_beta,
                                            w2 = best_params$w2, w3 = best_params$w3)

    ## Refit at the selected lambda and delta only
    beta_grid = weighted_lassosum2(corr, df_beta, penalty_table$penalty,
                                   delta = best_params$delta, nlambda = nlambda,
                                   ncores = NCORES)
    params = attr(beta_grid, "grid_param")
    closest = which.min(abs(log(params$lambda) - log(best_params$lambda)))

    beta_info = data.frame(df_beta[, c("chr", "pos", "a0", "a1")], beta = beta_grid[, closest])
  }

  ## Drop the variants that were shrunk to zero
  beta_info = beta_info[!is.na(beta_info$beta) & beta_info$beta != 0, ]

  if(!is.null(save_folder)){
    write.csv(beta_info, file.path(save_folder, "endoPRS_ss_beta.csv"), row.names = FALSE)
    write.csv(val_results, file.path(save_folder, "endoPRS_ss_validation_results.csv"),
              row.names = FALSE)
  }

  return(list(beta = beta_info, best_params = best_params, val_results = val_results))
}
