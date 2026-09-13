#' Function to clean GWAS summary statistics and align them to an LD reference panel
#'
#' This function is the summary statistic analogue of \code{format_gwas}. In addition to the
#' columns needed by the individual level version it also standardizes the effect size, standard
#' error and sample size columns, and it aligns the summary statistics to the LD reference panel
#' using \code{bigsnpr::snp_match}.
#'
#' Aligning to the LD reference (rather than to the genotype map) means \code{snp_match} performs
#' the allele flipping, so the correlation matrix itself never has to be sign flipped.
#'
#' @param df A data frame of GWAS results. It must contain columns for chromosome, position, allele0, allele1, effect size, standard error, and sample size. A P-value column is used if present and otherwise computed from beta and se.
#' @param ld_map A data frame describing the variants of the LD reference panel, with columns chr, pos, a0, a1. Its rows must be in the same order as the rows and columns of the LD correlation matrices.
#' @param filter_hapmap A logical value as to whether to filter to SNPs in hapmap.
#' @param hapmap A data frame with hapmap variants. Only provide if filter_hapmap is set to TRUE. It should contain columns CHR, POS, A1, A2.
#'
#' @return A data frame of summary statistics aligned to the LD reference, containing the columns
#'   chr, pos, a0, a1, beta, beta_se, n_eff, P, and `_NUM_ID_` (the row of ld_map that each variant
#'   corresponds to). Rows are sorted by `_NUM_ID_`.
#'
format_gwas_ss = function(df, ld_map, filter_hapmap = FALSE, hapmap = NULL){

  df = as.data.frame(df)

  ## Properly format the columns needed for a summary statistic model
  colnames(df)[toupper(colnames(df)) %in% c("CHROM", "CHR")] = "chr"
  colnames(df)[toupper(colnames(df)) %in% c("GENPOS", "POS", "BP", "PHYSICAL.POS", "POSITION")] = "pos"
  colnames(df)[toupper(colnames(df)) %in% c("ALLELE0", "A0", "ALLELE2", "A2", "OTHER_ALLELE")] = "a0"
  colnames(df)[toupper(colnames(df)) %in% c("ALLELE1", "A1", "EFFECT_ALLELE")] = "a1"
  colnames(df)[toupper(colnames(df)) %in% c("BETA", "EFFECT", "B")] = "beta"
  colnames(df)[toupper(colnames(df)) %in% c("SE", "BETA_SE", "STANDARD_ERROR")] = "beta_se"
  colnames(df)[toupper(colnames(df)) %in% c("N", "N_EFF", "OBS_CT")] = "n_eff"
  colnames(df)[toupper(colnames(df)) %in% c("P", "PVAL", "P_VALUE", "PVALUE")] = "P"

  ## Check that the required columns are present
  required = c("chr", "pos", "a0", "a1", "beta", "beta_se", "n_eff")
  missing_cols = required[!required %in% colnames(df)]
  if(length(missing_cols) > 0){
    stop(paste("The GWAS is missing columns that could not be inferred:",
               paste(missing_cols, collapse = ", "),
               "- please rename them before calling endoPRS."))
  }

  df$chr = as.numeric(df$chr)

  ## REGENIE and some other tools report very small p-values as NA, so recompute them from the
  ## chi-square statistic rather than dropping those variants (they are the most significant ones).
  if(!"P" %in% colnames(df)){
    df$P = NA_real_
  }
  needs_p = is.na(df$P)
  if(any(needs_p)){
    df$P[needs_p] = pchisq((df$beta[needs_p] / df$beta_se[needs_p])^2, df = 1, lower.tail = FALSE)
  }

  ## Drop variants that cannot be used by the model
  df = df[!is.na(df$beta) & !is.na(df$beta_se) & df$beta_se > 0 & !is.na(df$n_eff), ]

  ## Optionally restrict to hapmap variants
  if(filter_hapmap){
    if(is.null(hapmap)){
      stop("Must provide hapmap data frame if selected filter_hapmap option.")
    }
    colnames(hapmap)[toupper(colnames(hapmap)) %in% c("CHROM", "CHR")] = "CHR"
    colnames(hapmap)[toupper(colnames(hapmap)) %in% c("POS", "BP", "GENPOS")] = "POS"
    hapmap_pos = paste(hapmap$CHR, hapmap$POS, sep = "_")
    df = df[paste(df$chr, df$pos, sep = "_") %in% hapmap_pos, ]
  }

  ## Align to the LD reference. snp_match flips the sign of beta where the alleles are reversed,
  ## which is why the correlation matrix does not need to be sign flipped later on.
  df_matched = snp_match(df[, c(required, "P")], ld_map, join_by_pos = TRUE)

  ## build_corr_ss relies on the rows being in LD reference order
  df_matched = df_matched[order(df_matched$`_NUM_ID_`), ]

  return(df_matched)
}


#' Function to assign SNPs to the phenotype only, endophenotype only, and shared groups
#'
#' This function is the summary statistic analogue of \code{extract_snp_groups}. It differs in two
#' ways. First, SNPs are keyed by their row in the LD reference panel rather than by a string ID,
#' which removes any ambiguity about allele order. Second, the phenotype and the endophenotype are
#' allowed to use different p-value thresholds, which is the extension introduced in endoPRS-SS.
#'
#' @param pheno_gwas_matched Output of format_gwas_ss applied to the phenotype GWAS.
#' @param endo_gwas_matched Output of format_gwas_ss applied to the endophenotype GWAS.
#' @param thresh P-value threshold used to determine association with the phenotype.
#' @param thresh_endo P-value threshold used to determine association with the endophenotype. If NULL, thresh is used for both, which reproduces the behavior of the individual level version of endoPRS.
#'
#' @return A data frame with one row per selected variant, containing the column `_NUM_ID_` (its row
#'   in the LD reference) and group (pheno_only, endo_only, or both), sorted by `_NUM_ID_`.
#'
extract_snp_groups_ss = function(pheno_gwas_matched, endo_gwas_matched, thresh, thresh_endo = NULL){

  ## Default to a single threshold for both traits
  if(is.null(thresh_endo)) thresh_endo = thresh

  ## Variants passing each threshold, identified by their row in the LD reference
  pheno_hits = pheno_gwas_matched$`_NUM_ID_`[pheno_gwas_matched$P < thresh]
  endo_hits  = endo_gwas_matched$`_NUM_ID_`[endo_gwas_matched$P < thresh_endo]

  if(length(pheno_hits) == 0 & length(endo_hits) == 0){
    stop(paste("No variants pass the thresholds thresh =", thresh,
               "and thresh_endo =", thresh_endo, "- please loosen them."))
  }

  ## Effect sizes are always taken from the phenotype GWAS, including for endophenotype only
  ## variants, so a variant can only be used if the phenotype GWAS reports it.
  all_hits = sort(union(pheno_hits, endo_hits))
  usable = all_hits %in% pheno_gwas_matched$`_NUM_ID_`
  if(any(!usable)){
    warning(paste(sum(!usable), "of", length(all_hits),
                  "selected variants are absent from the phenotype GWAS and were dropped.",
                  "endoPRS-SS takes effect sizes from the phenotype GWAS for every variant,",
                  "including endophenotype only variants."))
    all_hits = all_hits[usable]
  }

  ## Label each variant by which traits it is associated with
  snp_groups = data.frame(`_NUM_ID_` = all_hits, group = "endo_only",
                          check.names = FALSE, stringsAsFactors = FALSE)
  snp_groups$group[all_hits %in% pheno_hits] = "pheno_only"
  snp_groups$group[all_hits %in% pheno_hits & all_hits %in% endo_hits] = "both"

  return(snp_groups)
}


#' Function to create the vector of per-variant penalties for endoPRS-SS
#'
#' This function is the summary statistic analogue of \code{create_penalty_table}. SNPs associated
#' with only the phenotype are given a penalty of 1, SNPs associated with only the endophenotype are
#' given a penalty of w2, and SNPs associated with both are given a penalty of w3.
#'
#' @param snp_groups Output of extract_snp_groups_ss.
#' @param df_beta The data frame of summary statistics the model will be fit on. Used to put the penalties in the same order as the rows of df_beta.
#' @param w2 Weight to assign SNPs in the endophenotype only group.
#' @param w3 Weight to assign SNPs associated with both the phenotype and the endophenotype.
#'
#' @return A data frame with one row per row of df_beta, containing `_NUM_ID_`, group, and penalty.
#'
create_penalty_table_ss = function(snp_groups, df_beta, w2, w3){

  ## Put the groups in the same order as the rows of df_beta and the columns of corr
  group = snp_groups$group[match(df_beta$`_NUM_ID_`, snp_groups$`_NUM_ID_`)]

  if(any(is.na(group))){
    stop("Some variants in df_beta were not assigned to a SNP group. Please check your input.")
  }

  penalty_table = data.frame(`_NUM_ID_` = df_beta$`_NUM_ID_`, group = group,
                             penalty = NA_real_, check.names = FALSE)
  penalty_table$penalty[group == "pheno_only"] = 1
  penalty_table$penalty[group == "endo_only"] = w2
  penalty_table$penalty[group == "both"] = w3

  ## The penalty is a multiplicative factor on lambda, so it has to be strictly positive
  if(any(is.na(penalty_table$penalty)) | any(penalty_table$penalty <= 0)){
    stop("All weights must be strictly positive and every SNP must be assigned one.")
  }

  return(penalty_table)
}
