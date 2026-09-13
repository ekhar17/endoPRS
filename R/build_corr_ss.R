#' Function to load an LD reference panel stored as one file per chromosome
#'
#' This is a convenience function for the common case where the LD reference panel is stored in the
#' LDpred2 format: one correlation matrix per chromosome saved as an .rds file, plus one file per
#' chromosome describing the variants in that matrix. The variant files must list the variants in
#' exactly the same order as the rows and columns of the corresponding correlation matrix.
#'
#' @param ld_dir Path to the directory containing the LD reference panel.
#' @param corr_pattern A sprintf pattern for the correlation matrix file names, where the chromosome is substituted in. Default is "chr\%d.rds".
#' @param info_pattern A sprintf pattern for the variant information file names. Default is "chr\%d_info.csv".
#' @param chrs Which chromosomes to load. Default is 1 to 22.
#'
#' @return A list with two elements, suitable for the ld_ref argument of fit_endoPRS_ss:
#' \itemize{
#' \item{map: A data frame of all reference variants with columns chr, pos, a0, a1, stacked in chromosome order.}
#' \item{corr_files: A named character vector of paths to the per-chromosome correlation matrices.}
#' }
#'
#' @export
load_ld_ref = function(ld_dir,
                       corr_pattern = "chr%d.rds",
                       info_pattern = "chr%d_info.csv",
                       chrs = 1:22){

  ld_map = data.frame()
  corr_files = character()

  for(chr in chrs){

    info_file = file.path(ld_dir, sprintf(info_pattern, chr))
    corr_file = file.path(ld_dir, sprintf(corr_pattern, chr))

    if(!file.exists(info_file)) stop(paste("Cannot find LD reference info file:", info_file))
    if(!file.exists(corr_file)) stop(paste("Cannot find LD reference matrix file:", corr_file))

    info = as.data.frame(data.table::fread(info_file))

    ## Accept either the LDpred2 reference naming or the bigSNP map naming. Note that bigsnpr's
    ## convention is allele1 = a1 (the effect allele) and allele2 = a0.
    colnames(info)[toupper(colnames(info)) %in% c("CHROMOSOME", "CHROM", "CHR")] = "chr"
    colnames(info)[toupper(colnames(info)) %in% c("PHYSICAL.POS", "POS", "BP")] = "pos"
    if(!"a1" %in% colnames(info)) colnames(info)[colnames(info) == "allele1"] = "a1"
    if(!"a0" %in% colnames(info)) colnames(info)[colnames(info) == "allele2"] = "a0"

    missing_cols = c("chr", "pos", "a0", "a1")[!c("chr", "pos", "a0", "a1") %in% colnames(info)]
    if(length(missing_cols) > 0){
      stop(paste("LD reference info file", info_file, "is missing columns:",
                 paste(missing_cols, collapse = ", ")))
    }

    info$chr = as.numeric(info$chr)
    ld_map = rbind(ld_map, info[, c("chr", "pos", "a0", "a1")])
    corr_files[as.character(chr)] = corr_file
  }

  return(list(map = ld_map, corr_files = corr_files))
}


#' Function to build the LD correlation matrix restricted to the selected variants
#'
#' This function subsets the per-chromosome LD reference matrices down to the variants selected by
#' the endoPRS p-value thresholds, and stacks them into a single block diagonal SFBM that can be
#' passed to weighted_lassosum2.
#'
#' Because format_gwas_ss aligns the summary statistics to the LD reference with
#' \code{bigsnpr::snp_match}, the effect sizes are already on the reference allele, so no sign
#' flipping of the correlation matrix is needed here.
#'
#' @param df_beta A data frame of summary statistics restricted to the selected variants, as produced by format_gwas_ss. It must contain the column `_NUM_ID_` and be sorted by it.
#' @param ld_ref A list with elements map and corr_files, as produced by load_ld_ref.
#' @param backingfile Path (without extension) for the backing file of the resulting SFBM. Defaults to a temporary file.
#'
#' @return An object of class SFBM whose rows and columns correspond, in order, to the rows of df_beta.
#'
build_corr_ss = function(df_beta, ld_ref, backingfile = tempfile()){

  ld_map = ld_ref$map
  corr_files = ld_ref$corr_files

  ## The columns of corr must line up with the rows of df_beta
  if(is.unsorted(df_beta$`_NUM_ID_`, strictly = TRUE)){
    stop("df_beta must be sorted by `_NUM_ID_` with no duplicates before building the LD matrix.")
  }

  corr = NULL

  for(chr in sort(unique(df_beta$chr))){

    ## Rows of df_beta, and rows of the reference panel, that are on this chromosome
    ind_chr = which(df_beta$chr == chr)
    ld_rows_chr = which(ld_map$chr == chr)

    ## Position of each selected variant within this chromosome's correlation matrix
    ind_local = match(df_beta$`_NUM_ID_`[ind_chr], ld_rows_chr)
    if(any(is.na(ind_local))){
      stop(paste("Some selected variants on chromosome", chr,
                 "could not be located in the LD reference panel. Please check your input."))
    }

    if(is.na(corr_files[as.character(chr)])){
      stop(paste("No LD reference matrix was supplied for chromosome", chr))
    }

    corr_chr = readRDS(corr_files[as.character(chr)])

    if(nrow(corr_chr) != length(ld_rows_chr) | ncol(corr_chr) != length(ld_rows_chr)){
      stop(paste("The LD matrix for chromosome", chr, "has", nrow(corr_chr),
                 "rows but the reference panel lists", length(ld_rows_chr),
                 "variants on that chromosome. The variant file must be in the same order",
                 "as the matrix."))
    }

    corr_chr = corr_chr[ind_local, ind_local, drop = FALSE]

    ## A single variant on a chromosome subsets down to a dense 1x1 matrix
    if(!methods::is(corr_chr, "sparseMatrix")){
      corr_chr = methods::as(corr_chr, "sparseMatrix")
    }

    if(is.null(corr)){
      corr = as_SFBM(corr_chr, backingfile, compact = TRUE)
    } else {
      corr$add_columns(corr_chr, nrow(corr))
    }
  }

  if(corr$ncol != nrow(df_beta)){
    stop("The LD matrix and the summary statistics ended up with different numbers of variants.")
  }

  return(corr)
}
