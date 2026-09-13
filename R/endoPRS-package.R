#' @keywords internal
"_PACKAGE"

## endoPRS calls a number of functions from the bigstatsr / bigsnpr family without qualifying them.
## They are imported here so that the package works whether or not the user has attached bigsnpr.

#' @importFrom bigsnpr snp_match
#' @importFrom bigstatsr big_spLinReg big_spLogReg big_prodMat covar_from_df AUC seq_log
#'   cols_along rows_along
#' @importFrom bigsparser as_SFBM
#' @importFrom bigassertr assert_df_with_names assert_lengths assert_pos
#' @importFrom bigparallelr nb_cores assert_cores register_parallel
#' @importFrom foreach foreach %dopar%
#' @importFrom stats pchisq cor lm residuals sd
#' @importFrom utils txtProgressBar setTxtProgressBar read.csv write.csv
#' @importFrom methods is as
NULL
