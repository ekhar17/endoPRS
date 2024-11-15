% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_endoPRS.R
\name{fit_endoPRS}
\alias{fit_endoPRS}
\title{Function to apply endoPRS method to generate a polygenic risk score model.}
\usage{
fit_endoPRS(
  G,
  map,
  fam,
  train_pheno,
  train_covar,
  val_pheno,
  val_covar,
  pheno_gwas,
  endo_gwas,
  filter_hapmap = F,
  hapmap = NULL,
  type = NULL,
  threshes = c(0.01, 1e-04, 1e-06),
  grid = NULL,
  NCORES = NULL,
  pheno_gwas_refit,
  endo_gwas_refit,
  save_folder = NULL
)
}
\arguments{
\item{G}{An object of class FBM from the bigsnpr package that contains genotypes of the individuals used for training and those used for validation.}

\item{map}{A data frame containing information about SNPs in G (genotype matrix). This can be the map object form the bigSNP class. It should contain columns labeled chromosome, marker ID, and allele1 and allele2.}

\item{fam}{A data frame containing information about the individuals in G (genotype matrix). Column one must correspond to FID and column 2 must correspond to IID.}

\item{train_pheno}{A data frame with three columns corresponding to the phenotypes of individuals used for training. Column one is FID, column two is IID, and column three is y the phenotype. If the phenotype is binary, y must consist of 0's and 1's.}

\item{train_covar}{A data frame corresponding with the covariates of the individuals used for training (ie genetic PC's, sex, age, etc). The first column must be FID, and the second column must be IID. The order of individuals must be the same as train_pheno. Covariates that are factors (i.e. assessment center) must be encoded as such. The covariates are not penalized in the endoPRS model.}

\item{val_pheno}{A data frame with three columns corresponding to the phenotypes of individuals used for validation. Column one is FID, column two is IID, and column three is y the phenotype. If the phenotype is binary, y must consist of 0's and 1's.}

\item{val_covar}{A data frame corresponding with the covariates of the individuals used for validation (ie genetic PC's, sex, age, etc). The first column must be FID, and the second column must be IID. The order of individuals must be the same as val_pheno. Covariates that are factors (i.e. assessment center) must be encoded as such. The covariates are not penalized in the endoPRS model.}

\item{pheno_gwas}{A data frame containing the results of the GWAS run on the phenotype. It must not include any individuals from the validation set. It must contain columns corresponding to chromosome, SNP, allele1, allele2, and P-value.}

\item{endo_gwas}{A data frame containing the results of the GWAS run on the endophenotype. It must not include any individuals from the validation set. It must contain columns corresponding to chromosome, SNP, allele1, allele2, and P-value.}

\item{filter_hapmap}{An optional logical character. It corresponds to whether to only run the model on variants in the hapmap3 set.}

\item{hapmap}{An optional data frame with hapmap variants. This must be provided if filter_hapmap is set to TRUE. The hapmap3 data frame should contain columns CHR, SNP, A1, A2.}

\item{type}{An optional character vector of "linear" or "logistic." It specifies which type of model to fit to the data. If not provided, it will be learned from the phenotype. If the phenotype consists of 0's and 1's penalized logistic regression will be used, otherwie penalized linear regression will be used.}

\item{threshes}{An optional character of p-value thresholds. This corresponds to the GWAS p-value thresholds used to determine association with the phenotype and endophenotype. The optimal p-value threshold is learned from tuning on the validation set. If not provided, the original values from the endoPRS manuscript will be used (0.01, 1e-4, and 1e-6).}

\item{grid}{An optional data frame of weights to use for the weighted penalty in the lasso model. The first column must be w2 and correspond to the weights applied to the group of SNPs associated with only the endophenotype and the second column must be w3 and correspond to the weights applied to the group of SNPs associated with both the phenotype and endophenotype. The set of SNPs associated with only the phenotype is given a weight of 1. If not provided, the original grid of weights from the endoPRS manuscript (0.1, 0.5, 1, 2, 10) will be used for both w2 and w3.}

\item{NCORES}{An optional value corresponding to the number of cores to use. Otherwise, the number of cores will be learned using nb_cores().}

\item{pheno_gwas_refit}{An optional data frame containing the results of the GWAS run on the phenotype using both the training and validation set. It must contain columns corresponding to chromosome, SNP, allele1, allele2, and P-value.}

\item{endo_gwas_refit}{An optional data frame containing the results of the GWAS run on the endophenotype using both the training and validation set. It must contain columns corresponding to chromosome, SNP, allele1, allele2, and P-value.}

\item{save_folder}{An optional argument. If included, it should specify a path to a directory that files can be written to. If specified, all the models that are fit are saved to that directory. This is particularly useful for larger data sets, as model fitting can take a while, so saving intermediate results can prevent the need for rerunning the same models again.}
}
\value{
A list with two elements: the betas of the final model and the model itself
\itemize{
\item{beta: A data frame consisting of the SNPs included in the final model and their corresponding coefficients. This can be used to apply the PRS to an external data set using software such as PLINK.}
\item{model: The final fitted endoPRS model as an object of class big_sp_list. This can be used to apply the final model to a test set using the predict() function from the bigsnpr package.  }
}
}
\description{
This function applies the endoPRS model. It runsa weighted lasso model that penalizes the SNPs differently based on
whether they are  associated with only the phenotype, only the endophenotype, or both. The optimal set of weights is
determined by validation set performance. The model is then refit with the chosen weights on the combined training and
validation set to obtain the final PRS model.
}
