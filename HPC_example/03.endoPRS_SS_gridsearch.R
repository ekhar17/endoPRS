## Example of running one cell of the endoPRS-SS tuning grid on a cluster.
## Each job fits one (thresh, thresh_endo, weight) combination and writes its results to save_folder.

library(bigsnpr)
library(data.table)
library(endoPRS)

args = commandArgs(trailingOnly = TRUE)
thresh      = as.numeric(args[1])
thresh_endo = as.numeric(args[2])
iteration   = as.numeric(args[3])
save_folder = args[4]

## Genotypes of the validation individuals
obj.bigsnp = snp_attach("/path/to/genotypes.rds")
G   = obj.bigsnp$genotypes
map = obj.bigsnp$map
fam = obj.bigsnp$fam

## Validation phenotype and covariates - FID, IID, then the phenotype / covariates
val_pheno = as.data.frame(fread("/path/to/val_pheno.txt"))
val_covar = as.data.frame(fread("/path/to/val_covar.txt"))

## GWAS summary statistics. They must not include the validation individuals.
pheno_gwas = as.data.frame(fread("/path/to/pheno_gwas.txt"))
endo_gwas  = as.data.frame(fread("/path/to/endo_gwas.txt"))

## LD reference panel - one correlation matrix and one variant file per chromosome
ld_ref = load_ld_ref("/path/to/ld_reference",
                     corr_pattern = "chr%d.rds",
                     info_pattern = "chr%d_info.csv")

fit_endoPRS_ss_single_iter(G, map, fam,
                           val_pheno, val_covar,
                           pheno_gwas, endo_gwas,
                           ld_ref = ld_ref,
                           thresh = thresh,
                           thresh_endo = thresh_endo,
                           iteration = iteration,
                           save_folder = save_folder)
