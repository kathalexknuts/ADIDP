#some example code for the paper entitled "Integrating Brain Imaging Endophenotypes with GWAS for Alzheimer’s Disease"
#requires version R/3.5.0 or higher
library("dplyr"); library("aSPU"); library("gsmr"); library("tidyr"); library("devtools"); library("TwoSampleMR"); library("MRInstruments")

#set p-value threshold
thresh <- 5*10^(-5)

#read in IGAP data
#IGAP summary statistics can be obtained at http://web.pasteur-lille.fr/en/recherche/u744/igap/igap_download.php
AD_GWAS <- readRDS("path-to-IGAP-Summary-Statistics-File")
AD_GWAS$Phenotype <- "AD"

#format outcome data for downstream use in MRBase (Hemani et al.) functions, http://www.mrbase.org
outcome_dat <- format_data(AD_GWAS, header = TRUE, snp_col = "MarkerName", beta_col = "Beta", se_col = "SE", effect_allele_col = "Effect_allele", 
                           other_allele_col = "Non_Effect_allele",pval_col = "Pvalue",type = "outcome",chr_col = "Chromosome",pos_col = "Position")

#read in exposure GWAS of interest
#compressed UKBB summary statistic files for IDPs can be downloaded using the command line prompts given at:
#https://www.dropbox.com/s/qhiftre33pi70xs/BIG_summary_stats_files.xls?dl=0
#these files will likely be too large to read into R, so the LD clumping and thresholding performed
#below should instead be done using plink with individual level data from an external reference panel (i.e. 1000G)
GWAS <- na.omit(read.table("path-to-exposure-GWAS", header = FALSE, fill = TRUE))
names(GWAS) <- c("CHR", "SNP", "BP", "A1", "A2", "MAF", "BETA", "SE", "P")
GWAS$Phenotype <- "IDPName"

#format exposure data for downstream use in MRBase functions
exposure_dat <- format_data(GWAS, header = TRUE, snp_col = "SNP", beta_col = "BETA",se_col = "SE",effect_allele_col = "A1",other_allele_col = "A2",
                            pval_col = "P",phenotype_col="Phenotype",type = "exposure",chr_col = "CHR",pos_col = "BP")

#MRBase function for LD clumping (uses the same procedure as plink with the EUR 1000G panel, n = 503)
exposure_dat_clumped <- clump_data(exposure_dat, clump_kb = 1000, clump_r2 = 0.1)
exposure_dat_clumped <- exposure_dat_clumped %>% filter(pval.exposure <= thresh)

exposure_dat_clumped_2SMR <- clump_data(exposure_dat, clump_kb = 1000, clump_r2 = 0.001)
exposure_dat_clumped_2SMR <- exposure_dat_clumped %>% filter(pval.exposure <= thresh)

#combine outcome and exposure data- this data can then be used directly in MRBase 2-SMR functions
dat <- harmonise_data(exposure_dat = exposure_dat_clumped, outcome_dat = outcome_dat)
dat_2SMR <- harmonise_data(exposure_dat = exposure_dat_clumped_2SMR, outcome_dat = outcome_dat)

#These give the same estimates as those using plink + the 1000G EUR data 
LD <- ld_matrix(dat$SNP, with_alleles = FALSE) 

#run all available MR meta analysis approaches for the given exposure phenotype
SMR <- mr(dat_2SMR)

#MR results in wide format, generally a nicer format to use for plots etc. 
mr_wide <- SMR %>%
  gather(variable, value, -c(id.exposure:nsnp)) %>%
  unite(temp, method, variable) %>%
  spread(temp, value)

#SPU/aSPU package maintained by Il-Youp Kwak
#make sure all vectors/matricies below have the same SNP order before multiplying
Z <- as.matrix(as.numeric(dat$beta.outcome)/as.numeric(dat$se.outcome), nrow = nrow(dat))
Z_mat <- diag(length(Z)); diag(Z_mat) <- Z
beta <- matrix(dat$beta.exposure, nrow = nrow(dat)); beta_mat <- diag(length(beta)); diag(beta_mat) <- beta; rownames(beta_mat) <- dat$SNP
table(rownames(LD) == rownames(beta_mat)) #should all be TRUE
Tsum <- t(beta)%*%Z; VarTsum <- t(beta)%*%LD%*%beta; #Tsum should be equal to the SPUs1 estimate in the output from the aSPUs function below
SPU <- aSPUs(beta_mat %*% Z, beta_mat %*%LD%*%beta_mat, n.perm = 1000, Ps = FALSE, prune = FALSE) 


#Run GSMR using the packages by Zhihong Zhu, Zhili Zheng, Futao Zhang and Jian Yang
bzx = dat$beta.exposure; bzx_se = dat$se.exposure; bzx_pval = dat$pval.exposure # SNP effects on the risk factor
bzy = dat$beta.outcome; bzy_se = dat$se.outcome; bzy_pval = dat$pval.outcome    # SNP effects on the disease
gsmr_results =  gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, LD, as.character(dat$SNP), 503, F, 1, 0.01, 0.01, 10, 0.1, 0.05, 0)


