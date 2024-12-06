# Estimate genetic correlation
# Load packages
library(bigreadr)
library(plyr)
library(dplyr)
library(ggplot2)

library(GECKO)

library(optparse)

# Input parameters
args_list = list(
  make_option("--info", type="character", default=NULL,
              help="INPUT: summary information",
              metavar="character")
)
opt_parser <- OptionParser(option_list = args_list)
opt <- parse_args(opt_parser)

# Fix parameter
REF_PANEL <- "/gwashug/reference/"
TEST = T
WIGHTIN = T 

# Function 1: process summary_info.txt
summaryInfo.process <- function(summaryInfo_file
){
  
  summary_info <- read.table(summaryInfo_file, sep = "\t")
  summary_list <- list(summ1 = summary_info[1, 2], 
                       pop1 = summary_info[2, 2],
                       summ2 = summary_info[3, 2], 
                       pop2 = summary_info[4, 2],
                       same_study = summary_info[5, 2])
  return(summary_list)
}

# Function 2: transfer to gemma format
GECKO.call <- function(summaryInfo_file
){
  
  ## process summary statistics
  summary_file <- summaryInfo.process(summaryInfo_file)
  summ_in1 <- fread2(summary_file$summ1)
  summ_out1 <- data.frame(chr = summ_in1$chr,
                          bp = summ_in1$ps, 
                          SNP = summ_in1$rs, 
                          A1 = summ_in1$allele0, 
                          A2 = summ_in1$allele1, 
                          N = summ_in1$n_obs, 
                          Z = summ_in1$beta/summ_in1$se, 
                          P = summ_in1$p_wald)
  n_in1 <- round(mean(summ_out1$N))
  summ_in2 <- fread2(summary_file$summ2)
  summ_out2 <- data.frame(chr = summ_in2$chr,
                          bp = summ_in2$ps, 
                          SNP = summ_in2$rs, 
                          A1 = summ_in2$allele0, 
                          A2 = summ_in2$allele1, 
                          N = summ_in2$n_obs, 
                          Z = summ_in2$beta/summ_in2$se, 
                          P = summ_in2$p_wald)
  n_in2 <- round(mean(summ_out2$N))
  
  ## set parameters
  if (summary_file[[5]] == "T"){
    
    nsin <- min(n_in1, n_in2)
    Fix_Vein <- F
  } else {
    
    nsin <- 0
    Fix_Vein <- T 
  }
  
  ## set ldsc
  ldscore <- fread2(paste0(REF_PANEL, summary_file$pop1, 
                           "_w_ld_chr/ldscore.txt"))
  snp_inter <- intersect(intersect(summ_out1$SNP, summ_out2$SNP), 
                          ldscore$SNP)
  ldscore <- ldscore %>% filter(SNP %in% snp_inter)
  summ_out1 <- summ_out1 %>% filter(SNP %in% snp_inter)
  summ_out2 <- summ_out2 %>% filter(SNP %in% snp_inter)
  
  ## run GECKO
  herit <- GECKO_R(summ_out1, 
                   summ_out2, 
                   n_in1, 
                   n_in2, 
                   nsin, 
                   ldscore,
                   WIGHTIN,
                   Fix_Vein,
                   TEST)
  cat(herit)
  v_g <- matrix(data = c(herit[1,5],herit[1,6],herit[1,6],herit[1,7]),
                nrow = 2, ncol = 2, byrow = T)
  
  ## output
  out_path <- gsub("summary_info.txt", "", summaryInfo_file)
  dir.create(paste0(out_path, "/output/"))
  write.table(v_g, file = paste0(out_path, "/output/herit.txt"), 
              row.names = F, quote = F)
  
  return(0)
}


GECKO.call(opt$info)
