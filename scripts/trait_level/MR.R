# Perform MR analysis
# Load packages
library(bigreadr)
library(plyr)
library(dplyr)
library(ggplot2)

library(TwoSampleMR)
library(ieugwasr)
library(plinkbinr)
library(MRPRESSO)

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
REF_PANEL <- "/gwashug/reference/1kg"
WIN_SIZE <- 10000
R2 <- 0.001
SIG_CUTOFF <- 1E-8

# Function 1: process summary_info.txt
summaryInfo.process <- function(summaryInfo_file
){
  
  summary_info <- read.table(summaryInfo_file, sep = "\t",quote="")
  summary_list <- list(summ1 = summary_info[1, 2], 
                       pop1 = summary_info[2, 2],
                       trait1 = summary_info[3, 2],
                       summ2 = summary_info[4, 2], 
                       pop2 = summary_info[5, 2],
                       trait2 = summary_info[6, 2],
                       MRp = as.numeric(summary_info[7, 2]))
  return(summary_list)
}

# Function 2: clump exposure data
exp.clump <- function(exp_data, 
                      pop,
                      MRP
){
  
  ## format exposure data
  exp_data <- exp_data %>%
    filter(p_wald < MRP)  %>%
    format_data(.,
                type='exposure',
                snp_col = "rs",
                beta_col = "beta",
                se_col = "se",
                eaf_col = "af",
                effect_allele_col = "allele1",
                other_allele_col = "allele0",
                pval_col = "p_wald")
  exp_data$id <- exp_data$id.exposure 
  exp_data$rsid <- exp_data$SNP
  exp_data$pval <- exp_data$pval.exposure
  
  ## clumping
  ref <- paste0(REF_PANEL, "/", pop, "/hm3_imp/merge")
  exp_clumped <- ld_clump(exp_data,
                          plink_bin = get_plink_exe(),
                          bfile = ref,
                          clump_kb = WIN_SIZE, 
                          clump_r2 = R2)
  return(exp_clumped)
}

# Function 3: process MR
MR.process <- function(exp_clumped,
                       out_file,
                       trait1,
                       trait2
){
  
  exp_clumped$id.exposure <- trait1
  out_dat <- read_outcome_data(snps = exp_clumped$SNP,
                               filename = out_file,
                               sep = " ",
                               snp_col = "rs",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "af",
                               effect_allele_col = "allele1",
                               other_allele_col = "allele0",
                               pval_col = "p_wald")
  
  out_dat$id.outcome <- trait2
  mr_data <- harmonise_data(exposure_dat = exp_clumped, 
                            outcome_dat = out_dat,
                            action= 3)
  mr_data <- mr_data %>% filter(pval.outcome > SIG_CUTOFF)
  return(mr_data)
}

## Function 4: forest plot
mr.forestplot <- function(single_result, 
                          trait1, 
                          trait2
){
  
  single_result$SNP[nrow(single_result)] <- "All-IVW"
  plt <- mr_forest_plot(single_result)
  plt <- plt[[1]] + theme_bw() +
    labs(x = paste0("MR Effect Size for ", trait1, " on ", trait2)) +
    guides(size = "none", colour = "none") +
    theme(legend.position = "none",
          plot.tag = element_text(size = 25, face = "bold"),
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 10), 
          axis.text.x = element_text(size = 8), 
          axis.title.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 8), 
          axis.title.y = element_text(size = 8, face = "bold"))
  return(plt)
}

## Function 5: funnel plot
mr.funnelplot <- function(single_result, 
                          trait1, 
                          trait2
){
  
  plt <- mr_funnel_plot(single_result)
  plt <- plt[[1]] + theme_bw() +
    labs(color = "MR Method") +
    scale_color_discrete(labels = c("Inverse variance weighted" = "IVW")) +
    theme(legend.position = "right",
          plot.tag = element_text(size = 25, face = "bold"),
          legend.title = element_text(size = 18, face = "bold"),
          legend.text = element_text(size = 15), 
          axis.text.x = element_text(size = 15), 
          axis.title.x = element_text(size = 18, face = "bold"),
          axis.text.y = element_text(size = 15), 
          axis.title.y = element_text(size = 18, face = "bold"))
  return(plt)
}

## Function 6: scatter plot
mr.scatterplot <- function(mr_result,
                           mr_dat, 
                           trait1, 
                           trait2
){
  
  mr_result <- mr_result[1:3, ]
  mr_result$method[3] <- "IVW"
  plt <- mr_scatter_plot(mr_result, mr_dat)
  plt <- plt[[1]] + theme_bw() +
    labs(x = paste0("SNP Effect on ", trait1),
         y = paste0("SNP Effect on ", trait2),
         color = "MR Method") +
    guides(colour = ggplot2::guide_legend(ncol = 1)) +
    theme(legend.position = "right",
          legend.box = "vertical",
          plot.tag = element_text(size = 25, face = "bold"),
          legend.title = element_text(size = 18, face = "bold"),
          legend.text = element_text(size = 10), 
          axis.text.x = element_text(size = 10), 
          axis.title.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10), 
          axis.title.y = element_text(size = 10, face = "bold"))
  return(plt)
}

## Function 6: LOO plot
mr.looplot <- function(res_loo,
                       trait1, 
                       trait2
){
  
  plt <- mr_leaveoneout_plot(res_loo)
  plt <- plt[[1]] + theme_bw() +
    labs(x = paste0("Sensitivity Analysis for ", trait1, 
                    " on ", trait2)) +
    guides(size = "none", colour = "none") +
    theme(legend.position = "none",
          plot.tag = element_text(size = 25, face = "bold"),
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 10), 
          axis.text.x = element_text(size = 10), 
          axis.title.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10), 
          axis.title.y = element_text(size = 10, face = "bold"))
  return(plt)
}

# Function 7: call MR
MR.call <- function(summaryInfo_file
){
  
  # summaryInfo_file = "/gwashug/data/binary_binary/UP138_UP457/MR/summary_info.txt"
  
  # summaryInfo_file=opt$summaryInfo
  summary_file <- summaryInfo.process(summaryInfo_file)
  ## forward MR
  ### MR analysis
  exp_dat <- fread2(summary_file[[1]])
  exp_clump <- exp.clump(exp_dat,
                         summary_file[[2]],
                         summary_file[[7]])
  MR_dat <- MR.process(exp_clump, summary_file[[4]], 
                       summary_file[[3]], summary_file[[6]])
  MR_result <- mr(MR_dat, method_list = c("mr_egger_regression", 
                                          "mr_weighted_median",
                                          "mr_ivw_mre", 
                                          "mr_ivw_fe"))
  MR_result$or <- exp(MR_result$b)
  MR_result$or_lci95 <- exp(MR_result$b - 1.96*MR_result$se)
  MR_result$or_uci95 <- exp(MR_result$b + 1.96*MR_result$se)
  ### sensitivity analysis
  hetero_result <- mr_heterogeneity(MR_dat)
  pleio_result <- mr_pleiotropy_test(MR_dat)
  single_result <- mr_singlesnp(MR_dat, all_method = c("mr_ivw"))
  loo_result <- mr_leaveoneout(MR_dat)
  ### Method selection
  if(pleio_result$pval <= 0.05){
    
    MR_result$best <- c("*", "-", "-", "-")
  } else {
    
    
    if(hetero_result$Q_pval[2] <= 0.05){
      
      MR_result$best <- c("-", "*", "*", "-")
    } else {
      
      MR_result$best <- c("-", "-", "-", "*")
    }
  }
  
  ## Output table
  out_path <- gsub("summary_info.txt", "", summaryInfo_file)
  system(paste0("mkdir ", out_path, "/output"))
  write.csv(MR_result,file = paste0(out_path, "/output/MR_result.csv"), 
            row.names = F, quote = F)
  write.csv(hetero_result, file = paste0(out_path, "/output/heterogeneity.csv"), 
            row.names = F, quote = F)
  write.csv(pleio_result, file = paste0(out_path, "/output/pleiotropy.csv"), 
            row.names = F, quote = F)
  write.csv(loo_result, file = paste0(out_path, "/output/loo.csv"), 
            row.names = F, quote = F)
  
  ## Plot
  forestplot <- mr.forestplot(single_result, summary_file[[3]], summary_file[[6]])
  funnelplot <- mr.funnelplot(single_result, summary_file[[3]], summary_file[[6]])
  scatterplot <- mr.scatterplot(MR_result, MR_dat, summary_file[[3]], summary_file[[6]])
  looplot <- mr.looplot(loo_result, summary_file[[3]], summary_file[[6]])
  ggsave(paste0(out_path, "/output/forestplot.png"), forestplot, 
         width = 7, height = 8, dpi = 300)
  ggsave(paste0(out_path, "/output/funnelplot.png"), funnelplot,
         width = 7, height = 8, dpi = 300)
  ggsave(paste0(out_path, "/output/scatterplot.png"), scatterplot,
         width = 7, height = 8, dpi = 300)
  ggsave(paste0(out_path, "/output/looplot.png"), looplot,
         width = 7, height = 8, dpi = 300)
  
  ## Reverse MR
  exp_dat <- fread2(summary_file[[4]])
  exp_clump <- exp.clump(exp_dat, 
                         summary_file[[5]],
                         summary_file[[7]])
  MR_dat <- MR.process(exp_clump, summary_file[[1]], 
                       summary_file[[6]], summary_file[[3]])
  MR_result <- mr(MR_dat)
  MR_result$or <- exp(MR_result$b)
  MR_result$or_lci95 <- exp(MR_result$b - 1.96*MR_result$se)
  MR_result$or_uci95 <- exp(MR_result$b + 1.96*MR_result$se)
  write.csv(MR_result,file = paste0(out_path, "/output/reverse_MR_result.csv"), 
            row.names = F, quote = F)
}

# Run
MR.call(opt$info)
