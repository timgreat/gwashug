# Perform coloc analysis
# Load packages
library(bigreadr)
library(dplyr)
library(glue)
library(ggplot2)
library(ggrepel)
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

SMR <- "/software/smr/smr-1.3.1"
REF_PATH <- "/gwashug/reference/1kg/"
GENE_PATH <- "/gwashug/reference/"
EQTL_PATH <- "/gwashug/reference/eqtl/GTEx_eGene/"
#EQTL_PATH <- "/gwashug/reference/eqtl/cell/"

WIN_SIZE <- 10000 # Kb
N_CORE <- 1

# Function 1: process summary_info.txt
summaryInfo.process <- function(summaryInfo_file
){
  
  summary_info <- read.table(summaryInfo_file, sep = "\t",quote="")
  summary_list <- list(summ1 = summary_info[1, 2], 
                       summ2 = summary_info[2, 2], 
                       use_model1 = summary_info[3, 2], 
                       use_model2 = summary_info[4, 2])
  return(summary_list)
}

# Function 2: transformat summ from gemma to SMR
gemma2SMR <- function(summ_gemma_path = NULL,
                      trait = NULL,
                      out_path = NULL){
  #
  summ_gemma <- fread2(summ_gemma_path)
  summ_smr <- summ_gemma[, c("rs", "allele1", "allele0", "af", 
                             "beta", "se", "p_wald", "n_obs")]
  colnames(summ_smr) <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
  summ_smr_path <- glue("{out_path}/summ_{trait}.ma")
  fwrite2(summ_smr, summ_smr_path, sep = "\t")
  
  return(summ_smr_path)
}

# Function 2: SMR test
smr.test <- function(summ_gemma_path = NULL,
                     eth = "EUR",
                     trait = NULL,
                     # MRPvalue = 5E-8,
                     eqtl_path = NULL,
                     eqtl_model = NULL,
                     out_path = NULL){
  
  ## transformat summ from gemma to SMR
  #out_path <- dirname(summ_gemma_path)
  summ_path <- gemma2SMR(summ_gemma_path = summ_gemma_path,
                         trait = trait,
                         out_path = out_path)
  print(summ_path)
  ## other parameters
  if (!file.exists(out_path)) {
    system(glue("mkdir -f {out_path}"))
  }
  ref_prefix <- paste0(REF_PATH,"/EUR/hm3_imp/merge")
  eqtl_model_s <- sub("\\GTEx_", "", eqtl_model)
  eqtl_prefix <- paste0(eqtl_path, "/", eqtl_model_s, "_hm3_cis")
  out_prefix <- paste0(out_path, "/SMR_", trait, "_", eqtl_model_s)
  ## run SMR
  
  cmd_smrx <- glue("{SMR} --bfile {ref_prefix} ",
                   "--gwas-summary {summ_path} ", 
                   "--eqtl-summary {eqtl_prefix} ",
                   "--peqtl-smr 5e-8 ",
                   "--ld-upper-limit 0.9 ",
                   "--ld-lower-limit 0.05 ",
                   "--heidi-mtd 1 ",
                   "--heidi-min-m 3 ",
                   "--heidi-max-m 20 ",
                   "--peqtl-heidi 1.57e-3 ",
                   "--diff-freq-prop 1 ",
                   "--cis-wind {WIN_SIZE} ",
                   "--out {out_prefix} ",
                   "--thread-num {N_CORE}")
  # print(ref_prefix)
  # print(cmd_smrx)
  system(cmd_smrx)
  
  ## format output
  smr_resultx <- fread2(paste0(out_prefix, ".smr"))
  smr_resultx <- smr_resultx[, c("Gene", "probeID", "ProbeChr", "Probe_bp",
                                 "topSNP", "topSNP_bp", "A1", "A2", "Freq",
                                 "b_GWAS", "se_GWAS", "p_GWAS",
                                 "b_eQTL", "se_eQTL", "p_eQTL",
                                 "b_SMR", "se_SMR", "p_SMR",
                                 "p_HEIDI", "nsnp_HEIDI")]
  colnames(smr_resultx) <- c("Gene", "Ensembl", "Chr", "START_gene",
                             "topSNP", "BP_SNP", "A1", "A2", "EAF",
                             "Beta_GWAS", "SE_GWAS", "P_GWAS",
                             "Beta_eQTL", "SE_eQTL", "P_eQTL",
                             "Beta_SMR", "SE_SMR", "P_SMR",
                             "P_HEIDI", "Nsnp_HEIDI")
  return(smr_resultx)
}

# Function 2: Manh.plt
Manh.smr.plt <- function(LOCI_info = NULL,
                         X_axis = NULL,
                         p_val1 = NULL,
                         p_val2 = NULL,
                         trait1 = NULL,
                         trait2 = NULL,
                         n_top = 10,
                         thresh = NULL
){
  #
  if (is.null(thresh)) {
    thresh1 = -log10(0.05/sum(!is.na(p_val1)))
    thresh2 = -log10(0.05/sum(!is.na(p_val2)))
  }
  LOCI_info$P1 <- p_val1
  LOCI_info$P2 <- p_val2
  LOCI_info <- subset(LOCI_info, rowSums(is.na(LOCI_info)) == 0)
  LOCI_info$CHR <- as.character(LOCI_info$CHR)
  LOCI_info_top <- subset(LOCI_info, -log10(LOCI_info$P1) > thresh1 &
                            -log10(LOCI_info$P2) > thresh2)
  if (nrow(LOCI_info_top) == 0) {
    
    LOCI_info_top1 <- dplyr::top_n(LOCI_info, n = n_top, wt = -LOCI_info$P1)
    LOCI_info_top1 <- subset(LOCI_info_top1, -log10(LOCI_info_top1$P1) >= thresh1)
    LOCI_info_top2 <- dplyr::top_n(LOCI_info, n = n_top, wt = -LOCI_info$P2)
    LOCI_info_top2 <- subset(LOCI_info_top2, -log10(LOCI_info_top2$P2) >= thresh2)
    LOCI_info_top <- rbind(LOCI_info_top1, LOCI_info_top2)
  }
  
  # Manhattan plot
  plt <- ggplot(LOCI_info) +
    geom_point(aes(x = BPcum, y = -log10(P1), color = CHR),
               alpha = 0.8, size = 1.5)+
    geom_point(aes(x = BPcum, y = log10(P2), color = CHR),
               alpha = 0.8, size = 1.5)+
    geom_point(data = LOCI_info_top, 
               aes(x = BPcum, y = -log10(P1)),
               color = "red", alpha = 0.8, size = 1.5)+
    geom_point(data = LOCI_info_top, 
               aes(x = BPcum, y = log10(P2)),
               color = "red", alpha = 0.8, size = 1.5)+
    scale_x_continuous(label = X_axis$CHR, 
                       breaks= X_axis$center, 
                       limits = range(LOCI_info$BPcum)) +
    scale_y_continuous(expand = c(0, 0), labels = abs) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 11)) +
    xlab("Chromosome") + ylab(bquote("-"~Log[10]~"(P value)"))+
    geom_hline(yintercept = c(thresh1, 0, -thresh2), 
               color = c('red', "grey10", 'red'), 
               size = 0.8, 
               linetype = c(2, 1, 2)) + 
    annotate("text",
             x = quantile(LOCI_info$BPcum, 0.1), y = 1.5,
             size = 5, colour = "grey10",
             label = paste0("SMR for ", trait1)) +
    annotate("text",
             x = quantile(LOCI_info$BPcum, 0.1), y = -1.5,
             size = 5, colour = "grey10",
             label = paste0("SMR for ", trait2)) +
    geom_label_repel(data = LOCI_info_top, 
                     aes(label = TERM, x = BPcum, y = -log10(P1)), 
                     segment.size = 0.25, size = 3)+
    geom_label_repel(data = LOCI_info_top, 
                     aes(label = TERM, x = BPcum, y = log10(P2)), 
                     segment.size = 0.25, size = 3)+
    theme_bw() +
    theme(legend.position="none",
          axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 12, face = "bold", color = "black"),
          panel.border = element_blank(),
          axis.line.y = element_line(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  
  #
  return(plt)
}

# Function 3: call SMR analysis
smr.call <- function(summaryInfo_file){
  
  # 
  summary_file <- summaryInfo.process(summaryInfo_file)
  out_path <- paste0(dirname(summaryInfo_file), "/output")
  if (!file.exists(out_path)) {
    
    system(glue("mkdir {out_path}"))
  }
  
  ## eqtl_path_list: eqtl_model  eqtl_path
  use_model1 <- summary_file$use_model1
  use_model2 <- summary_file$use_model2
  eqtl_path_list <- read.table(paste0(EQTL_PATH, "eqtl_model_list.txt"), header = T)
  # if (is.null(use_model)) {
  #   use_line <- 1: nrow(eqtl_path_list)
  # } else {
  use_line <- match(c(use_model1, use_model2), eqtl_path_list$eqtl_model)
  # }
  #@timi
  if (use_model1==use_model2){
      use_line <- match(c(use_model1), eqtl_path_list$eqtl_model)
  }
  #@timi
  smr_out <- lapply(use_line, function(x){
    #
    
    eqtl_pathx <- eqtl_path_list$eqtl_path[x]
    eqtl_modelx <- eqtl_path_list$eqtl_model[x]
    #
    smr_retultx1 <- smr.test(summ_gemma_path = summary_file$summ1,
                             eth = summary_file$pop1,
                             # trait = summary_file$trait1,
                             trait = "trait1",
                             # MRPvalue = summary_file$MRPvalue,
                             eqtl_path = eqtl_pathx,
                             eqtl_model = eqtl_modelx,
                             out_path = out_path)
    write.table(smr_retultx1, 
                file = paste0(out_path, "/smr_trait1_", eqtl_modelx, ".txt"),
                sep = "\t", row.names = F, quote = F)
    #
    smr_retultx2 <- smr.test(summ_gemma_path = summary_file$summ2,
                             eth = summary_file$pop2,
                             # trait = summary_file$trait2,
                             trait = "trait2",
                             # MRPvalue = summary_file$MRPvalue,
                             eqtl_path = eqtl_pathx,
                             eqtl_model = eqtl_modelx,
                             out_path = out_path)
    write.table(smr_retultx2, 
                file = paste0(out_path, "/smr_trait2_", eqtl_modelx, ".txt"),
                sep = "\t", row.names = F, quote = F)
    
    ## plot
    gene_manh <- readRDS(paste0(GENE_PATH, "/gene_loci_manh.rds"))
    p_val1 <- smr_retultx1$P_SMR[match(gene_manh[["LOCI_info"]]$TERM,
                                       smr_retultx1$Gene)]
    p_val2 <- smr_retultx2$P_SMR[match(gene_manh[["LOCI_info"]]$TERM,
                                       smr_retultx2$Gene)]
    
    Manh_plt <- Manh.smr.plt(LOCI_info = gene_manh[["LOCI_info"]],
                             X_axis = gene_manh[["X_axis"]],
                             n_top = 20,
                             p_val1 = p_val1,
                             p_val2 = p_val2,
                             trait1 = "trait1",
                             trait2 = "trait2",
                             thresh = NULL)
    ggsave(paste0(out_path, "/Manh_SMR_", eqtl_modelx, "_plt.png"),
           Manh_plt,
           height = 8, width = 10, dpi = 300)
    system(paste0("rm -rf ", out_path, "/*.ma"))
    system(paste0("rm -rf ", out_path, "/*.smr"))
    system(paste0("rm -rf ", out_path, "/*.list"))
    return(c(paste0(out_path, "/smr_trait1_", eqtl_modelx, ".txt"), 
             paste0(out_path, "/smr_trait2_", eqtl_modelx, ".txt"),
             paste0(out_path, "/Manh_SMR_", eqtl_modelx, "_plt.png")))
  })
  
}  

# #

smr.call(opt$info)
