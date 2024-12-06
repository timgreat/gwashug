# Estimate local genetic correlation

# Load packages
library(plyr)
library(dplyr)
library(bigreadr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(LAVA)
library(optparse)

# Fix parameter
REF_PATH <- "/gwashug/reference/1kg/"
BLK_PATH <- "/gwashug/reference/lava_block"


# Input parameters
args_list = list(
  make_option("--info", type="character", default=NULL,
              help="INPUT: summary information",
              metavar="character")
)
opt_parser <- OptionParser(option_list = args_list)
opt <- parse_args(opt_parser)

# opt <- list(info = "/gwashug/scripts/cogwas/test/summary_info.txt")

# Function 1: process summary_info.txt
summaryInfo.process <- function(summaryInfo_file
){
  
  summary_info <- read.table(summaryInfo_file, sep = "\t",quote="")
  
  summary_list <- list(summ1 = summary_info[1, 2], 
                       pop1 = summary_info[2, 2],
                       trait1 = summary_info[3, 2],
                       type1 = summary_info[4, 2],
                       summ2 = summary_info[5, 2], 
                       pop2 = summary_info[6, 2],
                       trait2 = summary_info[7, 2], 
                       type2 = summary_info[8, 2])

  if (any(grepl("Case1", summary_info[, 1]))){
    
    summary_list$case1 <- summary_info[grep("Case1", summary_info[, 1]), 2] %>% as.numeric
  }
  if (any(grepl("Case2", summary_info[, 1]))){
    
    summary_list$case2 <- summary_info[grep("Case2", summary_info[, 1]), 2] %>% as.numeric
  }

  return(summary_list)
}


# Function 2: Manh.plt
Manh.plt <- function(LOCI_info = NULL,
                     X_axis = NULL,
                     p_val = NULL,
                     Thresh = NULL
){
   
  LOCI_info$P <- p_val[, 1]
  LOCI_info$P2 <- p_val[, 2]
  if (is.null(Thresh)) {
    Thresh = -log10(0.05/sum(LOCI_info$P2 != -1, na.rm=T))
  }
  LOCI_info <- subset(LOCI_info, rowSums(is.na(LOCI_info)) == 0)
  LOCI_info$Sig <- 0
  LOCI_info$Sig[-log10(LOCI_info$P) >= Thresh] <- 1
  LOCI_info$Sig[LOCI_info$P2 == -1] <- 0
  LOCI_info$CHR <- factor(LOCI_info$CHR, levels = c(1:22))
  
  # Manhattan plot
  plt <- ggplot(LOCI_info) +
    geom_point(aes(x = BPcum, y = -log10(P), color = CHR),
               alpha = 0.8, size = 1.5)+
    scale_x_continuous(label = X_axis$CHR, 
                       breaks= X_axis$center, 
                       limits = range(LOCI_info$BPcum)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_color_manual(values = rep(c("#DE423D", "#3D5587"), 11)) +
    xlab("Chromosome") + ylab(bquote("-"~Log[10]~"(P value)"))+
    geom_hline(yintercept = Thresh, 
               color = 'red', size = 0.8, linetype = 2) + 
    geom_label_repel(data = LOCI_info[LOCI_info$Sig == 1,],
                     aes(label = TERM, x = BPcum, y = -log10(P)),
                     segment.size = 0.25, size = 3)+
    theme_bw() +
    theme(legend.position="none",
          axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 12, face = "bold", color = "black"),
          panel.border = element_blank(),
          axis.line.y = element_line(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())

  return(plt)
}

# Function 3: lava test
LAVA.test <- function(lava_input = NULL,
                      pop = "EUR",
                      out_path = NULL
){
  
  ## ref
  loci_df <- readRDS(paste0(BLK_PATH, "/", pop, "/LAVA_blocks_size1000.rds"))
  
  test_out <- paste0(out_path, "bivar_test.txt")
  plt_out <- paste0(out_path, "Manh_plt.txt")
  ## lava test for each locus
  # seq_along(loci_df$LOC)
  start_t <- proc.time()
  print(start_t)
  bivar_test_list <- lapply(seq_along(loci_df$LOC), function(x){
    
    if (x%%10==0){
      cat(x, "\n")
    }
    locusx <- process.locus(loci_df[x,], lava_input)
    
    gc_loc <- tryCatch({
      run.bivar(locusx)
    }, error = function(e) {
      print(NULL)
    })
     
    if (is.null(gc_loc) == F) {
      uni_loc <- run.univ(locusx)
      gc_loc$h2.pheno1 <- uni_loc[1, 2, 1]
      gc_loc$p.pheno1 <- uni_loc[1, 3, 1]
      gc_loc$h2.pheno2 <- uni_loc[2, 2, 1]
      gc_loc$p.pheno2 <- uni_loc[2, 3, 1]
    }
    return(gc_loc)
  })
  end_t <- proc.time()
  print(end_t)
  
  # names(bivar_test_list) <- loci_df$LOC[936:967]
  names(bivar_test_list) <- loci_df$LOC
  
  bivar_test <- lapply(names(bivar_test_list), function(locx){
    
    bivar_testx <- bivar_test_list[[locx]]
    if (!is.null(bivar_testx)) {
        bivar_testx$loc <- locx
    }
    return(bivar_testx)
    
  }) %>% Reduce("rbind", .)
  bivar_test <- bivar_test[!is.na(bivar_test$rho), ]
  
  ## format output
  sig_n <- bivar_test$p.pheno1 < 0.05/nrow(bivar_test) & bivar_test$p.pheno2 < 0.05/nrow(bivar_test)
  bivar_test$p_adj <- -1
  bivar_test$p_adj[sig_n] <- p.adjust(bivar_test$p[sig_n], 
                                      n = sum(sig_n), 
                                      method = "bonferroni")
  
  return(bivar_test)
}

# Function 4: process summary statistics
LAVA.call <- function(summaryInfo_file
){
  
  #
  summary_file <- summaryInfo.process(summaryInfo_file)
  tmp_path <- paste0(dirname(summaryInfo_file), "/tmpp")
  if(!file.exists(tmp_path)){
    
    system(paste0("mkdir ", tmp_path))
  }
  #
  summ_in1 <- fread2(summary_file$summ1)
  summ_out1 <- data.frame(SNP = summ_in1$rs, 
                          A1 = summ_in1$allele1, 
                          A2 = summ_in1$allele0,
                          N = summ_in1$n_obs, 
                          Z = summ_in1$beta/summ_in1$se)
  write.table(summ_out1, file = paste0(tmp_path, "/summ1.txt"), 
              row.names = F, quote = F)
  summ_in2 <- fread2(summary_file$summ2)
  summ_out2 <- data.frame(SNP = summ_in2$rs, 
                          A1 = summ_in2$allele1, 
                          A2 = summ_in2$allele0,
                          N = summ_in2$n_obs, 
                          Z = summ_in2$beta/summ_in2$se)
  write.table(summ_out2, file = paste0(tmp_path, "/summ2.txt"), 
              row.names = F, quote = F)
  
  case_vec <- control_vec <- c(NA, NA)
  if(!is.null(summary_file$case1)){
    
    case_vec[1] <- summary_file$case1
    control_vec[1] <- max(summ_out1$N) - summary_file$case1
  }
  if(!is.null(summary_file$case2)){
    
    case_vec[2] <- summary_file$case2
    control_vec[2] <- max(summ_out2$N) - summary_file$case2
  }

  input_info <- data.frame(phenotype = c(summary_file$trait1, summary_file$trait2), 
                           cases = case_vec, 
                           controls = control_vec, 
                           filename = c(paste0(tmp_path, "/summ1.txt"), 
                                        paste0(tmp_path, "/summ2.txt")))
  write.table(input_info, file = paste0(tmp_path, "/input.info.txt"), 
              row.names = F, quote = F, sep = " ")
  
  ### Read in summary statistics and related info
  input = process.input(input.info.file = paste0(tmp_path, "/input.info.txt"),
                        sample.overlap.file = NULL,
                        ref.prefix = paste0(REF_PATH, "/", summary_file$pop1, "/hm3_imp/merge"),
                        phenos = c(summary_file$trait1, summary_file$trait2))
  ## LAVA test
  out_path <- paste0(dirname(summaryInfo_file), "/output")
  if(!file.exists(out_path)){

    system(paste0("mkdir ", out_path))
  }
  bivar_test <- LAVA.test(lava_input = input,
                          pop = summary_file$pop1,
                          out_path = out_path)
  write.csv(bivar_test, file = paste0(out_path, "/bivar_test.csv"),
            row.names = F, quote = F)

  ## plot
  pop <- summary_file$pop1
  loci_manh <- readRDS(paste0(BLK_PATH, "/", pop,  "/block_loci_manh_size1000.rds"))
  p_val <- bivar_test[match(loci_manh[["LOCI_info"]]$TERM,
                            bivar_test$loc), c(9, 15)]
  Manh_plt <- Manh.plt(LOCI_info = loci_manh[["LOCI_info"]],
                       X_axis = loci_manh[["X_axis"]],
                       p_val = p_val)
  ggsave(paste0(out_path, "/Manh_plt.png"),
         Manh_plt,
         height = 5, width = 10, dpi = 300)
  system(paste0("rm -rf ", tmp_path,"/"))
}

# Run
LAVA.call(opt$info)

