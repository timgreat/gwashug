# Perform coloc analysis
# Load packages
library(bigreadr)
library(dplyr)
library(coloc)
library(optparse)
library(ggplot2)
library(ggrepel)
library(patchwork)
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
REF_PATH <- "/gwashug/reference/eqtl/GTEx/"
LD_PATH <- "/gwashug/reference/1kg/EUR/ld/"
COLOC_P1 <- 1E-04 
COLOC_P2 <- 1E-04 
COLOC_P3 <- 1E-05
WIN_PLOT <- 2E05
R2_LABELS <- c("R2 < 0.2", 
               "0.2 < R2 < 0.4",
               "0.4 < R2 < 0.6",
               "0.6 < R2 < 0.8",
               "0.8 < R2")
R2_LEVELS <- seq(0.2, 0.8, 0.2)
LOCUSZOOM_COLOR_TOP <- "purple"
LOCUSZOOM_COLOR <- scale_color_manual(values = c("R2 < 0.2" = "darkblue", 
                                                 "0.2 < R2 < 0.4" = "turquoise", 
                                                 "0.4 < R2 < 0.6" = "green", 
                                                 "0.6 < R2 < 0.8" = "orange", 
                                                 "0.8 < R2" = "red"))
# Function 1: process summary_info.txt
summaryInfo.process <- function(summaryInfo_file
){
  
  summary_info <- read.table(summaryInfo_file, sep = "\t",quote="")
  summary_list <- list(summ1 = summary_info[1, 2], 
                       type1 = summary_info[2, 2],
                       summ2 = summary_info[3, 2], 
                       type2 = summary_info[4, 2],
                       use_model1 = summary_info[5, 2],
                       use_model2 = summary_info[6, 2],
                       PPH4 = as.numeric(summary_info[7, 2]))
  return(summary_list)
}

# Function 2: perform coloc test with select eqtl model
coloc.eqtl.test <- function(summ1 = NULL,
                            type1 = "cc",
                            eqtl_path = NULL,
                            eqtl_model = NULL){
  
  gene_list <- read.table(paste0(eqtl_path, "/eGene_list.txt.gz"), header = F)[, 1, drop = T]
  N_eqtl <- 1000
  N1 <- summ1$n_obs[1] + summ1$n_mis[1]
  coloc_list <- lapply(gene_list, function(genex){
    
    ##
    summ_eqtlx <- fread2(paste0(eqtl_path, "/", genex, "_hm3_cis.txt.gz"))
    summ1x <- subset(summ1, summ1$rs %in% summ_eqtlx$SNP &
                       summ1$p_wald < 1 &
                       summ1$se != 0)
    summ1x <- summ1x[rowSums(is.na(summ1x)) == 0,]
    summ1x <- summ1x[!duplicated(summ1x$rs),]
    
    ## input format
    if (nrow(summ1x) > 0) {
      ## summ1x input format
      MAF1 <- ifelse(summ1x$af < 0.5,
                     summ1x$af, 1 - summ1x$af)
      dataset1 <- list(snp = summ1x$rs,
                       beta = summ1x$beta,
                       varbeta = (summ1x$se)^2,
                       N = N1,
                       MAF = MAF1,
                       type = type1)
      
      ## summ2 input format
      dataset2 <- list(snp = summ_eqtlx$SNP,
                       beta = summ_eqtlx$b,
                       varbeta = (summ_eqtlx$SE)^2,
                       N = N_eqtl,
                       MAF = summ_eqtlx$MAF,
                       type = "quant")
      
      ## colocation analysis
      coloc_result <- coloc.abf(dataset1, 
                                dataset2,
                                p1 = COLOC_P1, 
                                p2 = COLOC_P2, 
                                p12 = COLOC_P3)
      
      ## 95% credible set to contain a shared casual snp
      o <- order(coloc_result$results$SNP.PP.H4,decreasing=TRUE)
      cs <- cumsum(coloc_result$results$SNP.PP.H4[o])
      w <- which(cs > 0.95)[1]
      CI_snplist <- paste(coloc_result$results[o,][1:w,]$snp,
                          collapse = ", ")
      SNP_result <- coloc_result$results[, c("snp", "SNP.PP.H4")]
      SNP_result$Gene <- genex
      
      return(list("Gene" = genex, 
                  "summary" = coloc_result$summary, 
                  "CI_snplist" = CI_snplist,
                  "SNP_result" = SNP_result))
    } else {
      return(NULL)
    }
    
  })
  
  coloc_out <- lapply(coloc_list, function(coloc_outx){
    if (!is.null(coloc_outx)) {
      coloc_outx[["summary"]]
    }
  }) %>% Reduce("rbind", .) %>% as.data.frame()
  coloc_out$Model <- eqtl_model
  coloc_out$Gene = lapply(coloc_list, function(coloc_outx){
    if (!is.null(coloc_outx)) {
      coloc_outx[["Gene"]]
    }
  }) %>% unlist()
  coloc_out$CI_snplist = lapply(coloc_list, function(coloc_outx){
    if (!is.null(coloc_outx)) {
      coloc_outx[["CI_snplist"]]
    }
  }) %>% unlist()
  
  coloc_snp_out <- lapply(coloc_list, function(coloc_outx){
    if (!is.null(coloc_outx)) {
      coloc_outx[["SNP_result"]]
    }
  }) %>% Reduce("rbind", .) %>% as.data.frame()
  coloc_snp_out$Model <- eqtl_model
  
  return(list("coloc_out" = coloc_out,
              "coloc_snp_out" = coloc_snp_out))
}

# Function 3: gene label plot
gene.label.plt <- function(summ = NULL,
                           top_snp = NULL,
                           win_size = WIN_PLOT){
  ## set bp window around shared top SNP
  start_coord <- summ$ps[summ$rs == top_snp] - win_size
  end_coord <- summ$ps[summ$rs == top_snp] + win_size
  chrr <- summ$chr[summ$rs == top_snp]
  ## set information around target gene
  gene_bp <- fread2(paste0(REF_PATH, "/gene_bpf.txt.gz"), header = T, sep = "\t")
  gene_bpf <- gene_bp[gene_bp$CHR == chrr & 
                        gene_bp$START >= start_coord &
                        gene_bp$STOP <= end_coord,]
  gene_bpf$gene_pos <- (gene_bpf$START + gene_bpf$STOP)/2
  
  ## plot 
  plt <- ggplot(data = gene_bpf) + 
    geom_errorbar(aes(xmin = START, xmax = STOP, y = GENE),
                  size = 0.2, width = 0.3, color = "grey50") +
    geom_text(aes(x = gene_pos, y = GENE, label = GENE),
              vjust = 1, size = 2) + 
    scale_x_continuous(limits = c(start_coord, end_coord)) + 
    xlab(paste0("Chromosome ", gene_bpf$CHR)) +
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 10, color = "black"),
          axis.title.x = element_text(size = 15, face = "bold"))
  
  return(plt)
  
}

# Function 4: locuszoom plt 
locus.plt <- function(summ = NULL,
                      top_snp = NULL,
                      snp_list = NULL,
                      win_size = WIN_PLOT){
  ## set bp window around shared top SNP
  start_coord <- summ$ps[summ$rs == top_snp] - win_size
  end_coord <- summ$ps[summ$rs == top_snp] + win_size
  chrr <- summ$chr[summ$rs == top_snp]
  summ_sub <- subset(summ, summ$chr == chrr & 
                       summ$ps >= start_coord &
                       summ$ps <= end_coord &
                       summ$rs %in% snp_list)
  
  ld_chr <- fread2(paste0(LD_PATH, "/chr", chrr, ".ld.gz"))
  ld_chr_sub <- subset(ld_chr, ld_chr$SNP_A == top_snp & 
                         ld_chr$SNP_B %in% summ_sub$rs)
  #
  summ_sub$R2_raw <- 0
  summ_sub$R2_raw[match(ld_chr_sub$SNP_B, summ_sub$rs)] <- ld_chr_sub$R2
  
  ## set R2 group
  summ_sub$p_wald[summ_sub$p_wald == 0] <- min(summ_sub$p_wald[summ_sub$p_wald != 0])
  summ_sub$R2 <- cut(summ_sub$R2_raw, 
                     breaks = c(-Inf, R2_LEVELS, Inf),
                     labels = R2_LABELS) %>% 
    factor(., levels = R2_LABELS)
  
  ## plot
  plt <- ggplot(summ_sub[summ_sub$rs != top_snp, ]) + 
    # top SNP
    geom_point(data = summ_sub[summ_sub$rs == top_snp, ],
               mapping = aes(x = ps, y = -log(p_wald)),
               color = LOCUSZOOM_COLOR_TOP, shape = 18, size = 4) + 
    # other SNPs
    geom_point(aes(x = ps, y = -log(p_wald), color = R2),
               shape = 19, size = 2) + 
    LOCUSZOOM_COLOR + 
    annotate("text",
             x = summ_sub$ps[summ_sub$rs == top_snp] - 
               (end_coord - start_coord) * 0.05,
             y = -log(summ_sub$p_wald[summ_sub$rs == top_snp]), 
             size = 3, colour = "grey10", label = top_snp) +
    ylab(bquote("-"~Log[10]~"(P value)")) + 
    theme_bw() + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10))
  
  return(plt)
}

# Function 3: call coloc analysis
coloc.eqtl.call <- function(summaryInfo_file
){
  # 

  summary_file <- summaryInfo.process(summaryInfo_file)
  out_path <- paste0(dirname(summaryInfo_file), "/output/")
  if (!file.exists(out_path)) {
    system(paste0("mkdir ", out_path))
  }
  #
  summ_in1 <- fread2(summary_file$summ1)
  summ_in2 <- fread2(summary_file$summ2)
  type1 <- summary_file$type1
  type2 <- summary_file$type2
  use_model1 <- summary_file$use_model1
  use_model2 <- summary_file$use_model2
  PPH4 <- summary_file$PPH4
  
  ## eqtl_path_list: eqtl_model  eqtl_path
  use_model <- summary_file$use_model
  eqtl_path_list <- read.table(paste0(REF_PATH, "eqtl_model_list.txt"), header = T)
  # if (is.null(use_model)) {
  #   
  #   use_line <- 1:nrow(eqtl_path_list)
  # } else {
    
  use_line <- match(c(use_model1, use_model2), eqtl_path_list$eqtl_model)
  # }

  #@timi
  if (use_model1==use_model2){
      use_line <- match(c(use_model1), eqtl_path_list$eqtl_model)
  }
  #@timi

  coloc_out <- lapply(use_line, function(x){
    ##
    eqtl_pathx <- eqtl_path_list$eqtl_path[x]
    eqtl_modelx <- eqtl_path_list$eqtl_model[x]
    coloc_result1x <- coloc.eqtl.test(summ1 = summ_in1,
                                      type1 = type1,
                                      eqtl_path = eqtl_pathx,
                                      eqtl_model = eqtl_modelx)
    fwrite2(coloc_result1x[["coloc_out"]],
            file = paste0(out_path, "/coloc_trait1_", eqtl_modelx, ".txt"),
            sep = "\t", row.names = F, quote = F)
    fwrite2(coloc_result1x[["coloc_snp_out"]],
            file = paste0(out_path, "/coloc_snp_trait1_", eqtl_modelx, ".txt"),
            sep = "\t", row.names = F, quote = F)
    ##
    coloc_result2x <- coloc.eqtl.test(summ1 = summ_in2,
                                      type1 = type2,
                                      eqtl_path = eqtl_pathx,
                                      eqtl_model = eqtl_modelx)
    fwrite2(coloc_result2x[["coloc_out"]],
            file = paste0(out_path, "/coloc_trait2_", eqtl_modelx, ".txt"),
            sep = "\t", row.names = F, quote = F)
    fwrite2(coloc_result2x[["coloc_snp_out"]],
            file = paste0(out_path, "/coloc_snp_trait2_", eqtl_modelx, ".txt"),
            sep = "\t", row.names = F, quote = F)
    #coloc_result1x <- as.data.frame(fread2(file=paste0(out_path,"/coloc_test1_",eqtl_modelx,".txt"),sep="\t",quote=F))
    #coloc_result2x <- as.data.frame(fread2(file=paste0(out_path,"/coloc_test2_",eqtl_modelx,".txt"),sep="\t",quote=F))
    ##
    coloc_sig1 <- subset(coloc_result1x[["coloc_out"]],
                         coloc_result1x[["coloc_out"]]$PP.H4.abf > PPH4)
    coloc_sig2 <- subset(coloc_result2x[["coloc_out"]],
                         coloc_result2x[["coloc_out"]]$PP.H4.abf > PPH4)
    sig_gene_list <- intersect(coloc_sig1$Gene, coloc_sig2$Gene)
    ##
    if (length(sig_gene_list) > 0) {
      
      locuszoom_out <- lapply(sig_gene_list, function(sig_genex){
        # load eqtl summ for sig_genex
        summ_eqtl <- fread2(paste0(eqtl_pathx, "/", sig_genex, "_hm3_cis.txt.gz"),
                            select = c("Chr", "SNP", "BP", "p"))
        colnames(summ_eqtl) <- c("chr","rs",  "ps", "p_wald")
        inter_snp <- Reduce("intersect",
                            list(summ_in1$rs,
                                 summ_in2$rs,
                                 summ_eqtl$rs))
        
        ## plot for summ1 and corresponding eQTL
        CI_snplist1 <- strsplit(coloc_sig1$CI_snplist[coloc_sig1$Gene == sig_genex],
                                ",") %>% unlist()
        top_snp1 <- CI_snplist1[1]
        gene_label_plt1 <- gene.label.plt(summ = summ_in1,
                                          top_snp = top_snp1)
        locus_plt1 <- locus.plt(summ = summ_in1,
                                top_snp = top_snp1,
                                snp_list = inter_snp)
        locus_plt1_eqtl <- locus.plt(summ = summ_eqtl,
                                     top_snp = top_snp1,
                                     snp_list = inter_snp)
        
        ## plot for summ1 and corresponding eQTL
        CI_snplist2 <- strsplit(coloc_sig2$CI_snplist[coloc_sig2$Gene == sig_genex],
                                ",") %>% unlist()
        top_snp2 <- CI_snplist2[1]
        gene_label_plt2 <- gene.label.plt(summ = summ_in2,
                                          top_snp = top_snp2)
        locus_plt2 <- locus.plt(summ = summ_in2,
                                top_snp = top_snp2,
                                snp_list = inter_snp)
        locus_plt2_eqtl <- locus.plt(summ = summ_eqtl,
                                     top_snp = top_snp2,
                                     snp_list = inter_snp)
        ## output
        out_listx <- paste0(paste0(out_path, "/Locuszoom_", eqtl_modelx, "_", sig_genex),
                            c("_summ1.png", "_eqtl1.png", "_summ2.png", "_eqtl2.png"))
        ggsave(out_listx[1], 
               (locus_plt1 / gene_label_plt1) + 
                 plot_layout(heights = c(2, 1.2)),
               height = 8, width = 10)
        ggsave(out_listx[2], 
               (locus_plt1_eqtl / gene_label_plt1) + 
                 plot_layout(heights = c(2, 1.2)),
               height = 8, width = 10)
        ggsave(out_listx[3], 
               (locus_plt2 / gene_label_plt2) + 
                 plot_layout(heights = c(2, 1.2)),
               height = 8, width = 10)
        ggsave(out_listx[4], 
               (locus_plt2_eqtl / gene_label_plt2) + 
                 plot_layout(heights = c(2, 1.2)),
               height = 8, width = 10)
        
        return(out_listx)
      }) %>% unlist
    } else {
      
      locuszoom_out <- NA
    }
    
    return(c(paste0(out_path, "/coloc_trait1_", eqtl_modelx, ".txt"),
             paste0(out_path, "/coloc_snp_trait1_", eqtl_modelx, ".txt"),
             paste0(out_path, "/coloc_trait2_", eqtl_modelx, ".txt"),
             paste0(out_path, "/coloc_snp_trait2_", eqtl_modelx, ".txt"),
             locuszoom_out))
  })
  
}

# Run
coloc.eqtl.call(opt$info)
