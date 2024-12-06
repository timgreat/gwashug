# Perform coloc analysis
# Load packages
library(bigreadr)
library(dplyr)
library(coloc)
library(stringr)
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
REF_PATH <- "/gwashug/reference/lava_block/"
COLOC_P1 <- 1E-04 
COLOC_P2 <- 1E-04 
COLOC_P3 <- 1E-05
WIN_SIZE <- 0

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
                       type2 = summary_info[8, 2], 
                       PPH4 = as.numeric(summary_info[9, 2]))
  return(summary_list)
}


# Function 2: perform coloc test
coloc.test <- function(summ1 = NULL,
                       summ2 = NULL,
                       pop = "EUR",
                       type1 = "cc",
                       type2 = "cc", 
                       PPH4 = 0.95
){
  ## block file
  loci_df <- readRDS(paste0(REF_PATH, pop, "/LAVA_blocks.rds"))
  
  ## chr rs ps n_mis n_obs allele1 allele0 af beta se p_wald
  coloc_list <- lapply(1:nrow(loci_df), function(x){
    
    ## set summ subset for analysis
    blockx <- loci_df[x,]
    summ1x <- subset(summ1, summ1$chr == blockx$CHR& 
                       summ1$p_wald < 1 &
                       summ1$ps >= blockx$START - WIN_SIZE &
                       summ1$ps <= blockx$STOP + WIN_SIZE)
    summ1x <- summ1x[rowSums(is.na(summ1x)) == 0,]
    summ1x <- summ1x[!duplicated(summ1x$rs),]
    
    #
    summ2x <- subset(summ2, summ2$chr == blockx$CHR& 
                       summ2$p_wald < 1 &
                       summ2$ps >= blockx$START - WIN_SIZE &
                       summ2$ps <= blockx$STOP + WIN_SIZE)
    summ2x <- summ2x[rowSums(is.na(summ2x)) == 0,]
    summ2x <- summ2x[!duplicated(summ2x$rs),]
    
    #
    if (nrow(summ1x) > 0 & nrow(summ2x) > 0) {
      
      ## summ1 input format
      inter_snp <- intersect(summ1x$rs, summ2x$rs)
      if (length(inter_snp) > 0) {
        
        N1 <- summ1x$n_obs[1] + summ1x$n_mis[1]
        MAF1 <- ifelse(summ1x$af < 0.5,
                       summ1x$af, 1 - summ1x$af)
        dataset1 <- list(snp = summ1x$rs,
                         beta = summ1x$beta,
                         varbeta = (summ1x$se)^2,
                         N = N1,
                         MAF = MAF1,
                         type = type1)
        
        ## summ2 input format
        N2 <- summ2x$n_obs[1] + summ2x$n_mis[1]
        MAF2 <- ifelse(summ2x$af < 0.5,
                       summ2x$af, 1 - summ2x$af)
        dataset2 <- list(snp = summ2x$rs,
                         beta = summ2x$beta,
                         varbeta = (summ2x$se)^2,
                         N = N2,
                         MAF = MAF2,
                         type = type2)
        
        ## colocation analysis
        coloc_result <- coloc.abf(dataset1, 
                                  dataset2,
                                  p1 = COLOC_P1, 
                                  p2 = COLOC_P2, 
                                  p12 = COLOC_P3)
        
        # return(coloc_result)
        
        ## 95% credible set to contain a shared casual snp
        o <- order(coloc_result$results$SNP.PP.H4, decreasing = TRUE)
        cs <- cumsum(coloc_result$results$SNP.PP.H4[o])
        w <- which(cs > PPH4)[1]
        CI_snplist <- paste(coloc_result$results[o, ][1:w, ]$snp,
                            collapse = ", ")
        SNP_result <- coloc_result$results[, c("snp", "SNP.PP.H4")]
        SNP_result$LOC <- blockx$LOC
        
        return(list("LOC" = blockx$LOC, 
                    "summary" = coloc_result$summary, 
                    "CI_snplist" = CI_snplist,
                    "SNP_result" = SNP_result))
      } else {
        
        return(NULL)
        
      }
      
    } else {
      
      return(NULL)
      
    }
  })
  
  coloc_out <- lapply(coloc_list, function(coloc_outx){
    if (!is.null(coloc_outx)) {
      coloc_outx[["summary"]]
    }
  }) %>% Reduce("rbind", .) %>% as.data.frame()
  coloc_out$LOC = lapply(coloc_list, function(coloc_outx){
    if (!is.null(coloc_outx)) {
      coloc_outx[["LOC"]]
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
  
  return(list("coloc_out" = coloc_out,
              "coloc_snp_out" = coloc_snp_out))
}

manh.coloc.plt <- function(coloc_out = NULL,
                           pop = "EUR",
                           n_top = 10){
  ## load ref information for manh plot
  loci_manh <- readRDS(paste0(REF_PATH, pop, "/block_loci_manh.rds"))
  LOCI_info <- loci_manh$LOCI_info[match(coloc_out$LOC, loci_manh$LOCI_info$TERM), ]
  LOCI_info$CHR <- as.factor(LOCI_info$CHR)
  ## add PP.H4 and inverse of snp number
  LOCI_info$PP.H4 <- coloc_out$PP.H4.abf
  LOCI_info$inver_snp <- str_split(coloc_out$CI_snplist, ",") %>%
    lapply(., function(listx){
      1/length(listx)
    }) %>% unlist
  LOCI_info_top <- top_n(LOCI_info, n = n_top, wt = PP.H4)
  
  ## plot
  plt <- ggplot(LOCI_info, aes(x = BPcum, y = PP.H4))+
    geom_point(aes(size = inver_snp, color = CHR))+
    scale_radius()+
    scale_x_continuous(label = loci_manh$X_axis$CHR, 
                       breaks= loci_manh$X_axis$center, 
                       limits = range(LOCI_info$BPcum)) +
    scale_color_manual(values = rep(c("royalblue4", "darksalmon"), 11)) +
    xlab("Chromosome") + ylab("PP.H4")+
    geom_hline(yintercept = 0.5, 
               color = 'gold', size = 0.8, linetype = 2) + 
    scale_y_continuous(limits = c(0, 1))+
    geom_text_repel(data = LOCI_info_top,
                    aes(label = TERM), max.iter = 1E8,
                    segment.size = 0.25, size = 3.5, color = "black")+
    guides(size = guide_legend(title = "Inverse of SNP number",
                               override.aes = list(alpha = 1, shape = 21)),
           color = "none")+
    theme_bw()+
    theme(legend.position = "right",
          legend.title = element_text(size = 11, face = "bold", color = "black"),
          legend.text = element_text(size = 10, color = "black"),
          axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 12, face = "bold", color = "black"),
          panel.border = element_blank(),
          axis.line.y = element_line(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  
  
  return(plt)
  
}


# Function 3: call coloc analysis
coloc.call <- function(summaryInfo_file){
  
  # 
  summary_file <- summaryInfo.process(summaryInfo_file)
  out_path <- paste0(dirname(summaryInfo_file), "/output/")
  #
  summ_in1 <- fread2(summary_file$summ1)
  summ_in2 <- fread2(summary_file$summ2)
  
  #
  coloc_result <- coloc.test(summ1 = summ_in1,
                             summ2 = summ_in2,
                             pop = summary_file$pop1,
                             type1 = summary_file$type1,
                             type2 = summary_file$type2, 
                             PPH4 = summary_file$PPH4)
  # plot
  Manh_coloc_plt <- manh.coloc.plt(coloc_out = coloc_result[["coloc_out"]], 
                                   pop = summary_file$pop1,
                                   n_top = 10)
  #
  if(!file.exists(out_path)){
    
    system(paste0("mkdir ", out_path))
  }
  write.table(coloc_result[["coloc_out"]], 
              file = paste0(out_path, "/coloc_out.txt"),
              sep = "\t", row.names = F, quote = F)
  write.table(coloc_result[["coloc_snp_out"]], 
              file = paste0(out_path, "/coloc_snp_out.txt"),
              sep = "\t", row.names = F, quote = F)
  ggsave(paste0(out_path, "/Manh_coloc_plt.png"),
         Manh_coloc_plt,
         height = 6, width = 13, dpi = 300)
}

# # Run
coloc.call(opt$info)
