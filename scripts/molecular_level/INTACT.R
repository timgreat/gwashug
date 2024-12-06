# Perform coloc analysis
# Load packages
library(bigreadr)
library(dplyr)
library(glue)
library(INTACT)
library(optparse)

# Input parameters
args_list = list(
  make_option("--info", type="character", default=NULL,
              help="INPUT: summary information",
              metavar="character")
)
opt_parser <- OptionParser(option_list = args_list)
opt <- parse_args(opt_parser)


EQTL_PATH <- "/gwashug/reference/eqtl/GTEx/"

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

# Function 3: INTACT test
intact.test <- function(coloc_path = NULL,
                        smr_path = NULL
){
  ##
  if (file.exists(coloc_path) & file.exists(smr_path)) {
    print("step1")
    
    ## format output from coloc and SMR
    coloc_out <- fread2(coloc_path, select = c("Gene", "PP.H4.abf"))
    coloc_out <- subset(coloc_out, rowSums(is.na(coloc_out)) == 0)
    
    smr_out <- fread2(smr_path, select = c("Gene", "Beta_SMR", "SE_SMR"))
    smr_out <- subset(smr_out, rowSums(is.na(smr_out)) == 0 &
                        smr_out$SE_SMR > 0)
    smr_out$Z_SMR <- smr_out$Beta_SMR / smr_out$SE_SMR
    inter_gene <- intersect(coloc_out$Gene, smr_out$Gene)
    
    if (length(inter_gene) > 0) {
      print("step2")
      ## INTACT analysis
      intact_df <- data.frame(Gene = inter_gene,
                              PPH4_coloc = coloc_out$PP.H4.abf[match(inter_gene, coloc_out$Gene)],
                              Z_SMR = smr_out$Z_SMR[match(inter_gene, smr_out$Gene)])
      intact_df$PP_INTACT <- intact(GLCP_vec = intact_df$PPH4_coloc,
                                    z_vec = intact_df$Z_SMR)
      
    } else {
      
      intact_df <- data.frame(Gene = NA,
                              PPH4_coloc = NA,
                              Z_SMR = NA)
      
    }
    
  } else {
    
    intact_df <- data.frame(Gene = NA,
                            PPH4_coloc = NA,
                            Z_SMR = NA)
  }
  
  return(intact_df)
}

# Function 3: call INTACT analysis
intact.call <- function(summaryInfo_file){
  
  # 
  summary_file <- summaryInfo.process(summaryInfo_file)
  in_path <- dirname(summaryInfo_file)
  coloceqtl_path <- gsub("INTACT", "COLOCEQTL", in_path)
  smr_path <- gsub("INTACT", "SMR", in_path)
  ## eqtl_path_list: eqtl_model  eqtl_path
  # use_model <- summary_file$use_model
  eqtl_path_list <- read.table(paste0(EQTL_PATH, "eqtl_model_list.txt"), header = T)
  # if (is.null(use_model)) {
  #   use_line <- 1:nrow(eqtl_path_list)
  # } else {
    use_line <- match(c(summary_file$use_model1, summary_file$use_model2), eqtl_path_list$eqtl_model)
  # }
  #@timi
  if (summary_file$use_model1==summary_file$use_model2){
      use_line <- match(c(summary_file$use_model1), eqtl_path_list$eqtl_model)
  }
  #@timi
  intact_out <- lapply(use_line, function(x){
   
    #
    eqtl_modelx <- eqtl_path_list$eqtl_model[x]
    coloc_pathx1 <- paste0(coloceqtl_path, "/output/coloc_trait1_", eqtl_modelx, ".txt")
    coloc_pathx2 <- paste0(coloceqtl_path, "/output/coloc_trait2_", eqtl_modelx, ".txt")
    smr_pathx1 <- paste0(smr_path, "/output/smr_trait1_", eqtl_modelx, ".txt")
    smr_pathx2 <- paste0(smr_path, "/output/smr_trait2_", eqtl_modelx, ".txt")
    ##
    intact_resultx1 <- intact.test(coloc_path = coloc_pathx1,
                                   smr_path = smr_pathx1)
    out_path <- paste0(in_path, "/output")
    if (!file.exists(out_path)){
      
      system(paste0("mkdir ", out_path))
    }
    write.table(intact_resultx1, 
                file = paste0(out_path, "/intact_trait1_", eqtl_modelx, ".txt"),
                sep = "\t", row.names = F, quote = F)
    ##
    intact_resultx2 <- intact.test(coloc_path = coloc_pathx2,
                                   smr_path = smr_pathx2)
    write.table(intact_resultx2, 
                file = paste0(out_path, "/intact_trait2_", eqtl_modelx, ".txt"),
                sep = "\t", row.names = F, quote = F)
    
    return(c(paste0(out_path, "/intact_trait1_", eqtl_modelx, ".txt"), 
             paste0(out_path, "/intact_trait2_", eqtl_modelx, ".txt")))
  })
  
}  

# #
intact.call(opt$info)
