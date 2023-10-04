
suppressMessages(library("plyr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("data.table", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("crayon", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggplot2", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("farver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("labeling", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("optparse", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("dplyr", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("backports", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("broom", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rstudioapi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cli", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tzdb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggeasy", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("sandwich", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("digest", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tidyverse", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cowplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ActivePathways", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))


opt = NULL

options(warn = 1)

AP_function = function(option_list)
{

  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ edgeR_results ----
  
  edgeR_results<-as.data.frame(fread(file=opt$edgeR_results, sep="\t", header=T), stringsAsFactors=F)
  
  cat("edgeR_results_0\n")
  cat(str(edgeR_results))
  cat("\n")
  
  edgeR_results_NO_NA<-edgeR_results[!is.na(edgeR_results$Symbol),]
  
  cat("edgeR_results_NO_NA_0\n")
  cat(str(edgeR_results_NO_NA))
  cat("\n")
  
  indx.DUP<-which(duplicated(edgeR_results_NO_NA$Symbol) == TRUE)
  
  cat("indx.DUP_0\n")
  cat(str(indx.DUP))
  cat("\n")
  
  edgeR_results_NO_NA_NO_DUP<-edgeR_results_NO_NA[-indx.DUP,]
  
  cat("edgeR_results_NO_NA_NO_DUP_0\n")
  cat(str(edgeR_results_NO_NA_NO_DUP))
  cat("\n")
  
  scores<-matrix(edgeR_results_NO_NA_NO_DUP$FDR)
  row.names(scores)<-edgeR_results_NO_NA_NO_DUP$Symbol
  colnames(scores)<-"Adj_pval"
  
  cat("scores_0\n")
  cat(str(scores))
  cat("\n")
  
  #### READ and transform out ----
  
  background_genes = opt$background_genes
  
  HPA<-read.GMT(background_genes)
  
  # cat("HPA_0\n")
  # cat(str(HPA))
  # cat("\n")
  # 
  background<-makeBackground(HPA)
  
  cat("background_0\n")
  cat(str(background))
  cat("\n")
  
  #### READ and transform out ----
  
  path_to_GMT = opt$path_to_GMT
  
  cat("path_to_GMT_\n")
  cat(sprintf(as.character(path_to_GMT)))
  cat("\n")
  
  #### READ and transform search_terms ----
  
  search_terms = unlist(strsplit(opt$search_terms, split=","))
  
  cat("search_terms_\n")
  cat(sprintf(as.character(search_terms)))
  cat("\n")
  
  
  #### list all files in path_to_GMT ----
  
  file_list <- list.files(path=path_to_GMT, include.dirs = FALSE)
  
  
  cat("file_list\n")
  cat(str(file_list))
  cat("\n")
  
  
  indexes_sel <- grep("Hs\\.symbols\\.gmt$",file_list)
  
  file_list_sel <- as.data.frame(file_list[indexes_sel], stringsAsFactors=F)
  colnames(file_list_sel)<-"file"
  
  cat("file_list_sel\n")
  cat(str(file_list_sel))
  cat("\n")
  
  ### Loop to open gmt files ---
  
  START<-grep('c3.tft.gtrd.v2023.1.Hs.symbols.gmt',file_list_sel$file)
  
  START<-1
  
  cat("START\n")
  cat(str(START))
  cat("\n")
  
  tested<-list()
  
  DEBUG<-0
  
  setwd(path_to_GMT)
  
  Final_result<-list()
  
  for(i in START:dim(file_list_sel)[1])
  {
    file_sel<-file_list_sel$file[i]
    
    cat("------>\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(file_sel)))
    cat("\n")
    
    
    gmt_sel<-read.GMT(file_sel)
    
    # cat(sprintf(as.character(length(gmt_sel))))
    # cat("\n")
    
    # if(DEBUG == 1)
    # {
    #   cat("gmt_sel_0\n")
    #   cat(str(gmt_sel))
    #   cat("\n")
    # }
    
    id_vector<-vector()
    
    for(k in 1:length(gmt_sel))
    {
      id_list<-unique(gmt_sel[[k]]$id)
      
      # if(DEBUG == 1)
      # {
      #   cat("id_list_0\n")
      #   cat(str(id_list))
      #   cat("\n")
      #   
      # }
      
      id_vector[k]<-paste(id_list, collapse='__')
      
      # if(DEBUG == 1)
      # {
      #   cat("id_vector_0\n")
      #   cat(str(id_vector))
      #   cat("\n")
      #   
      # }
      
    }#k in 1:length(gmt_sel)
    
    if(DEBUG == 1)
    {
      cat("id_vector_0\n")
      cat(str(id_vector))
      cat("\n")
    }
    
    # toMatch<-c("Chr1p13","cHr1p21")
    
    toMatch<-search_terms
    
    matches <- grep(paste(toMatch,collapse="|"),id_vector)
    
    if(DEBUG == 1)
    {
      cat("matches_0\n")
      cat(str(matches))
      cat("\n")
    }
    
    toMatch<-tolower(search_terms)
    # toMatch<-tolower(c("Chr1p13","cHr1p21"))
    
    matches_lc <- grep(paste(toMatch,collapse="|"),id_vector)
    
    if(DEBUG == 1)
    {
      cat("matches_lc_0\n")
      cat(str(matches_lc))
      cat("\n")
    }
    
    
    total_matches<-unique(c(matches,matches_lc))
    
    gmt_sel_GREP<-gmt_sel[total_matches]
    
    tested[[i]]<-gmt_sel_GREP
    
    if(DEBUG == 1)
    {
      cat("gmt_sel_GREP_0\n")
      cat(str(gmt_sel_GREP))
      cat("\n")
    }
    
    if(length(gmt_sel_GREP) >0)
    {
      AP_result<-ActivePathways(scores, gmt_sel_GREP , background= background)
      
      if(DEBUG == 1)
      {
        cat("AP_result_0\n")
        cat(str(AP_result))
        cat("\n")
      }
      
      FLAG_null<-sum(is.null(AP_result))
      
      if(FLAG_null == 0)
      {
        AP_result_df<-as.data.frame(AP_result, stringsAsFactors=F)
        
        AP_result_df$overlap<-as.character(paste(unlist(AP_result_df$overlap), collapse=";"))
        
        AP_result_df$collection<-file_sel
        
        if(DEBUG == 1)
        {
          cat("AP_result_df_0\n")
          cat(str(AP_result_df))
          cat("\n")
        }
        
        Final_result[[i]]<-AP_result_df
        
      }#FLAG_null == 0
    }#length(gmt_sel_GREP) >1
  }# i in 1:dim(file_list_sel)[1]
  
  if(length(Final_result)>0)
  {
    FINAL_df = unique(as.data.frame(data.table::rbindlist(Final_result, fill = T)))
    
    cat("FINAL_df_0\n")
    cat(str(FINAL_df))
    cat("\n")
    
    ALL_tested = unique(as.data.frame(data.table::rbindlist(tested, fill = T)))
    
    cat("ALL_tested_0\n")
    cat(str(ALL_tested))
    cat("\n")
    
    path_MySigDb<-paste(out,'MySigDb_analysis','/',sep='')
    
    if (file.exists(path_MySigDb)){
      
    }else{
      
      dir.create(file.path(path_MySigDb))
      
    }#path_MySigDb
    
    setwd(path_MySigDb)
    
    write.table(FINAL_df,file=paste("result",type,".tsv", sep=''),sep="\t", quote=F, row.names = F)
    saveRDS(ALL_tested,file='ALL_tested_gmt.rds')
    
  }#length(Final_result)>0
}




printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}


#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--path_to_GMT"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--search_terms"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--edgeR_results"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--background_genes"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  AP_function(opt)
  
  
  
}


###########################################################################

system.time( main() )