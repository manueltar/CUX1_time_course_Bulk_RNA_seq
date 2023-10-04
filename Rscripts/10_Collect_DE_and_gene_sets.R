
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


opt = NULL

options(warn = 1)

collects_DE = function(option_list)
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
  
  #### READ and transform Master_path_analysis ----
  
  Master_path_analysis = opt$Master_path_analysis
  
  cat("OUT_\n")
  cat(sprintf(as.character(Master_path_analysis)))
  cat("\n")
  
  #### READ and transform analysis_array ----
  
  analysis_array = unlist(strsplit(opt$analysis_array, split=","))
  
  cat("analysis_array_\n")
  cat(sprintf(as.character(analysis_array)))
  cat("\n")
  
  #### LOOP comparisons ----
  
  DEBUG <- 0
  
  list_results<-list()
  
  for(i in 1:length(analysis_array))
  {
    analysis_array_sel<-analysis_array[i]
    
    cat("--------comparison-------------->\t")
    cat(sprintf(as.character(analysis_array_sel)))
    cat("\n")
    
    
    path_comparisons<-paste(Master_path_analysis,analysis_array_sel,'/',sep='')
    
    if (file.exists(path_comparisons)){
      
      
    }else{
      
      
      
    }#path_comparisons
    
    #### Open DE file ----
    
    setwd(path_comparisons)
    
    DE_file<-as.data.frame(fread(file=paste("results_",analysis_array_sel,'.tsv', sep=''), sep="\t", header=T), stringsAsFactors=F)
    
    if(DEBUG == 1)
    {
      cat("DE_file_0\n")
      cat(str(DE_file))
      cat("\n")
    }
    
    LogFC_indx<-grep("logFC\\.",colnames(DE_file))
    
    if(DEBUG == 1)
    {
      cat("LogFC_indx_0\n")
      cat(sprintf(as.character((LogFC_indx))))
      cat("\n")
    }
    
    indx.int<-c(LogFC_indx,
                which(colnames(DE_file) == 'FDR'),
                which(colnames(DE_file) == 'Symbol'),
                which(colnames(DE_file) == 'ensembl_gene_id'))
    if(DEBUG == 1)
    {
      cat("indx.int_0\n")
      cat(sprintf(as.character((indx.int))))
      cat("\n")
    }
    
    DE_file_subset<-unique(DE_file[,indx.int])
    
    if(DEBUG == 1)
    {
      cat("DE_file_subset_0\n")
      cat(str(DE_file_subset))
      cat("\n")
    }
    
    
    DE_file_subset$Minus_logpval<--1*log10(DE_file_subset$FDR)
    
    if(DEBUG == 1)
    {
      cat("DE_file_subset_0\n")
      cat(str(DE_file_subset))
      cat("\n")
    }
    
    indx.dup<-which(duplicated(DE_file_subset$ensembl_gene_id) == "TRUE")
    
    if(DEBUG == 1)
    {
      cat("indx.dup_0\n")
      cat(str(indx.dup))
      cat("\n")
    }
    
    check<-DE_file_subset[indx.dup,]
    
    if(DEBUG == 1)
    {
      cat("check_0\n")
      cat(str(check))
      cat("\n")
    }
    
    
    DE_file_subset.m<-melt(DE_file_subset[,-which(colnames(DE_file_subset) == 'FDR')], id.vars=c("ensembl_gene_id","Symbol","Minus_logpval"),
                    variable.name='variable', value.name='logFC')
    
    DE_file_subset.m$time_point<-as.character(DE_file_subset.m$variable)
    
    if(DEBUG == 1)
    {
      cat("DE_file_subset.m_0\n")
      cat(str(DE_file_subset.m))
      cat("\n")
    }
    
    DE_file_subset.m$time_point<-gsub("VS","\\.",DE_file_subset.m$time_point)
    DE_file_subset.m$time_point<-gsub("logFC\\.","",DE_file_subset.m$time_point)
    DE_file_subset.m$time_point<-gsub("homALT\\.","",DE_file_subset.m$time_point)
    DE_file_subset.m$time_point<-gsub("wt\\.","",DE_file_subset.m$time_point)
    DE_file_subset.m$time_point<-gsub("Del80\\.","",DE_file_subset.m$time_point)
    DE_file_subset.m$time_point<-gsub("Time_course\\.","",DE_file_subset.m$time_point)
    DE_file_subset.m$time_point<-gsub("Basal\\.Basal","Basal",DE_file_subset.m$time_point)
    DE_file_subset.m$time_point<-gsub("16_hrs\\.16_hrs","16_hrs",DE_file_subset.m$time_point)
    DE_file_subset.m$time_point<-gsub("24_hrs\\.24_hrs","24_hrs",DE_file_subset.m$time_point)
    DE_file_subset.m$time_point<-gsub("48_hrs\\.48_hrs","48_hrs",DE_file_subset.m$time_point)
    DE_file_subset.m$time_point<-gsub("72_hrs\\.72_hrs","72_hrs",DE_file_subset.m$time_point)
    DE_file_subset.m$time_point<-gsub("\\.","__",DE_file_subset.m$time_point)
    
    if(DEBUG == 1)
    {
      cat("DE_file_subset.m_1\n")
      cat(str(DE_file_subset.m))
      cat("\n")
    }
    
    
    DE_file_subset.m$time_point<-factor(DE_file_subset.m$time_point,
                                                  levels=c('Basal','16_hrs','24_hrs','48_hrs','72_hrs',
                                                           'Basal__16_hrs','Basal__24_hrs','Basal__48_hrs','Basal__72_hrs',
                                                           '16_hrs__24_hrs','16_hrs__48_hrs','16_hrs__72_hrs',
                                                           '24_hrs__48_hrs','24_hrs__72_hrs',
                                                           '48_hrs__72_hrs'),
                                                  ordered=T)
    
    DE_file_subset.m<-droplevels(unique(DE_file_subset.m[,-which(colnames(DE_file_subset.m) == 'variable')]))
    DE_file_subset.m$comparison<-analysis_array_sel
    
    if(DEBUG == 1)
    {
      cat("DE_file_subset.m_2\n")
      cat(str(DE_file_subset.m))
      cat("\n")
      cat(sprintf(as.character(names(summary(DE_file_subset.m$time_point)))))
      cat("\n")
      cat(sprintf(as.character(summary(DE_file_subset.m$time_point))))
      cat("\n")
    }
    
    list_results[[i]]<-DE_file_subset.m
    
  }# i in 1:length(analysis_array)
  
  
  if(length(list_results) >0)
  {
    FINAL_df = unique(as.data.frame(data.table::rbindlist(list_results, fill = T)))
    
    FINAL_df$comparison<-factor(FINAL_df$comparison,
                                levels=c("homALT_vs_wt","Del80_vs_wt","homALT_vs_Del80","Time_course_wt"),
                                ordered=T)
    
  
    cat("FINAL_df_0\n")
    cat(str(FINAL_df))
    cat("\n")
    cat(sprintf(as.character(names(summary(FINAL_df$time_point)))))
    cat("\n")
    cat(sprintf(as.character(summary(FINAL_df$time_point))))
    cat("\n")
    cat(sprintf(as.character(names(summary(FINAL_df$comparison)))))
    cat("\n")
    cat(sprintf(as.character(summary(FINAL_df$comparison))))
    cat("\n")
    
    
   
    #### Open DE file ----
    
    setwd(Master_path_analysis)
    
    write.table(FINAL_df,file="DE_edgeR_results.tsv", sep="\t", quote=F, row.names = F)
    saveRDS(FINAL_df,file="DE_edgeR_results.rds")
    
    
    
  }#length(list_results) >0
 
}

collects_GSEA = function(option_list)
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
  
  #### READ and transform minGS_size ----
  
  minGS_size = opt$minGS_size
  
  cat("minGS_size_\n")
  cat(sprintf(as.character(minGS_size)))
  cat("\n")
  
  #### READ and transform maxGS_size ----
  
  maxGS_size = opt$maxGS_size
  
  cat("maxGS_size_\n")
  cat(sprintf(as.character(maxGS_size)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform Master_path_analysis ----
  
  Master_path_analysis = opt$Master_path_analysis
  
  cat("OUT_\n")
  cat(sprintf(as.character(Master_path_analysis)))
  cat("\n")
  
  #### READ and transform analysis_array ----
  
  analysis_array = unlist(strsplit(opt$analysis_array, split=","))
  
  cat("analysis_array_\n")
  cat(sprintf(as.character(analysis_array)))
  cat("\n")
  
  #### LOOP comparisons ----
  
  DEBUG <- 0
  
  list_results<-list()
  
  for(i in 1:length(analysis_array))
  {
    analysis_array_sel<-analysis_array[i]
    
    cat("--------comparison-------------->\t")
    cat(sprintf(as.character(analysis_array_sel)))
    cat("\n")
    
    
    path_comparisons<-paste(Master_path_analysis,analysis_array_sel,'/','MySigDb_analysis','/',sep='')
    
    if (file.exists(path_comparisons)){
      
      
    }else{
      
      
      
    }#path_comparisons
    
    #### Open MySigDB_results file ----
    
    setwd(path_comparisons)
    
    MySigDB_results<-as.data.frame(fread(file=paste("result",analysis_array_sel,'_MySigDB_GSEA.tsv', sep=''), sep="\t", header=T), stringsAsFactors=F)
    
    if(DEBUG == 1)
    {
      cat("MySigDB_results_0\n")
      cat(str(MySigDB_results))
      cat("\n")
    }
    
    if(dim(MySigDB_results)[1] >0)
    {
      MySigDB_results_subset<-unique(MySigDB_results[,-which(colnames(MySigDB_results) == 'collection')])
      
      if(DEBUG == 1)
      {
        cat("MySigDB_results_subset_0\n")
        cat(str(MySigDB_results_subset))
        cat("\n")
      }
      
      MySigDB_results_subset.dt<-data.table(MySigDB_results_subset, key=c('term.id','term.name'))
      
      
      MySigDB_results_subset_MIN<-as.data.frame(MySigDB_results_subset.dt[,.SD[which.min(adjusted.p.val)], by=key(MySigDB_results_subset.dt)], stringsAsFactos=F)
      
      if(DEBUG == 1)
      {
        cat("MySigDB_results_subset_MIN_0\n")
        cat(str(MySigDB_results_subset_MIN))
        cat("\n")
      }
      
      MySigDB_results_subset_MIN_subset<-MySigDB_results_subset_MIN[which(MySigDB_results_subset_MIN$term.size >= minGS_size &
                                                                            MySigDB_results_subset_MIN$term.size <= maxGS_size),]
      
      if(DEBUG == 1)
      {
        cat("MySigDB_results_subset_MIN_subset_0\n")
        cat(str(MySigDB_results_subset_MIN_subset))
        cat("\n")
      }
      
      MySigDB_results_subset_MIN_subset$comparison<-analysis_array_sel
      
      MySigDB_results_subset_MIN_subset$Minus_logpval<-round(-1*log10(MySigDB_results_subset_MIN_subset$adjusted.p.val),2)
      
      if(DEBUG == 1)
      {
        cat("MySigDB_results_subset_MIN_subset_1\n")
        cat(str(MySigDB_results_subset_MIN_subset))
        cat("\n")
      }
      
      
      FINAL<-MySigDB_results_subset_MIN_subset[,c(which(colnames(MySigDB_results_subset_MIN_subset) == 'term.id'),which(colnames(MySigDB_results_subset_MIN_subset) == 'overlap'),
                                                  which(colnames(MySigDB_results_subset_MIN_subset) == 'Minus_logpval'),
                                                  which(colnames(MySigDB_results_subset_MIN_subset) == 'comparison'))]
      
      if(DEBUG == 1)
      {
        cat("FINAL_1\n")
        cat(str(FINAL))
        cat("\n")
      }
      
     
      list_results[[i]]<-FINAL
      
    }#dim(MySigDB_results)[1] >0
  }# i in 1:length(analysis_array)
  
  
  if(length(list_results) >0)
  {
    FINAL_df = unique(as.data.frame(data.table::rbindlist(list_results, fill = T)))
    
    FINAL_df$comparison<-factor(FINAL_df$comparison,
                                levels=c("homALT_vs_wt","Del80_vs_wt","homALT_vs_Del80","Time_course_wt"),
                                ordered=T)
    
    
    cat("FINAL_df_0\n")
    cat(str(FINAL_df))
    cat("\n")
    cat(sprintf(as.character(names(summary(FINAL_df$time_point)))))
    cat("\n")
    cat(sprintf(as.character(summary(FINAL_df$time_point))))
    cat("\n")
    cat(sprintf(as.character(names(summary(FINAL_df$comparison)))))
    cat("\n")
    cat(sprintf(as.character(summary(FINAL_df$comparison))))
    cat("\n")
    
    
    
    #### Open DE file ----
    
    setwd(Master_path_analysis)
    
    write.table(FINAL_df,file="GSEA_ActivePathways_results.tsv", sep="\t", quote=F, row.names = F)
    saveRDS(FINAL_df,file="GSEA_ActivePathways_results.rds")
    
    
    
  }#length(list_results) >0
  
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
    make_option(c("--analysis_array"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Master_path_analysis"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--minGS_size"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--maxGS_size"), type="numeric", default=NULL, 
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
  
  collects_DE(opt)
  collects_GSEA(opt)

  
}


###########################################################################

system.time( main() )