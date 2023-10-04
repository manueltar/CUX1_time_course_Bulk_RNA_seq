
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

collect_results_build_expression_object_genes = function(option_list)
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
  
  #### READ and transform path_to_results ----
  
  path_to_results = opt$path_to_results
  
  cat("path_to_results_\n")
  cat(sprintf(as.character(path_to_results)))
  cat("\n")
  
  #### Read manifest file ----
  
  manifest<-as.data.frame(fread(file=opt$manifest, sep="\t", header=T), stringsAsFactors=F)
  
  manifest$treatment[which(manifest$sample == 'hESC_MK')]<-'untreated'
  manifest$time_point[which(manifest$sample == 'hESC_MK')]<-'Basal'
  
  
  cat("manifest_0\n")
  cat(str(manifest))
  cat("\n")
  
  #### build a df with the files from Rsem quantification ----
  
  setwd(path_to_results)
  
  files_df<-list.files()
  
  cat("files_df_0\n")
  cat(str(files_df))
  cat("\n")
  
  files_df_sel<-files_df[grep("genes\\.results$",files_df)]
  
  cat("files_df_sel_0\n")
  cat(str(files_df_sel))
  cat("\n")
  
  
 
  
  df_files <- data.frame(matrix(NA,
                                nrow = length(files_df_sel),
                                ncol = 2))
  colnames(df_files)[which(colnames(df_files) == 'X1')]<-'Sample_label'
  colnames(df_files)[which(colnames(df_files) == 'X2')]<-'path'
  
  cat("df_files_0\n")
  cat(str(df_files))
  cat("\n")
  
  
  df_files$path<-paste(path_to_results,files_df_sel,sep='')
  df_files$Sample_label<-gsub("\\.genes\\.results$","",files_df_sel)
  
  
  cat("df_files_1\n")
  cat(str(df_files))
  cat("\n")
  
  
  #### Merge manifest with df_files by Sample_label ----
  
  df_files<-merge(df_files,manifest, by="Sample_label")
  
  cat("df_files_2\n")
  cat(str(df_files))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(df_files$sample))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(df_files$sample)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(df_files$treatment))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(df_files$treatment)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(df_files$time_point))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(df_files$time_point)))))
  cat("\n")
  
  #### Create and order factors----
  
  df_files$sample<-factor(df_files$sample,
                          levels=c('WT_A','WT_B','WT_C','clone_13','clone_27','clone_29','del_233','del_235','del_287','hESC_MK'),
                          ordered=T)
  
  
  df_files$treatment<-factor(df_files$treatment,
                          levels=c('untreated','5nM_PMA','50MicroM_Hemin'),
                          ordered=T)
  
  df_files$time_point<-factor(df_files$time_point,
                             levels=c('Basal','16_hrs','24_hrs','48_hrs','72_hrs','96_hrs'),
                             ordered=T)
  
  df_files<-df_files[order(df_files$sample,df_files$treatment,df_files$time_point,df_files$Sample_name),]
  
  Sample_Label_order<-df_files$Sample_label
  
  df_files$Sample_label<-factor(df_files$Sample_label,
                                levels=Sample_Label_order,
                                ordered=T)
  
  cat("df_files_3\n")
  cat(str(df_files))
  cat("\n")
  cat(sprintf(as.character(names(summary(df_files$sample)))))
  cat("\n")
  cat(sprintf(as.character(summary(df_files$sample))))
  cat("\n")
  cat(sprintf(as.character(names(summary(df_files$treatment)))))
  cat("\n")
  cat(sprintf(as.character(summary(df_files$treatment))))
  cat("\n")
  cat(sprintf(as.character(names(summary(df_files$time_point)))))
  cat("\n")
  cat(sprintf(as.character(summary(df_files$time_point))))
  cat("\n")
  cat(sprintf(as.character(names(summary(df_files$Sample_label)))))
  cat("\n")
  cat(sprintf(as.character(summary(df_files$Sample_label))))
  cat("\n")
  
  ###### Open the results files and round to integer ----
  
  
  Sample_label_array<-df_files$Sample_label
  
  
  list_results<-list()
  
  
  DEBUG<-0
  
  for(i in 1:length(Sample_label_array))
  {
    Sample_label_array_sel<-Sample_label_array[i]
    
    cat("------>\t")
    cat(sprintf(as.character(Sample_label_array_sel)))
    cat("\n")
    
    df_files_sel<-df_files[which(df_files$Sample_label == Sample_label_array_sel),]
    
    if(DEBUG == 1)
    {
      cat("df_files_sel_0\n")
      cat(str(df_files_sel))
      cat("\n")
      
    }
    
    
    
    if(dim(df_files_sel)[1] >0)
    {
      sel_file<-df_files_sel$path
      
      cat("--->\t")
      cat(sprintf(as.character(sel_file)))
      cat("\n")
      
      SIZE_gate<-file.info(sel_file)$size
      
      if(DEBUG == 1)
      {
        cat("SIZE_gate_0\n")
        cat(str(SIZE_gate))
        cat("\n")
      }
      
      if(SIZE_gate> 0)
      {
        LINE_gate<-length(readLines(sel_file))
        
        if(DEBUG == 1)
        {
          cat("LINE_gate_0\n")
          cat(str(LINE_gate))
          cat("\n")
        }
        
        if(LINE_gate> 0)
        {
          results<-as.data.frame(fread(file=sel_file, sep="\t", header = T, fill=TRUE), stringsAsFactors = F)
          
          if(DEBUG == 1)
          {
            cat("results_0\n")
            cat(str(results))
            cat("\n")
          }
          
          indx.int<-c(which(colnames(results) == 'gene_id'),which(colnames(results) == 'expected_count'))
          
          results_subset<-unique(results[,indx.int])
          
          colnames(results_subset)[which(colnames(results_subset) == 'gene_id')]<-'ensembl_gene_id'
          results_subset$expected_count<-round(results_subset$expected_count,0)
          results_subset$Sample_label<-Sample_label_array_sel
            
          if(DEBUG == 1)
          {
            cat("results_subset_0\n")
            cat(str(results_subset))
            cat("\n")
          }
          
          list_results[[i]]<-results_subset
          
        }#LINE_gate> 0
      }# SIZE_gate> 0
    }# dim(df_files_sel)[1] >0
  }# i in 1:length(Sample_label_array)
  
  if(length(list_results) >0)
  {
    df_LONG = unique(as.data.frame(data.table::rbindlist(list_results, fill = T)))
    
    df_LONG$expected_count<-as.integer(df_LONG$expected_count)
   
    cat("df_LONG_0\n")
    cat(str(df_LONG))
    cat("\n")
    
    
    
    df_wide<-as.data.frame(pivot_wider(df_LONG,
                                       id_cols=ensembl_gene_id,
                                       names_from=Sample_label,
                                       values_from=expected_count), stringsAsFactors=F)



    cat("df_wide_0\n")
    cat(str(df_wide))
    cat("\n")
    
    df_wide[is.na(df_wide)]<-0
    
    cat("df_wide_1\n")
    cat(str(df_wide))
    cat("\n")
    
    
    #### SAVE ----
    
    setwd(out)
    
    write.table(df_wide,file="Gene_expression_results.tsv", sep="\t", quote=F, row.names = F)
    saveRDS(df_files, file="manifest_FINAL_with_results.rds")
    
    
    
  }#list_results
}

collect_results_build_expression_object_transcripts = function(option_list)
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
  
  #### READ and transform path_to_results ----
  
  path_to_results = opt$path_to_results
  
  cat("path_to_results_\n")
  cat(sprintf(as.character(path_to_results)))
  cat("\n")
  
  #### Read manifest_FINAL_with_results file  ----
  

  setwd(out)
  
  manifest_FINAL_with_results<-readRDS(file="manifest_FINAL_with_results.rds")
  
  colnames(manifest_FINAL_with_results)[which(colnames(manifest_FINAL_with_results) == 'path')]<-'path_to_gene_results'
  
  
  cat("manifest_FINAL_with_results_0\n")
  cat(str(manifest_FINAL_with_results))
  cat("\n")
  
  #### build a df with the files from Rsem quantification ----
  
  setwd(path_to_results)
  
  files_df<-list.files()
  
  cat("files_df_0\n")
  cat(str(files_df))
  cat("\n")
  
  files_df_sel<-files_df[grep("isoforms\\.results$",files_df)]
  
  cat("files_df_sel_0\n")
  cat(str(files_df_sel))
  cat("\n")
  
  
  
  
  df_files <- data.frame(matrix(NA,
                                nrow = length(files_df_sel),
                                ncol = 2))
  colnames(df_files)[which(colnames(df_files) == 'X1')]<-'Sample_label'
  colnames(df_files)[which(colnames(df_files) == 'X2')]<-'path_to_transcript_results'
  
  cat("df_files_0\n")
  cat(str(df_files))
  cat("\n")
  
  
  df_files$path_to_transcript_results<-paste(path_to_results,files_df_sel,sep='')
  df_files$Sample_label<-gsub("\\.isoforms\\.results$","",files_df_sel)
  
  
  cat("df_files_1\n")
  cat(str(df_files))
  cat("\n")
  
  
  #### Merge manifest_FINAL_with_results with df_files by Sample_label ----
  
  df_files<-merge(df_files,manifest_FINAL_with_results, by="Sample_label")
  
  
  df_files<-df_files[order(df_files$sample,df_files$treatment,df_files$time_point,df_files$Sample_name),]
  
  Sample_Label_order<-df_files$Sample_label
  
  df_files$Sample_label<-factor(df_files$Sample_label,
                                levels=Sample_Label_order,
                                ordered=T)
  
  cat("df_files_0\n")
  cat(str(df_files))
  cat("\n")
  cat(sprintf(as.character(names(summary(df_files$sample)))))
  cat("\n")
  cat(sprintf(as.character(summary(df_files$sample))))
  cat("\n")
  cat(sprintf(as.character(names(summary(df_files$treatment)))))
  cat("\n")
  cat(sprintf(as.character(summary(df_files$treatment))))
  cat("\n")
  cat(sprintf(as.character(names(summary(df_files$time_point)))))
  cat("\n")
  cat(sprintf(as.character(summary(df_files$time_point))))
  cat("\n")
  cat(sprintf(as.character(names(summary(df_files$Sample_label)))))
  cat("\n")
  cat(sprintf(as.character(summary(df_files$Sample_label))))
  cat("\n")
  
  ###### Open the results files and round to integer ----
  
  
  Sample_label_array<-df_files$Sample_label
  
  
  list_results<-list()
  
  
  DEBUG<-0
  
  for(i in 1:length(Sample_label_array))
  {
    Sample_label_array_sel<-Sample_label_array[i]
    
    cat("------>\t")
    cat(sprintf(as.character(Sample_label_array_sel)))
    cat("\n")
    
    df_files_sel<-df_files[which(df_files$Sample_label == Sample_label_array_sel),]
    
    if(DEBUG == 1)
    {
      cat("df_files_sel_0\n")
      cat(str(df_files_sel))
      cat("\n")
      
    }
    
    
    
    if(dim(df_files_sel)[1] >0)
    {
      sel_file<-df_files_sel$path_to_transcript_results
      
      cat("--->\t")
      cat(sprintf(as.character(sel_file)))
      cat("\n")
      
      SIZE_gate<-file.info(sel_file)$size
      
      if(DEBUG == 1)
      {
        cat("SIZE_gate_0\n")
        cat(str(SIZE_gate))
        cat("\n")
      }
      
      if(SIZE_gate> 0)
      {
        LINE_gate<-length(readLines(sel_file))
        
        if(DEBUG == 1)
        {
          cat("LINE_gate_0\n")
          cat(str(LINE_gate))
          cat("\n")
        }
        
        if(LINE_gate> 0)
        {
          results<-as.data.frame(fread(file=sel_file, sep="\t", header = T, fill=TRUE), stringsAsFactors = F)
          
          if(DEBUG == 1)
          {
            cat("results_0\n")
            cat(str(results))
            cat("\n")
          }
          
          indx.int<-c(which(colnames(results) == 'gene_id'),which(colnames(results) == 'transcript_id'),which(colnames(results) == 'expected_count'))
          
          results_subset<-unique(results[,indx.int])
          
          colnames(results_subset)[which(colnames(results_subset) == 'gene_id')]<-'ensembl_gene_id'
          results_subset$expected_count<-round(results_subset$expected_count,0)
          results_subset$Sample_label<-Sample_label_array_sel
          
          if(DEBUG == 1)
          {
            cat("results_subset_0\n")
            cat(str(results_subset))
            cat("\n")
          }
          
          list_results[[i]]<-results_subset
          
          
        }#LINE_gate> 0
      }# SIZE_gate> 0
    }# dim(df_files_sel)[1] >0
  }# i in 1:length(Sample_label_array)
  
  if(length(list_results) >0)
  {
    df_LONG = unique(as.data.frame(data.table::rbindlist(list_results, fill = T)))
    
    df_LONG$expected_count<-as.integer(df_LONG$expected_count)
    
    cat("df_LONG_0\n")
    cat(str(df_LONG))
    cat("\n")
    
    
    
    df_wide<-as.data.frame(pivot_wider(df_LONG,
                                       id_cols=c("ensembl_gene_id","transcript_id"),
                                       names_from=Sample_label,
                                       values_from=expected_count), stringsAsFactors=F)
    
    
    
    cat("df_wide_0\n")
    cat(str(df_wide))
    cat("\n")
    
    df_wide[is.na(df_wide)]<-0
    
    cat("df_wide_1\n")
    cat(str(df_wide))
    cat("\n")
    
    df_wide<-df_wide[order(df_wide$ensembl_gene_id),]
    
    cat("df_wide_2\n")
    cat(str(df_wide))
    cat("\n")
    
    transcript_wide<-unique(df_wide[,-c(which(colnames(df_wide) == 'ensembl_gene_id'))])
    
    cat("transcript_wide_0\n")
    cat(str(transcript_wide))
    cat("\n")
    
    key_gene_ENST<-unique(df_wide[,c(which(colnames(df_wide) == 'ensembl_gene_id'),which(colnames(df_wide) == 'transcript_id'))])
    
    cat("key_gene_ENST_0\n")
    cat(str(key_gene_ENST))
    cat("\n")
    
    #### SAVE ----
    
    setwd(out)
    
    write.table(transcript_wide,file="Transcript_expression_results.tsv", sep="\t", quote=F, row.names = F)
    write.table(key_gene_ENST,file="Gene_transcript_correspondence.tsv", sep="\t", quote=F, row.names = F)
    saveRDS(df_files, file='manifest_FINAL_with_results.rds')
    

    
  }#list_results
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
    make_option(c("--path_to_results"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--manifest"), type="character", default=NULL, 
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
  
  collect_results_build_expression_object_genes(opt)
  collect_results_build_expression_object_transcripts(opt)
  
  
}


###########################################################################

system.time( main() )