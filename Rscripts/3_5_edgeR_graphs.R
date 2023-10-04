
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
# suppressMessages(library("ggforce", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))


opt = NULL

# options(warn = 1)

Graphs_counts = function(option_list)
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
  
  #### Master_path_analysis----
  
  Master_path_analysis = opt$Master_path_analysis
  
  cat("OUT_\n")
  cat(sprintf(as.character(Master_path_analysis)))
  cat("\n")
  
  #### Read manifest_FINAL_with_results file ----
  
  setwd(Master_path_analysis)
  
  design_df<-readRDS(file='design_df.RDS')
  
  cat("design_df_0\n")
  cat(str(design_df))
  cat("\n")
  cat(sprintf(as.character(names(summary(design_df$sample)))))
  cat("\n")
  cat(sprintf(as.character(summary(design_df$sample))))
  cat("\n")
  cat(sprintf(as.character(names(summary(design_df$treatment)))))
  cat("\n")
  cat(sprintf(as.character(summary(design_df$treatment))))
  cat("\n")
  cat(sprintf(as.character(names(summary(design_df$time_point)))))
  cat("\n")
  cat(sprintf(as.character(summary(design_df$time_point))))
  cat("\n")
  cat(sprintf(as.character(names(summary(design_df$Genotype)))))
  cat("\n")
  cat(sprintf(as.character(summary(design_df$Genotype))))
  cat("\n")
  
  
  #### Master_path_analysis----
  
  setwd(Master_path_analysis)
 
  
  logCPM.obs<-as.data.frame(fread(file='CPM_obs.tsv', sep="\t"), stringsAsFactors=F)
  colnames(logCPM.obs)[which(colnames(logCPM.obs) == 'V1')]<-'ensembl_gene_id'

  # cat("logCPM.obs_0\n")
  # cat(str(logCPM.obs))
  # cat("\n")

  logCPM.fit<-as.data.frame(fread(file='CPM_fit.tsv', sep="\t"), stringsAsFactors=F)
  colnames(logCPM.fit)<-colnames(logCPM.obs) #[which(colnames(logCPM.fit) == 'V1')]<-'ensembl_gene_id'
  
  # cat("logCPM.fit_0\n")
  # cat(str(logCPM.fit))
  # cat("\n")

  POST_NOR_subset<-as.data.frame(fread(file='POST_NOR_subset.tsv', sep="\t"), stringsAsFactors=F)
  colnames(POST_NOR_subset)[which(colnames(POST_NOR_subset) == 'V1')]<-'ensembl_gene_id'
  
  multiVals <- function(x) paste(x,collapse=";")
  POST_NOR_subset$Symbol <- mapIds(org.Hs.eg.db, keys=POST_NOR_subset$ensembl_gene_id, keytype="ENSEMBL",
                           column="SYMBOL", multiVals=multiVals)
  
  # cat("POST_NOR_subset_0\n")
  # cat(str(POST_NOR_subset))
  # cat("\n")
  
  
  #### LOOP TO PRINT ----
  
  path_graphs<-paste(Master_path_analysis,'graphs','/',sep='')
  
  if (file.exists(path_graphs)){
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
  path_graphs<-paste(Master_path_analysis,'graphs','/','Per_gene_graphs','/',sep='')
  
  if (file.exists(path_graphs)){
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
  ENSG_array<-unique(POST_NOR_subset$ensembl_gene_id)
  
  cat("ENSG_array_0\n")
  cat(str(ENSG_array))
  cat("\n")
  
  
  test_index<-which(ENSG_array == 'ENSG00000005961')
  
  START<-test_index
  START<-1
  
  DEBUG<-0
  
  for(i in START:length(ENSG_array))
  {
    ENSG_array_sel<-ENSG_array[i]
    
    # if(DEBUG == 1)
    # {
      cat("----------------------------->\t")
      cat(sprintf(as.character(i)))
      cat("\t")
      cat(sprintf(as.character(ENSG_array_sel)))
      cat("\t")
      
    # }
    
    POST_NOR_subset_sel<-POST_NOR_subset[which(POST_NOR_subset$ensembl_gene_id == ENSG_array_sel),]
    
    # if(DEBUG == 1)
    # {
    #   cat("POST_NOR_subset_sel_0\n")
    #   cat(str(POST_NOR_subset_sel))
    #   cat("\n")
    # }
    
    Symbol_sel<-unique(POST_NOR_subset_sel$Symbol)
    
    # if(DEBUG == 1)
    # {
      cat(sprintf(as.character(Symbol_sel)))
      cat("\n")
      
    # }
    
    POST_NOR_subset_sel.m<-melt(POST_NOR_subset_sel, id.vars=c('ensembl_gene_id','Symbol'), variable.name='Sample_label', value.name='counts')
    POST_NOR_subset_sel.m$counts<-as.numeric(POST_NOR_subset_sel.m$counts)
    POST_NOR_subset_sel.m$Sample_label<-as.character(POST_NOR_subset_sel.m$Sample_label)
    POST_NOR_subset_sel.m$log2_counts<-log2(POST_NOR_subset_sel.m$counts+1)
    
    if(DEBUG == 1)
    {
      cat("POST_NOR_subset_sel.m_0\n")
      cat(str(POST_NOR_subset_sel.m))
      cat("\n")
    }
    
    logCPM.obs_sel<-logCPM.obs[which(logCPM.obs$ensembl_gene_id == ENSG_array_sel),]
    
    # if(DEBUG == 1)
    # {
    #   cat("logCPM.obs_sel_0\n")
    #   cat(str(logCPM.obs_sel))
    #   cat("\n")
    # }
    
    logCPM.obs_sel.m<-melt(logCPM.obs_sel, id.vars='ensembl_gene_id', variable.name='Sample_label', value.name='CPM')
    logCPM.obs_sel.m$CPM<-as.numeric(logCPM.obs_sel.m$CPM)
    logCPM.obs_sel.m$Sample_label<-as.character(logCPM.obs_sel.m$Sample_label)

    if(DEBUG == 1)
    {
      cat("logCPM.obs_sel.m_0\n")
      cat(str(logCPM.obs_sel.m))
      cat("\n")
    }
    
    logCPM.fit_sel<-logCPM.fit[which(logCPM.fit$ensembl_gene_id == ENSG_array_sel),]
    
    # if(DEBUG == 1)
    # {
    #   cat("logCPM.fit_sel_0\n")
    #   cat(str(logCPM.fit_sel))
    #   cat("\n")
    # }
    
    logCPM.fit_sel.m<-melt(logCPM.fit_sel, id.vars='ensembl_gene_id', variable.name='Sample_label', value.name='CPM_fit')
    logCPM.fit_sel.m$CPM_fit<-as.numeric(logCPM.fit_sel.m$CPM_fit)
    logCPM.fit_sel.m$Sample_label<-as.character(logCPM.fit_sel.m$Sample_label)
    
    if(DEBUG == 1)
    {
      cat("logCPM.fit_sel_0\n")
      cat(str(logCPM.fit_sel.m))
      cat("\n")
    }
    
    
    #### REP merge ----
    
    REP<-merge(POST_NOR_subset_sel.m,
               logCPM.obs_sel.m,
               by=c("Sample_label","ensembl_gene_id"))
    
    if(DEBUG == 1)
    {
      cat("REP_0\n")
      cat(str(REP))
      cat("\n")
    }
    
    REP<-merge(REP,
               logCPM.fit_sel.m,
               by=c("Sample_label","ensembl_gene_id"))
    
    if(DEBUG == 1)
    {
      cat("REP_1\n")
      cat(str(REP))
      cat("\n")
    }
    
    REP<-merge(REP,
               design_df,
               by="Sample_label")
    
    if(DEBUG == 1)
    {
      cat("REP_2\n")
      cat(str(REP))
      cat("\n")
    }
    
    indx.int<-c(which(colnames(REP) == 'ensembl_gene_id'),which(colnames(REP) == 'Symbol'),
                which(colnames(REP) == 'log2_counts'),which(colnames(REP) == 'CPM'),which(colnames(REP) == 'CPM_fit'),
                which(colnames(REP) == 'Sample_label'),which(colnames(REP) == 'sample'),which(colnames(REP) == 'treatment'),
                which(colnames(REP) == 'time_point'),which(colnames(REP) == 'Genotype'))
    
    REP_subset<-unique(REP[order(REP$sample,REP$treatment,REP$time_point),indx.int])

   if(DEBUG ==1)
   {
     cat("REP_subset_0\n")
     cat(str(REP_subset))
     cat("\n")

   }
    
    #### Graphs ----
    
    ### log2_counts
    
    A<-round(summary(REP_subset$log2_counts[!is.na(REP_subset$log2_counts)]),2)
    
    step<-abs(A[6]-A[1])/10
    
    if(step == 0)
    {
      
      step<-1
    }
    
    if(DEBUG ==1)
    {
      cat("Summary_log2_counts\n")
      cat(sprintf(as.character(names(A))))
      cat("\n")
      cat(sprintf(as.character(A)))
      cat("\n")
      cat(sprintf(as.character(step)))
      cat("\n")
    }
    
   
    
    breaks.Rank<-unique(sort(c(0,5,seq(from= A[1], to=A[6]+step,by=step))))
    labels.Rank<-as.character(round(breaks.Rank,1))
    
    
    if(DEBUG == 1)
    {
      cat("breaks.Rank:\t")
      cat(sprintf(as.character(breaks.Rank)))
      cat("\n")
      
      cat("labels.Rank:\t")
      cat(sprintf(as.character(labels.Rank)))
      cat("\n")
      
    }
    
    
    
    graph_log2_counts<-ggplot(data=REP_subset,aes(x=time_point, y=log2_counts, color=sample)) +
      geom_jitter(width = 0.4,size=4)+
      theme_bw()+
      scale_y_continuous(name="log2(counts+1)",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
      scale_color_brewer(palette = "Set3", drop=F)
    
    graph_log2_counts<-graph_log2_counts+
      facet_grid(. ~ REP_subset$Genotype, scales='free_x', space='free_x') +
      theme_cowplot(font_size = 14)+
      theme( strip.background = element_blank(),
             strip.placement = "inside",
             strip.text = element_text(size=14),
             panel.spacing = unit(0.2, "lines"),
             panel.background=element_rect(fill="white"),
             panel.border=element_rect(colour="black",size=1),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())+
      theme(axis.title.y=element_text(size=18, family="sans"),
            axis.title.x=element_text(size=18, family="sans"),
            axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
            axis.text.x=element_blank(),
            legend.title=element_text(size=16,color="black", family="sans"),
            legend.text=element_text(size=12,color="black", family="sans"))+
      theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
            legend.key.height = unit(1.5, 'cm'), #change legend key height
            legend.key.width = unit(1, 'cm'), #change legend key width
            legend.title = element_text(size=14), #change legend title font size
            legend.text = element_text(size=14))+ #change legend text font size
      theme(legend.position = "hidden")+
      scale_x_discrete(name=NULL, drop=T)+
      ggeasy::easy_center_title()
    
    if(DEBUG == 1)
    {
      cat("Graph_log2_counts:\t")
      
    }
    
    ### CPM
    
    A<-round(summary(REP_subset$CPM[!is.na(REP_subset$CPM)]),2)
    
    step<-abs(A[6]-A[1])/10
    
    if(step == 0)
    {
      
      step<-1
    }
    
    if(DEBUG ==1)
    {
      cat("Summary_CPM\n")
      cat(sprintf(as.character(names(A))))
      cat("\n")
      cat(sprintf(as.character(A)))
      cat("\n")
      cat(sprintf(as.character(step)))
      cat("\n")
    }
    
    
    
    breaks.Rank<-unique(sort(c(0,1.5,seq(from= A[1], to=A[6]+step,by=step))))
    labels.Rank<-as.character(round(breaks.Rank,1))
    
    
    if(DEBUG == 1)
    {
      cat("breaks.Rank:\t")
      cat(sprintf(as.character(breaks.Rank)))
      cat("\n")
      
      cat("labels.Rank:\t")
      cat(sprintf(as.character(labels.Rank)))
      cat("\n")
      
    }
    
    
    
    graph_CPM<-ggplot(data=REP_subset,aes(x=time_point, y=CPM, color=sample)) +
      geom_jitter(width = 0.4,size=4)+
      theme_bw()+
      scale_y_continuous(name="CPM",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
      scale_color_brewer(palette = "Set3", drop=F)
    
    graph_CPM<-graph_CPM+
      facet_grid(. ~ REP_subset$Genotype, scales='free_x', space='free_x') +
      theme_cowplot(font_size = 14)+
      theme( strip.background = element_blank(),
             strip.placement = "inside",
             strip.text = element_text(size=14),
             panel.spacing = unit(0.2, "lines"),
             panel.background=element_rect(fill="white"),
             panel.border=element_rect(colour="black",size=1),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())+
      theme(axis.title.y=element_text(size=18, family="sans"),
            axis.title.x=element_text(size=18, family="sans"),
            axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
            axis.text.x=element_blank(),
            legend.title=element_text(size=16,color="black", family="sans"),
            legend.text=element_text(size=12,color="black", family="sans"))+
      theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
            legend.key.height = unit(1.5, 'cm'), #change legend key height
            legend.key.width = unit(1, 'cm'), #change legend key width
            legend.title = element_text(size=14), #change legend title font size
            legend.text = element_text(size=14))+ #change legend text font size
      theme(legend.position = "hidden")+
      scale_x_discrete(name=NULL, drop=T)+
      ggeasy::easy_center_title()
    
    if(DEBUG == 1)
    {
      cat("Graph_CPM:\t")
      
    }
    
    ### CPM_fit
    
    A<-round(summary(REP_subset$CPM_fit[!is.na(REP_subset$CPM_fit)]),2)
    
    step<-abs(A[6]-A[1])/10
    
    if(step == 0)
    {
      
      step<-1
    }
    
    if(DEBUG ==1)
    {
      cat("Summary_CPM_fit\n")
      cat(sprintf(as.character(names(A))))
      cat("\n")
      cat(sprintf(as.character(A)))
      cat("\n")
      cat(sprintf(as.character(step)))
      cat("\n")
    }
    
    
    
    breaks.Rank<-unique(sort(c(0,1.5,seq(from= A[1], to=A[6]+step,by=step))))
    labels.Rank<-as.character(round(breaks.Rank,1))
    
    
    if(DEBUG == 1)
    {
      cat("breaks.Rank:\t")
      cat(sprintf(as.character(breaks.Rank)))
      cat("\n")
      
      cat("labels.Rank:\t")
      cat(sprintf(as.character(labels.Rank)))
      cat("\n")
      
    }
    
    
    
    graph_CPM_fit<-ggplot(data=REP_subset,aes(x=time_point, y=CPM_fit, color=sample)) +
      geom_jitter(width = 0.4,size=4)+
      theme_bw()+
      scale_y_continuous(name="CPM_fit",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
      scale_color_brewer(palette = "Set3", drop=F)
    
    graph_CPM_fit<-graph_CPM_fit+
      facet_grid(. ~ REP_subset$Genotype, scales='free_x', space='free_x') +
      theme_cowplot(font_size = 14)+
      theme( strip.background = element_blank(),
             strip.placement = "inside",
             strip.text = element_text(size=14),
             panel.spacing = unit(0.2, "lines"),
             panel.background=element_rect(fill="white"),
             panel.border=element_rect(colour="black",size=1),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())+
      theme(axis.title.y=element_text(size=18, family="sans"),
            axis.title.x=element_text(size=18, family="sans"),
            axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
            axis.text.x=element_text(angle=45,size=10,hjust=1,vjust=1, color="black", family="sans"),
            legend.title=element_text(size=16,color="black", family="sans"),
            legend.text=element_text(size=12,color="black", family="sans"))+
      theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
            legend.key.height = unit(1.5, 'cm'), #change legend key height
            legend.key.width = unit(1, 'cm'), #change legend key width
            legend.title = element_text(size=14), #change legend title font size
            legend.text = element_text(size=14))+ #change legend text font size
      theme(legend.position="bottom",legend.title=element_blank(), legend.text = element_text(size=12))+
      guides(colour = guide_legend(nrow = 1))+
      scale_x_discrete(name=NULL, drop=T)+
      ggeasy::easy_center_title()
    
    if(DEBUG == 1)
    {
      cat("Graph_CPM_fit:\t")
      
    }
    
    graph_DEF<-plot_grid(graph_log2_counts,graph_CPM,graph_CPM_fit,
                         ncol = 1,
                         nrow=3,
                         rel_heights=c(0.66,0.66,0.95))
    
    setwd(path_graphs)
    
    svgname<-paste(paste(Symbol_sel,ENSG_array_sel,"Model_plot", sep='_'),".svg",sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      ggsave(svgname, plot= graph_DEF,
             device="svg",
             height=10, width=12)
    }
    
    if(DEBUG == 1)
    {
      cat("Graph_Part_END:\t")
      
    }
  }#i in 1:length(ENSG_array)
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
    make_option(c("--manifest_FINAL_with_results"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--gene_expression_results"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--selection_samples"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--selection_treatment"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--selection_time_point"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Master_path_analysis"), type="character", default=NULL, 
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
  
 
  Graphs_counts(opt)


  
}


###########################################################################

system.time( main() )