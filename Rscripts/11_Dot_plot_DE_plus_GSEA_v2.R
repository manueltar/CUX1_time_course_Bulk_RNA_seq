
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
# suppressMessages(library("ggrepel", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

opt = NULL

options(warn = 1)

Dot_plot_function_CTRL_gene = function(option_list)
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
  
 
  #### READ and transform ctrl_gene ----
  
  ctrl_gene = unlist(strsplit(opt$ctrl_gene, split=","))
  
  cat("ctrl_gene_\n")
  cat(sprintf(as.character(ctrl_gene)))
  cat("\n")
  
  #### READ and transform plot_together ----
  
  plot_together = unlist(strsplit(opt$plot_together, split=","))
  
  cat("plot_together_\n")
  cat(sprintf(as.character(plot_together)))
  cat("\n")
  
  #### READ and transform plot_ctrl ----
  
  plot_ctrl = unlist(strsplit(opt$plot_ctrl, split=","))
  
  cat("plot_ctrl_\n")
  cat(sprintf(as.character(plot_ctrl)))
  cat("\n")
  
  #### Read the edgeR results file ----
  
  MySigDB_results_FINAL<-readRDS(file=opt$MySigDB_results_FINAL)
  
  MySigDB_results_FINAL$representation<-NA
  
  MySigDB_results_FINAL$representation[which(MySigDB_results_FINAL$comparison%in%plot_together)]<-"plot_together"
  MySigDB_results_FINAL$representation[which(MySigDB_results_FINAL$comparison%in%plot_ctrl)]<-"plot_ctrl"
  
  MySigDB_results_FINAL$representation<-factor(MySigDB_results_FINAL$representation,
                                               levels=c("plot_together","plot_ctrl"),
                                               ordered=T)
  
  
  cat("MySigDB_results_FINAL_0\n")
  cat(str(MySigDB_results_FINAL))
  cat("\n")
  
  #### Read the edgeR results file ----
  
  edgeR_results_FINAL<-readRDS(file=opt$edgeR_results_FINAL)
  
  cat("edgeR_results_FINAL_0\n")
  cat(str(edgeR_results_FINAL))
  cat("\n")
  
  edgeR_results_FINAL$representation<-NA
  
  edgeR_results_FINAL$representation[which(edgeR_results_FINAL$comparison%in%plot_together)]<-"plot_together"
  edgeR_results_FINAL$representation[which(edgeR_results_FINAL$comparison%in%plot_ctrl)]<-"plot_ctrl"
  
  edgeR_results_FINAL$representation<-factor(edgeR_results_FINAL$representation,
                                               levels=c("plot_together","plot_ctrl"),
                                               ordered=T)
  
  cat("edgeR_results_FINAL_1\n")
  cat(str(edgeR_results_FINAL))
  cat("\n")
 
  
 
  #### PLOT CTRL GENE ----
 
  REP_sel<-droplevels(edgeR_results_FINAL[which(edgeR_results_FINAL$ensembl_gene_id%in%ctrl_gene),])
  
  cat("REP_sel_0\n")
  cat(str(REP_sel))
  cat("\n")
  
  array_representations<-levels(edgeR_results_FINAL$representation)
  
  DEBUG<-0
  
  for(i in 1:length(array_representations))
  {
    array_representations_sel<-array_representations[i]
    
    cat("------------------------>\t")
    cat(sprintf(as.character((i))))
    cat("\t")
    cat(sprintf(as.character((array_representations_sel))))
    cat("\n")
    
    

    REP_sel_subset<-droplevels(unique(REP_sel[which(REP_sel$representation == array_representations_sel),]))

    REP_sel_subset$Significance<-NA
    
    REP_sel_subset$Significance[which(REP_sel_subset$Minus_logpval >= 1.3)]<-'YES'
    REP_sel_subset$Significance[which(REP_sel_subset$Minus_logpval < 1.3)]<-'NO'
    
    REP_sel_subset$Significance<-factor(REP_sel_subset$Significance,
                                        levels=c('NO','YES'))
    
    if(DEBUG == 1)
    {
      cat("REP_sel_subset_0\n")
      cat(str(REP_sel_subset))
      cat("\n")
    }
    
    ### graph parameters_FC
    
    indx_FC<-which(colnames(REP_sel_subset) == 'logFC')
    
    A_FC<-summary(REP_sel_subset[,indx_FC])
    
    
    if(DEBUG == 1)
    {
      cat("A_FC\n")
      cat(sprintf(as.character(names(A_FC))))
      cat("\n")
      cat(sprintf(as.character(A_FC)))
      cat("\n")
    }
    
    max_value<-A_FC[6]
    min_value<-A_FC[1]
    
    
    step<-round(abs(max_value-min_value)/10,0)
    
    if(step == 0)
    {
      
      step<-1
    }
    breaks_FC<-unique(sort(unique(c(0,seq(min_value,max_value+step, by=step)))))
    labels_FC<-as.character(round(breaks_FC),1)
    
    if(DEBUG == 1)
    {
      cat("step_FC\n")
      cat(sprintf(as.character(step)))
      cat("\n")
      cat("labels_FC\n")
      cat(sprintf(as.character(labels_FC)))
      cat("\n")
    }
    
    ### graph parameters_Minus_logpval
    
    indx_Minus_logpval<-which(colnames(REP_sel_subset) == 'Minus_logpval')
    
    A_Minus_logpval<-summary(REP_sel_subset[,indx_Minus_logpval])
    
    
    if(DEBUG == 1)
    {
      cat("A_Minus_logpval\n")
      cat(sprintf(as.character(names(A_Minus_logpval))))
      cat("\n")
      cat(sprintf(as.character(A_Minus_logpval)))
      cat("\n")
    }
    
    max_value<-A_Minus_logpval[6]
    min_value<-A_Minus_logpval[1]
    
    
    step<-round(abs(max_value-min_value)/10,0)
    
    if(step == 0)
    {
      
      step<-1
    }
    breaks_Minus_logpval<-unique(sort(unique(c(0,seq(min_value,max_value+step, by=step)))))
    labels_Minus_logpval<-as.character(round(breaks_Minus_logpval),1)
    
    if(DEBUG == 1)
    {
      cat("step_Minus_logpval\n")
      cat(sprintf(as.character(step)))
      cat("\n")
      cat("labels_Minus_logpval\n")
      cat(sprintf(as.character(labels_Minus_logpval)))
      cat("\n")
    }
    
    #### dotplot
    
    if(DEBUG == 1)
    {
      cat("Graph_Part_START:\n")
      
    }
    
    dotplot<-ggplot(data=REP_sel_subset,
                    aes(y=Symbol,
                        x=time_point)) +
      geom_point(aes(color=Significance,
                     fill=logFC,
                     size=Minus_logpval), stroke=1, shape=21)+
      scale_color_manual(values=c("gray","black"),name='Significant', drop=F)+
      scale_size(range = c(0,20), name='-log10pval',
                 breaks=breaks_Minus_logpval, labels=labels_Minus_logpval, limits=c(breaks_Minus_logpval[1],breaks_Minus_logpval[length(breaks_Minus_logpval)]))+
      scale_fill_gradient2(
        low = "blue", 
        mid = "white", 
        high = "red", 
        midpoint = 0,
        breaks=breaks_FC,labels=labels_FC,
        limits=c(breaks_FC[1]-0.01,breaks_FC[length(breaks_FC)]+0.01),name=paste('LogFC',sep="\n"),na.value = "gray")+
      scale_y_discrete(name=NULL, drop=F)+
      theme_classic()+
      ggeasy::easy_center_title()
    
    if(DEBUG == 1)
    {
      cat("Graph_Part_MIDDLE:\n")
      
      cat("REP_sel_subset_MIDDLE\n")
      cat(str(REP_sel_subset))
      cat("\n")
      
    }
    
    dotplot<-dotplot+
      facet_grid(. ~ comparison, scales='free_x', space='free_x',
                 drop=F) +
      theme_cowplot(font_size = 14)+
      theme( strip.background = element_blank(),
             strip.placement = "inside",
             strip.text = element_text(size=14, angle=90),
             panel.spacing = unit(0.2, "lines"),
             panel.background=element_rect(fill="white"),
             panel.border=element_rect(colour="black",size=1),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())+
      theme(axis.title.y=element_blank(),
            axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
            axis.text.x=element_text(angle=45,size=14, vjust=1, hjust=1,color="black", family="sans"))+
      theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
            legend.key.height = unit(1.5, 'cm'), #change legend key height
            legend.key.width = unit(1, 'cm'), #change legend key width
            legend.title = element_text(size=14), #change legend title font size
            legend.text = element_text(size=14))+ #change legend text font size
      scale_x_discrete(name=NULL, drop=T)

    if(DEBUG == 1)
    {
      cat("Graph_Part_END_I:\n")
      
    }
  
    
    path_graphs<-paste(out,'Dotplots','/',sep='')
    
    if (file.exists(path_graphs)){
      
      
    }else{
      
      dir.create(file.path(path_graphs))
      
    }#path_graphs
    
    setwd(path_graphs)
    
    svgname<-paste(paste("Dotplot",'CTRL_genes',paste(array_representations_sel, collapse="__"), sep='_'),".svg",sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      ggsave(svgname, plot= dotplot,
             device="svg",
             height=10, width=12)
    }

    if(DEBUG == 1)
    {
      cat("Graph_Part_END_II:\n")
      
      

    }
   
  }#i in 1:length(LogFC_indx)
}

Dot_plot_function_MSigDB = function(option_list)
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
  
 
  
  #### READ and transform ctrl_gene ----
  
  ctrl_gene = unlist(strsplit(opt$ctrl_gene, split=","))
  
  cat("ctrl_gene_\n")
  cat(sprintf(as.character(ctrl_gene)))
  cat("\n")
  
  #### READ and transform plot_together ----
  
  plot_together = unlist(strsplit(opt$plot_together, split=","))
  
  cat("plot_together_\n")
  cat(sprintf(as.character(plot_together)))
  cat("\n")
  
  #### READ and transform plot_ctrl ----
  
  plot_ctrl = unlist(strsplit(opt$plot_ctrl, split=","))
  
  cat("plot_ctrl_\n")
  cat(sprintf(as.character(plot_ctrl)))
  cat("\n")
  
  #### Read the edgeR results file ----
  
  MySigDB_results_FINAL<-readRDS(file=opt$MySigDB_results_FINAL)
  
  MySigDB_results_FINAL$representation<-NA
  
  MySigDB_results_FINAL$representation[which(MySigDB_results_FINAL$comparison%in%plot_together)]<-"plot_together"
  MySigDB_results_FINAL$representation[which(MySigDB_results_FINAL$comparison%in%plot_ctrl)]<-"plot_ctrl"
  
  MySigDB_results_FINAL$representation<-factor(MySigDB_results_FINAL$representation,
                                               levels=c("plot_together","plot_ctrl"),
                                               ordered=T)
  
  
  cat("MySigDB_results_FINAL_0\n")
  cat(str(MySigDB_results_FINAL))
  cat("\n")
  
  #### Read the edgeR results file ----
  
  edgeR_results_FINAL<-readRDS(file=opt$edgeR_results_FINAL)
  
  cat("edgeR_results_FINAL_0\n")
  cat(str(edgeR_results_FINAL))
  cat("\n")
  
  edgeR_results_FINAL$representation<-NA
  
  edgeR_results_FINAL$representation[which(edgeR_results_FINAL$comparison%in%plot_together)]<-"plot_together"
  edgeR_results_FINAL$representation[which(edgeR_results_FINAL$comparison%in%plot_ctrl)]<-"plot_ctrl"
  
  edgeR_results_FINAL$representation<-factor(edgeR_results_FINAL$representation,
                                             levels=c("plot_together","plot_ctrl"),
                                             ordered=T)
  
  cat("edgeR_results_FINAL_1\n")
  cat(str(edgeR_results_FINAL))
  cat("\n")
  
  
  
  #### PLOT MSig_DB ----
  
  
  array_representations<-levels(edgeR_results_FINAL$representation)
  
  DEBUG<-0
  
  path_graphs<-paste(out,'Dotplots','/','gene_sets','/',sep='')
  
  if (file.exists(path_graphs)){
    
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
 
  
  
  
  for(i in 1:length(array_representations))
  {
    array_representations_sel<-array_representations[i]
    
    cat("------------------------>\t")
    cat(sprintf(as.character((i))))
    cat("\t")
    cat(sprintf(as.character((array_representations_sel))))
    cat("\n")
    
    MySigDB_results_FINAL_sel<-droplevels(unique(MySigDB_results_FINAL[which(MySigDB_results_FINAL$representation == array_representations_sel),]))
    
    if(DEBUG == 1)
    {
      cat("MySigDB_results_FINAL_sel_0\n")
      cat(str(MySigDB_results_FINAL_sel))
      cat("\n")
    }
    
    
    term_array<-unique(MySigDB_results_FINAL_sel$term.id)
    
    for(k in 1:length(term_array))
    {
      term_array_sel<-term_array[k]
      
      cat("--TERM--->\t")
      cat(sprintf(as.character((k))))
      cat("\t")
      cat(sprintf(as.character((term_array_sel))))
      cat("\n")
      
      M.REP<-unique(MySigDB_results_FINAL_sel[which(MySigDB_results_FINAL_sel$term.id == term_array_sel),])
      
      if(DEBUG == 1)
      {
        cat("M.REP_0\n")
        cat(str(M.REP))
        cat("\n")
      }
      
      
      candidate_genes<-unique(as.character(unlist(strsplit(M.REP$overlap, split=";"))))
      
      if(DEBUG == 1)
      {
        cat("candidate_genes_0\n")
        cat(str(candidate_genes))
        cat("\n")
      }
      
      if(length(candidate_genes) > 42)
      {
        
        candidate_genes<-candidate_genes[1:42]
      }
      
    
      REP_sel<-droplevels(edgeR_results_FINAL[which(edgeR_results_FINAL$Symbol%in%candidate_genes &
                                                      edgeR_results_FINAL$representation==  array_representations_sel),])
      REP_sel$gene_set<-term_array_sel
      
      if(DEBUG == 1)
      {
        cat("REP_sel_0\n")
        cat(str(REP_sel))
        cat("\n")
        cat(str(unique(REP_sel$Symbol)))
        cat("\n")
      }
      
      
     
      REP_sel$Significance<-NA
      
      REP_sel$Significance[which(REP_sel$Minus_logpval >= 1.3)]<-'YES'
      REP_sel$Significance[which(REP_sel$Minus_logpval < 1.3)]<-'NO'
      
      REP_sel$Significance<-factor(REP_sel$Significance,
                                          levels=c('NO','YES'))
      
      if(DEBUG == 1)
      {
        cat("REP_sel_0\n")
        cat(str(REP_sel))
        cat("\n")
      }
      
      ### graph parameters_FC
      
      indx_FC<-which(colnames(REP_sel) == 'logFC')
      
      A_FC<-summary(REP_sel[,indx_FC])
      
      
      if(DEBUG == 1)
      {
        cat("A_FC\n")
        cat(sprintf(as.character(names(A_FC))))
        cat("\n")
        cat(sprintf(as.character(A_FC)))
        cat("\n")
      }
      
      max_value<-A_FC[6]
      min_value<-A_FC[1]
      
      
      step<-round(abs(max_value-min_value)/5,0)
      
      if(step == 0)
      {
        
        step<-1
      }
      breaks_FC<-unique(sort(unique(c(0,seq(min_value,max_value+step, by=step)))))
      labels_FC<-as.character(round(breaks_FC),1)
      
      if(DEBUG == 1)
      {
        cat("step_FC\n")
        cat(sprintf(as.character(step)))
        cat("\n")
        cat("labels_FC\n")
        cat(sprintf(as.character(labels_FC)))
        cat("\n")
      }
      
      ### graph parameters_Minus_logpval
      
      indx_Minus_logpval<-which(colnames(REP_sel) == 'Minus_logpval')
      
      A_Minus_logpval<-summary(REP_sel[,indx_Minus_logpval])
      
      
      if(DEBUG == 1)
      {
        cat("A_Minus_logpval\n")
        cat(sprintf(as.character(names(A_Minus_logpval))))
        cat("\n")
        cat(sprintf(as.character(A_Minus_logpval)))
        cat("\n")
      }
      
      max_value<-A_Minus_logpval[6]
      min_value<-A_Minus_logpval[1]
      
      
      step<-round(abs(max_value-min_value)/5,0)
      
      if(step == 0)
      {
        
        step<-1
      }
      breaks_Minus_logpval<-unique(sort(unique(c(0,seq(min_value,max_value+step, by=step)))))
      labels_Minus_logpval<-as.character(round(breaks_Minus_logpval),1)
      
      if(DEBUG == 1)
      {
        cat("step_Minus_logpval\n")
        cat(sprintf(as.character(step)))
        cat("\n")
        cat("labels_Minus_logpval\n")
        cat(sprintf(as.character(labels_Minus_logpval)))
        cat("\n")
      }
      
      REP_sel<-REP_sel[order(REP_sel$logFC, decreasing=T),]
      
      #### dotplot
      
      if(DEBUG == 1)
      {
        cat("Graph_Part_START:\n")
        
      }
      
      
      
      dotplot<-ggplot(data=REP_sel,
                      aes(y=Symbol,
                          x=time_point)) +
        geom_point(aes(color=Significance,
                       fill=logFC,
                       size=Minus_logpval), stroke=1, shape=21)+
        scale_color_manual(values=c("gray","black"),name='Significant', drop=F)+
        scale_size(range = c(0,10), name='-log10pval',
                   breaks=breaks_Minus_logpval, labels=labels_Minus_logpval, limits=c(breaks_Minus_logpval[1],breaks_Minus_logpval[length(breaks_Minus_logpval)]))+
        scale_fill_gradient2(
          low = "blue", 
          mid = "white", 
          high = "red", 
          midpoint = 0,
          breaks=breaks_FC,labels=labels_FC,
          limits=c(breaks_FC[1]-0.01,breaks_FC[length(breaks_FC)]+0.01),name=paste('LogFC',sep="\n"),na.value = "gray")+
        scale_y_discrete(name=NULL, drop=F)+
        theme_classic()+
        ggeasy::easy_center_title()
      
      if(DEBUG == 1)
      {
        cat("Graph_Part_MIDDLE:\n")
        
        cat("REP_sel_MIDDLE\n")
        cat(str(REP_sel))
        cat("\n")
        
      }
      
      dotplot<-dotplot+
        facet_grid(. ~ comparison, scales='free_x', space='free_x',
                   drop=F) +
        theme_cowplot(font_size = 10)+
        theme( strip.background = element_blank(),
               strip.placement = "inside",
               strip.text = element_text(size=10, angle=0),
               panel.spacing = unit(0.2, "lines"),
               panel.background=element_rect(fill="white"),
               panel.border=element_rect(colour="black",size=1),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())+
        theme(axis.title.y=element_blank(),
              axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
              axis.text.x=element_text(angle=45,size=14, vjust=1, hjust=1,color="black", family="sans"))+
        theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
              legend.key.height = unit(1.5, 'cm'), #change legend key height
              legend.key.width = unit(1, 'cm'), #change legend key width
              legend.title = element_text(size=14), #change legend title font size
              legend.text = element_text(size=14))+ #change legend text font size
        scale_x_discrete(name=NULL, drop=T)
      
      if(DEBUG == 1)
      {
        cat("Graph_Part_END_I:\n")
        
      }
      
      
      
      
      path_graphs<-paste(out,'Dotplots','/','gene_sets','/',array_representations_sel,'/',sep='')
      
      if (file.exists(path_graphs)){
        
        
      }else{
        
        dir.create(file.path(path_graphs))
        
      }#path_graphs
      
      setwd(path_graphs)
      
      svgname<-paste(paste("Dotplot",term_array_sel,array_representations_sel,sep='_'),".svg",sep='')
      makesvg = TRUE
      
      if (makesvg == TRUE)
      {
        ggsave(svgname, plot= dotplot,
               device="svg",
               height=10, width=12)
      }
      
      if(DEBUG == 1)
      {
        cat("Graph_Part_END_II:\n")
        quit(status = 1)
      }
    }#k in 1:dim(MySigDB_results_FINAL_sel)[1]
  }#i in 1:length(LogFC_indx)
}

Dot_plot_figure_graph = function(option_list)
{
  
  suppressMessages(library("ggrepel", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  
  
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
  
  #### READ and transform selected_genes ----
  
  selected_genes = unlist(strsplit(opt$selected_genes, split=","))
  
  cat("selected_genes_\n")
  cat(sprintf(as.character(selected_genes)))
  cat("\n")
  
  #### READ and transform RMV_pathways ----
  
  RMV_pathways = unlist(strsplit(opt$RMV_pathways, split=","))
  
  cat("RMV_pathways_\n")
  cat(sprintf(as.character(RMV_pathways)))
  cat("\n")
  
  #### READ and transform KEEP_pathways ----
  
  KEEP_pathways = unlist(strsplit(opt$KEEP_pathways, split=","))
  
  cat("KEEP_pathways_\n")
  cat(sprintf(as.character(KEEP_pathways)))
  cat("\n")
  
  #### READ and transform plot_together ----
  
  plot_together = unlist(strsplit(opt$plot_together, split=","))
  
  cat("plot_together_\n")
  cat(sprintf(as.character(plot_together)))
  cat("\n")
  
  #### READ and transform plot_ctrl ----
  
  plot_ctrl = unlist(strsplit(opt$plot_ctrl, split=","))
  
  cat("plot_ctrl_\n")
  cat(sprintf(as.character(plot_ctrl)))
  cat("\n")
  
  #### Read the edgeR results file ----
  
  MySigDB_results_FINAL<-readRDS(file=opt$MySigDB_results_FINAL)
  
  MySigDB_results_FINAL$representation<-NA
  
  MySigDB_results_FINAL$representation[which(MySigDB_results_FINAL$comparison%in%plot_together)]<-"plot_together"
  MySigDB_results_FINAL$representation[which(MySigDB_results_FINAL$comparison%in%plot_ctrl)]<-"plot_ctrl"
  
  MySigDB_results_FINAL$representation<-factor(MySigDB_results_FINAL$representation,
                                               levels=c("plot_together","plot_ctrl"),
                                               ordered=T)
  
  
  cat("MySigDB_results_FINAL_0\n")
  cat(str(MySigDB_results_FINAL))
  cat("\n")
  
  #### Read the edgeR results file ----
  
  edgeR_results_FINAL<-readRDS(file=opt$edgeR_results_FINAL)
  
  cat("edgeR_results_FINAL_0\n")
  cat(str(edgeR_results_FINAL))
  cat("\n")
  
  edgeR_results_FINAL$representation<-NA
  
  edgeR_results_FINAL$representation[which(edgeR_results_FINAL$comparison%in%plot_together)]<-"plot_together"
  edgeR_results_FINAL$representation[which(edgeR_results_FINAL$comparison%in%plot_ctrl)]<-"plot_ctrl"
  
  edgeR_results_FINAL$representation<-factor(edgeR_results_FINAL$representation,
                                             levels=c("plot_together","plot_ctrl"),
                                             ordered=T)
  
  cat("edgeR_results_FINAL_1\n")
  cat(str(edgeR_results_FINAL))
  cat("\n")
  
  #### ----
  
  REP_sel<-droplevels(unique(edgeR_results_FINAL[which(edgeR_results_FINAL$Symbol%in%selected_genes),]))
  
  REP_sel$Symbol2<-factor(REP_sel$Symbol,
                          levels=rev(selected_genes),
                          ordered=T)
  
  cat("REP_sel_0\n")
  cat(str(REP_sel))
  cat("\n")
  
  
  #### PLOT MSig_DB ----
  
  
  array_representations<-levels(edgeR_results_FINAL$representation)
  
  DEBUG<-1
  
  for(i in 1:length(array_representations))
  {
    array_representations_sel<-array_representations[i]
    
    cat("------------------------>\t")
    cat(sprintf(as.character((i))))
    cat("\t")
    cat(sprintf(as.character((array_representations_sel))))
    cat("\n")
    
    REP_sel_subset<-droplevels(unique(REP_sel[which(REP_sel$representation == array_representations_sel),]))
    
    REP_sel_subset$Significance<-NA
    
    REP_sel_subset$Significance[which(REP_sel_subset$Minus_logpval >= 1.3)]<-'YES'
    REP_sel_subset$Significance[which(REP_sel_subset$Minus_logpval < 1.3)]<-'NO'
    
    REP_sel_subset$Significance<-factor(REP_sel_subset$Significance,
                                        levels=c('NO','YES'))
    
    if(DEBUG == 1)
    {
      cat("REP_sel_subset_0\n")
      cat(str(REP_sel_subset))
      cat("\n")
    }
    
    ### graph parameters_FC
    
    indx_FC<-which(colnames(REP_sel_subset) == 'logFC')
    
    A_FC<-summary(REP_sel_subset[,indx_FC])
    
    
    if(DEBUG == 1)
    {
      cat("A_FC\n")
      cat(sprintf(as.character(names(A_FC))))
      cat("\n")
      cat(sprintf(as.character(A_FC)))
      cat("\n")
    }
    
    max_value<-A_FC[6]
    min_value<-A_FC[1]
    
    
    step<-round(abs(max_value-min_value)/5,1)
    
    if(step == 0)
    {
      
      step<-1
    }
    breaks_FC<-unique(sort(unique(c(0,seq(min_value,max_value+step, by=step)))))
    labels_FC<-as.character(round(breaks_FC),2)
    
    if(DEBUG == 1)
    {
      cat("step_FC\n")
      cat(sprintf(as.character(step)))
      cat("\n")
      cat("labels_FC\n")
      cat(sprintf(as.character(labels_FC)))
      cat("\n")
    }
    
    ### graph parameters_Minus_logpval
    
    indx_Minus_logpval<-which(colnames(REP_sel_subset) == 'Minus_logpval')
    
    A_Minus_logpval<-summary(REP_sel_subset[,indx_Minus_logpval])
    
    
    if(DEBUG == 1)
    {
      cat("A_Minus_logpval\n")
      cat(sprintf(as.character(names(A_Minus_logpval))))
      cat("\n")
      cat(sprintf(as.character(A_Minus_logpval)))
      cat("\n")
    }
    
    max_value<-A_Minus_logpval[6]
    min_value<-A_Minus_logpval[1]
    
    
    step<-round(abs(max_value-min_value)/5,1)
    
    if(step == 0)
    {
      
      step<-1
    }
    breaks_Minus_logpval<-unique(sort(unique(c(0,seq(min_value,max_value+step, by=step)))))
    labels_Minus_logpval<-as.character(round(breaks_Minus_logpval),2)
    
    if(DEBUG == 1)
    {
      cat("step_Minus_logpval\n")
      cat(sprintf(as.character(step)))
      cat("\n")
      cat("labels_Minus_logpval\n")
      cat(sprintf(as.character(labels_Minus_logpval)))
      cat("\n")
    }
    
    
    #### dotplot
    
    if(DEBUG == 1)
    {
      cat("Graph_Part_START:\n")
      
    }
    
    dotplot<-ggplot(data=REP_sel_subset,
                    aes(y=Symbol2,
                        x=time_point)) +
      geom_point(aes(color=Significance,
                     fill=logFC,
                     size=Minus_logpval), stroke=1, shape=21)+
      scale_color_manual(values=c("black","black"),name='Significant', drop=F)+
      scale_size(range = c(0,8), name='-log10pval',
                 breaks=breaks_Minus_logpval, labels=labels_Minus_logpval, limits=c(breaks_Minus_logpval[1],breaks_Minus_logpval[length(breaks_Minus_logpval)]))+
      scale_fill_gradient2(
        low = "blue", 
        mid = "white", 
        high = "red", 
        midpoint = 0,
        breaks=breaks_FC,labels=labels_FC,
        limits=c(breaks_FC[1]-0.01,breaks_FC[length(breaks_FC)]+0.01),name=paste('LogFC',sep="\n"),na.value = "gray")+
      scale_y_discrete(name=NULL, drop=F)+
      theme_classic()+
      ggeasy::easy_center_title()
    
    if(DEBUG == 1)
    {
      cat("Graph_Part_MIDDLE:\n")
      
      cat("REP_sel_subset_MIDDLE\n")
      cat(str(REP_sel_subset))
      cat("\n")
      
    }
    
    dotplot<-dotplot+
      facet_grid(. ~ comparison, scales='free_x', space='free_x',
                 drop=F) +
      theme_cowplot(font_size = 10)+
      theme( strip.background = element_blank(),
             strip.placement = "inside",
             strip.text = element_text(size=10, angle=0),
             panel.spacing = unit(0.2, "lines"),
             panel.background=element_rect(fill="white"),
             panel.border=element_rect(colour="black",size=1),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())+
      theme(axis.title.y=element_blank(),
            axis.text.y=element_text(angle=0,size=16, color="black", family="sans"),
            axis.text.x=element_text(angle=45,size=14, vjust=1, hjust=1,color="black", family="sans"))+
      theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
            legend.key.height = unit(1.5, 'cm'), #change legend key height
            legend.key.width = unit(1, 'cm'), #change legend key width
            legend.title = element_text(size=14), #change legend title font size
            legend.text = element_text(size=14))+ #change legend text font size
      scale_x_discrete(name=NULL, drop=T)
    
    if(DEBUG == 1)
    {
      cat("Graph_Part_END_I:\n")
      
    }
    
    ##### MySigDB_results_FINAL part----
    
    
    MySigDB_results_FINAL_sel<-droplevels(unique(MySigDB_results_FINAL[which(MySigDB_results_FINAL$representation == array_representations_sel),]))
    
    MySigDB_results_FINAL_sel$Significance<-NA
    
    MySigDB_results_FINAL_sel$Significance[which(MySigDB_results_FINAL_sel$Minus_logpval >= 1.3)]<-'YES'
    MySigDB_results_FINAL_sel$Significance[which(MySigDB_results_FINAL_sel$Minus_logpval < 1.3)]<-'NO'
    
    MySigDB_results_FINAL_sel$Significance<-factor(MySigDB_results_FINAL_sel$Significance,
                                        levels=c('NO','YES'))
    
    if(DEBUG == 1)
    {
      cat("MySigDB_results_FINAL_sel_0\n")
      cat(str(MySigDB_results_FINAL_sel))
      cat("\n")
    }
    
    indx.RMV<-grep(paste(RMV_pathways,collapse="|"),MySigDB_results_FINAL_sel$term.id)
    
    if(DEBUG == 1)
    {
      cat("indx.RMV_0\n")
      cat(str(indx.RMV))
      cat("\n")
    }
    
    
    
    MySigDB_results_FINAL_sel<-unique(MySigDB_results_FINAL_sel[-indx.RMV,])
    
    if(DEBUG == 1)
    {
      cat("MySigDB_results_FINAL_sel_1\n")
      cat(str(MySigDB_results_FINAL_sel))
      cat("\n")
    }
    
    indx.display<-grep(paste(KEEP_pathways,collapse="|"),MySigDB_results_FINAL_sel$term.id)
    
    if(DEBUG == 1)
    {
      cat("indx.display_0\n")
      cat(str(indx.display))
      cat("\n")
    }
    
    MySigDB_results_FINAL_sel<-unique(MySigDB_results_FINAL_sel[indx.display,])
    
    
    # MySigDB_results_FINAL_sel$comparison<-factor(MySigDB_results_FINAL_sel$comparison,
    #                                                levels=rev(levels(MySigDB_results_FINAL_sel$comparison)),
    #                                                ordered=T)
    
    if(DEBUG == 1)
    {
      cat("MySigDB_results_FINAL_sel_2\n")
      cat(str(MySigDB_results_FINAL_sel))
      cat("\n")
    }
    
    
    
    MySigDB_results_FINAL_sel.dt<-data.table(MySigDB_results_FINAL_sel, key=c("term.id","comparison"))
    
    
    MySigDB_results_FINAL_sel_2<-as.data.frame(MySigDB_results_FINAL_sel.dt[,.(n_genes_in_overlap=length(unlist(strsplit(overlap, split=';'))),
                                                                              Minus_logpval=Minus_logpval,
                                                                              overlap=overlap,
                                                                              Significance=Significance), by=key(MySigDB_results_FINAL_sel.dt)], stringsAsFactors=F)
    
    
    
    MySigDB_results_FINAL_sel_2$log2_n_overlap<-log2(MySigDB_results_FINAL_sel_2$n_genes_in_overlap)
    levels_comp<-levels(MySigDB_results_FINAL_sel_2$comparison)
    
    MySigDB_results_FINAL_sel_2$comparison<-factor(MySigDB_results_FINAL_sel_2$comparison,
                                                  levels=rev(levels_comp),
                                                  ordered=T)
    
    if(DEBUG == 1)
    {
      cat("MySigDB_results_FINAL_sel_2_0\n")
      cat(str(MySigDB_results_FINAL_sel_2))
      cat("\n")
    }
    
    ### graph parameters_n_overlap
    
    indx_log2_n_overlap<-which(colnames(MySigDB_results_FINAL_sel_2) == 'log2_n_overlap')
    
    A_log2_n_overlap<-summary(MySigDB_results_FINAL_sel_2[,indx_log2_n_overlap])
    
    
    if(DEBUG == 1)
    {
      cat("A_log2_n_overlap\n")
      cat(sprintf(as.character(names(A_log2_n_overlap))))
      cat("\n")
      cat(sprintf(as.character(A_log2_n_overlap)))
      cat("\n")
    }
    
    max_value<-A_log2_n_overlap[6]
    min_value<-A_log2_n_overlap[1]
    
    
    step<-round(abs(max_value-min_value)/5,0)
    
    if(step == 0)
    {
      
      step<-1
    }
    breaks_log2_n_overlap<-unique(sort(unique(c(seq(min_value,max_value+step, by=step)))))
    labels_log2_n_overlap<-as.character(round(2^breaks_log2_n_overlap),0)
    
    if(DEBUG == 1)
    {
      cat("step_log2_n_overlap\n")
      cat(sprintf(as.character(step)))
      cat("\n")
      cat("labels_log2_n_overlap\n")
      cat(sprintf(as.character(labels_log2_n_overlap)))
      cat("\n")
    }
    
    ### graph parameters_Minus_logpval
    
    indx_Minus_logpval<-which(colnames(MySigDB_results_FINAL_sel_2) == 'Minus_logpval')
    
    A_Minus_logpval<-summary(MySigDB_results_FINAL_sel_2[,indx_Minus_logpval])
    
    
    if(DEBUG == 1)
    {
      cat("A_Minus_logpval\n")
      cat(sprintf(as.character(names(A_Minus_logpval))))
      cat("\n")
      cat(sprintf(as.character(A_Minus_logpval)))
      cat("\n")
    }
    
    max_value<-A_Minus_logpval[6]
    min_value<-A_Minus_logpval[1]
    
    
    step<-round(abs(max_value-min_value)/3,1)
    
    if(step == 0)
    {
      
      step<-1
    }
    breaks_Minus_logpval<-unique(sort(unique(c(seq(min_value,max_value+step, by=step)))))
    labels_Minus_logpval<-as.character(round(breaks_Minus_logpval),2)
    
    if(DEBUG == 1)
    {
      cat("step_Minus_logpval\n")
      cat(sprintf(as.character(step)))
      cat("\n")
      cat("labels_Minus_logpval\n")
      cat(sprintf(as.character(labels_Minus_logpval)))
      cat("\n")
    }
    
    jitter_pos <- position_jitter(width=0.5, height = 0.5, seed = 1)
    
    if(DEBUG == 1)
    {
      cat("jitter_pos\n")
      cat(str(jitter_pos))
      cat("\n")
    }
    
    vector_colors<-brewer.pal(3, "Dark2")[c(3,2,1)]
    
    if(DEBUG == 1)
    {
      cat("vector_colors\n")
      cat(str(vector_colors))
      cat("\n")
    }
    
    #### dotplot
    
    if(DEBUG == 1)
    {
      cat("Gene_set_dotplot_START:\n")
      
    }
    
   
    Gene_set_dotplot<-ggplot(data=MySigDB_results_FINAL_sel_2,
                    aes(y=as.numeric(comparison),
                        x=log2_n_overlap,
                        color=comparison))+
      geom_point(aes(size=Minus_logpval, fill=comparison),position=jitter_pos, stroke=1, shape=21)+
      scale_size(range = c(0,8), name='-log10pval',
                 breaks=breaks_Minus_logpval, labels=labels_Minus_logpval, limits=c(breaks_Minus_logpval[1]-0.1,breaks_Minus_logpval[length(breaks_Minus_logpval)]+0.1))+
      scale_y_continuous(name=NULL, breaks=NULL)+
      scale_x_continuous(name="n_genes_overlap_in_gene_set",
                         breaks=breaks_log2_n_overlap,
                         labels=labels_log2_n_overlap,
                         limits=c(breaks_log2_n_overlap[1],breaks_log2_n_overlap[length(breaks_log2_n_overlap)]))+
      scale_color_manual(values=vector_colors, drop=F)+
      scale_fill_manual(values=vector_colors, drop=F)+
      geom_text_repel(aes(label = term.id), 
                      position=jitter_pos,
                      family = "sans",
                      size = 4,
                      min.segment.length = 0, 
                      seed = 42, 
                      box.padding = 0.5,
                      arrow = arrow(length = unit(0.015, "npc")),
                      max.overlaps = Inf)+
      theme_classic()+
      theme(axis.title.y=element_blank(),
            axis.title.x=element_text(size=16, color="black", family="sans"),
            axis.text.y=element_blank(),
            axis.text.x=element_text(angle=0,size=14, color="black", family="sans"))+
      theme(legend.position="right",
            legend.key.size = unit(1.5, 'cm'), #change legend key size
            legend.key.height = unit(1.5, 'cm'), #change legend key height
            legend.key.width = unit(1, 'cm'), #change legend key width
            legend.title = element_text(size=16, family ="sans"), #change legend title font size
            legend.text = element_text(size=14, family ="sans"))+ #change legend text font size
      ggeasy::easy_center_title()
    
    if(DEBUG == 1)
    {
      cat("Gene_set_dotplot_END:\n")
    }
    
    
   
    graph_DEF<-plot_grid(dotplot,Gene_set_dotplot,
                         ncol = 2,
                         nrow=1,
                         rel_widths=c(1,1))
    
    path_graphs<-paste(out,'Dotplots','/',sep='')
    
    if (file.exists(path_graphs)){
      
      
    }else{
      
      dir.create(file.path(path_graphs))
      
    }#path_graphs
    
    setwd(path_graphs)
    
    svgname<-paste(paste("Dotplot",'Figure',paste(array_representations_sel, collapse="__"), sep='_'),".svg",sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      ggsave(svgname, plot= graph_DEF,
             device="svg",
             height=10, width=12)
    }
    
    if(DEBUG == 1)
    {
      cat("THE_END:\n")
      
    }
  }#i in 1:length(LogFC_indx)
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
    make_option(c("--MySigDB_results_FINAL"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--edgeR_results_FINAL"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ctrl_gene"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--plot_together"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--plot_ctrl"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--selected_genes"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--RMV_pathways"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--KEEP_pathways"), type="character", default=NULL, 
                metavar="filename", 
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
  
  
  # Dot_plot_function_CTRL_gene(opt)
  # Dot_plot_function_MSigDB(opt)
  Dot_plot_figure_graph(opt)
  
}


###########################################################################

system.time( main() )