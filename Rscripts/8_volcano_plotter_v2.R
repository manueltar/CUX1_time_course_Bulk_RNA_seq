
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
suppressMessages(library("clusterProfiler", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("enrichplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("DOSE", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

opt = NULL

# options(warn = 1)



volcano_plotter = function(option_list)
{
  suppressMessages(library("edgeR", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
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
  
  #### Read the edgeR results file ----
  
  MySigDB_results<-as.data.frame(fread(file=opt$MySigDB_results, sep="\t", header=T), stringsAsFactors=F)
  
  cat("MySigDB_results_0\n")
  cat(str(MySigDB_results))
  cat("\n")
  
  MySigDB_results_subset<-unique(MySigDB_results[,-which(colnames(MySigDB_results) == 'collection')])
  
  cat("MySigDB_results_subset_0\n")
  cat(str(MySigDB_results_subset))
  cat("\n")
  
  MySigDB_results_subset.dt<-data.table(MySigDB_results_subset, key=c('term.id','term.name'))
  
  
  MySigDB_results_subset_MIN<-as.data.frame(MySigDB_results_subset.dt[,.SD[which.min(adjusted.p.val)], by=key(MySigDB_results_subset.dt)], stringsAsFactos=F)
  
  cat("MySigDB_results_subset_MIN_0\n")
  cat(str(MySigDB_results_subset_MIN))
  cat("\n")
  
  
  #### Read the edgeR results file ----
  
  edgeR_results<-as.data.frame(fread(file=opt$edgeR_results, sep="\t", header=T), stringsAsFactors=F)
  
  cat("edgeR_results_0\n")
  cat(str(edgeR_results))
  cat("\n")
  
  
  LogFC_indx<-grep("logFC\\.",colnames(edgeR_results))
  
  cat("LogFC_indx_0\n")
  cat(sprintf(as.character((LogFC_indx))))
  cat("\n")
  
  candidate_colnames<-colnames(edgeR_results)[LogFC_indx]
  
  cat("candidate_colnames_0\n")
  cat(str(candidate_colnames))
  cat("\n")
  
  DEBUG<-0
  
  for(i in 1:length(candidate_colnames))
  {
    candidate_colnames_sel<-candidate_colnames[i]
    
    cat("------------------------>\t")
    cat(sprintf(as.character((i))))
    cat("\t")
    cat(sprintf(as.character((candidate_colnames_sel))))
    cat("\t")
    
    LogFC_sel_2<-gsub("logFC\\.","",candidate_colnames_sel)
    
    cat(sprintf(as.character((LogFC_sel_2))))
    cat("\n")
    
    indx.int<-c(which(colnames(edgeR_results) == candidate_colnames_sel),
                which(colnames(edgeR_results) == 'FDR'),
                which(colnames(edgeR_results) == 'Symbol'),
                which(colnames(edgeR_results) == 'ensembl_gene_id'))
    
    if(DEBUG == 1)
    {
      cat("indx.int\n")
      cat(sprintf(as.character((indx.int))))
      cat("\n")

    }

    edgeR_results_subset<-unique(edgeR_results[,indx.int])

    if(DEBUG == 1)
    {
      cat("edgeR_results_subset_0\n")
      cat(str(edgeR_results_subset))
      cat("\n")
    }

    edgeR_results_subset$ENTREZID_id = mapIds(org.Hs.eg.db,
                                   keys=edgeR_results_subset$ensembl_gene_id,
                                   column="ENTREZID",
                                   keytype="ENSEMBL",
                                   multiVals="first")

    if(DEBUG == 1)
    {
      cat("edgeR_results_subset_1\n")
      cat(str(edgeR_results_subset))
      cat("\n")
    }

    edgeR_results_subset.NO.NA<-edgeR_results_subset[!is.na(edgeR_results_subset$ENTREZID_id),]
    edgeR_results_subset.NO.NA$Minus_logpval<--1*log10(edgeR_results_subset.NO.NA$FDR)

    if(DEBUG == 1)
    {
      cat("edgeR_results_subset.NO.NA_0\n")
      cat(str(edgeR_results_subset.NO.NA))
      cat("\n")
    }
    
    indx.FC<-which(colnames(edgeR_results_subset.NO.NA) == candidate_colnames_sel)
    indx.selection_UP<-which(edgeR_results_subset.NO.NA$Minus_logpval >=1.3 &
                            edgeR_results_subset.NO.NA[,indx.FC] > 0)
    indx.selection_DOWN<-which(edgeR_results_subset.NO.NA$Minus_logpval >=1.3 &
                               edgeR_results_subset.NO.NA[,indx.FC] < 0)
    
    
    selected_genes<-c("ENSG00000005961","ENSG00000257923","ENSG00000259207","ENSG00000170180")
    indx.selected_genes<-which(edgeR_results_subset.NO.NA$ensembl_gene_id%in%selected_genes)
    
    if(DEBUG == 1)
    {
      cat("indx.FC_0\n")
      cat(str(indx.FC))
      cat("\n")
      cat("indx.selection_UP_0\n")
      cat(str(indx.selection_UP))
      cat("\n")
      cat("indx.selection_DOWN_0\n")
      cat(str(indx.selection_DOWN))
      cat("\n")
      cat("indx.selected_genes_0\n")
      cat(str(indx.selected_genes))
      cat("\n")
    }
    
    edgeR_results_subset.NO.NA$Volcano_label<-""
    
    indx.labelling<-which(edgeR_results_subset.NO.NA$Minus_logpval >=3)
    
    edgeR_results_subset.NO.NA$Volcano_label[indx.labelling]<-edgeR_results_subset.NO.NA$Symbol[indx.labelling]

    if(DEBUG == 1)
    {
      cat("edgeR_results_subset.NO.NA_1\n")
      cat(str(edgeR_results_subset.NO.NA))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(edgeR_results_subset.NO.NA$Volcano_label))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(edgeR_results_subset.NO.NA$Volcano_label)))))
      cat("\n")
    }
    
    A_MINUS_LOGPVAL<-summary(edgeR_results_subset.NO.NA$Minus_logpval)
    
    
    if(DEBUG == 1)
    {
      cat("A_MINUS_LOGPVAL\n")
      cat(sprintf(as.character(names(A_MINUS_LOGPVAL))))
      cat("\n")
      cat(sprintf(as.character(A_MINUS_LOGPVAL)))
      cat("\n")
    }
    
    max_value<-A_MINUS_LOGPVAL[6]
    min_value<-A_MINUS_LOGPVAL[1]
    
    
    step<-round(abs(max_value-min_value)/10,0)
    
    if(step == 0)
    {
      
      step<-1
    }
    breaks.MINUS_LOGPVAL<-sort(unique(c(0,seq(min_value,max_value+step, by=step))))
    labels.MINUS_LOGPVAL<-as.character(round(breaks.MINUS_LOGPVAL),1)
    
    if(DEBUG == 1)
    {
      cat("step_MINUS_LOGPVAL\n")
      cat(sprintf(as.character(step)))
      cat("\n")
      cat("breaks.MINUS_LOGPVAL\n")
      cat(sprintf(as.character(breaks.MINUS_LOGPVAL)))
      cat("\n")
      cat("labels.MINUS_LOGPVAL\n")
      cat(sprintf(as.character(labels.MINUS_LOGPVAL)))
      cat("\n")
    }
    
    
    A_FC<-summary(edgeR_results_subset.NO.NA[,indx.FC])
    
    
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
    breaks.FC<-sort(unique(c(0,seq(min_value,max_value+step, by=step))))
    labels.FC<-as.character(round(breaks.FC),1)
    
    if(DEBUG == 1)
    {
      cat("step_FC\n")
      cat(sprintf(as.character(step)))
      cat("\n")
      cat("labels.FC\n")
      cat(sprintf(as.character(labels.FC)))
      cat("\n")
    }
    
    if(DEBUG == 1)
    {
      cat("Graph_Part_START:\n")
      
    }
    
    volcano<-ggplot(data=edgeR_results_subset.NO.NA,
                    aes(x=edgeR_results_subset.NO.NA[,indx.FC], 
                        y=Minus_logpval,
                        label=Volcano_label)) +
      geom_point(size=4, color="gray")+
      geom_point(data=edgeR_results_subset.NO.NA[indx.selection_UP,],
                 aes(x=edgeR_results_subset.NO.NA[indx.selection_UP,indx.FC], 
                     y=Minus_logpval),color="red",size=4)+
      geom_point(data=edgeR_results_subset.NO.NA[indx.selection_DOWN,],
                 aes(x=edgeR_results_subset.NO.NA[indx.selection_DOWN,indx.FC], 
                     y=Minus_logpval), color="blue",size=4)+
      geom_text_repel(aes(label = Volcano_label), 
                      box.padding = 1,
                      max.overlaps = Inf,
                      show.legend = FALSE)+
      scale_x_continuous(name=candidate_colnames_sel, 
                         breaks=breaks.FC,
                         labels=labels.FC, 
                         limits=c(breaks.FC[1]-0.1,breaks.FC[length(breaks.FC)]+0.1))+
      scale_y_continuous(name="-log10pval adjusted", 
                         breaks=breaks.MINUS_LOGPVAL,
                         labels=labels.MINUS_LOGPVAL, 
                         limits=c(breaks.MINUS_LOGPVAL[1],breaks.MINUS_LOGPVAL[length(breaks.MINUS_LOGPVAL)]+0.1))+
      theme_bw()+
      theme(axis.title.y=element_text(size=24, family="sans"),
            axis.title.x=element_text(size=24, family="sans"),
            axis.text.y=element_text(angle=0,size=24, color="black", family="sans"),
            axis.text.x=element_text(angle=0, size=24, color="black", family="sans"),
            legend.title=element_text(size=16,color="black", family="sans"),
            legend.text=element_text(size=12,color="black", family="sans"))+
      theme(legend.position="hidden",legend.title=element_blank(), legend.text = element_text(size=10))+
      geom_vline(xintercept=0,linetype="dashed")+
      ggeasy::easy_center_title()
    
    path_graphs<-paste(out,'volcano_plots','/',sep='')
    
    if (file.exists(path_graphs)){
      
      
    }else{
      
      dir.create(file.path(path_graphs))
      
    }#path_graphs
    
    path_graphs<-paste(out,'volcano_plots','/',LogFC_sel_2,'/',sep='')
    
    if (file.exists(path_graphs)){
      
      
    }else{
      
      dir.create(file.path(path_graphs))
      
    }#path_graphs
    
    setwd(path_graphs)
    
    svgname<-paste(paste("Volcano_general",LogFC_sel_2, sep='_'),".svg",sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      ggsave(svgname, plot= volcano,
             device="svg",
             height=10, width=12)
    }

    if(DEBUG == 1)
    {
      cat("Graph_Part_END:\n")
      
      

    }
    
   
    
    for(k in 1:dim(MySigDB_results_subset_MIN)[1])
    {
      MySigDB_results_subset_MIN_sel<-MySigDB_results_subset_MIN[k,]
      
      if(DEBUG == 1)
      {
        cat("MySigDB_results_subset_MIN_sel_0\n")
        cat(str(MySigDB_results_subset_MIN_sel))
        cat("\n")
        
      }
      
      gene_set_ID<-unique(MySigDB_results_subset_MIN_sel$term.id)
      
      cat("---gene_set-->\t")
      cat(sprintf(as.character(k)))
      cat("\t")
      cat(sprintf(as.character(gene_set_ID)))
      cat("\n")
      
      
      overlap_Symbol<-unlist(strsplit(MySigDB_results_subset_MIN_sel$overlap, split=';'))
      
      if(DEBUG == 1)
      {
        cat("overlap_Symbol_0\n")
        cat(str(overlap_Symbol))
        cat("\n")
        
      }
      
    
        
        indx.FC<-which(colnames(edgeR_results_subset.NO.NA) == candidate_colnames_sel)
        indx.selection_UP<-which(edgeR_results_subset.NO.NA$Minus_logpval >=1.3 &
                                   edgeR_results_subset.NO.NA[,indx.FC] > 0)
        indx.selection_DOWN<-which(edgeR_results_subset.NO.NA$Minus_logpval >=1.3 &
                                     edgeR_results_subset.NO.NA[,indx.FC] < 0)
        
        
      
        
        if(DEBUG == 1)
        {
          cat("indx.FC_0\n")
          cat(str(indx.FC))
          cat("\n")
          cat("indx.selection_UP_0\n")
          cat(str(indx.selection_UP))
          cat("\n")
          cat("indx.selection_DOWN_0\n")
          cat(str(indx.selection_DOWN))
          cat("\n")
         
        }
        
        edgeR_results_subset.NO.NA_sel<-edgeR_results_subset.NO.NA[which(edgeR_results_subset.NO.NA$Symbol%in%overlap_Symbol),]
        
        if(DEBUG == 1)
        {
          cat("edgeR_results_subset.NO.NA_sel_0\n")
          cat(str(edgeR_results_subset.NO.NA_sel))
          cat("\n")
          
        }
        
        indx.FC<-which(colnames(edgeR_results_subset.NO.NA_sel) == candidate_colnames_sel)
        indx.selection_UP<-which(edgeR_results_subset.NO.NA_sel$Minus_logpval >=1.3 &
                                   edgeR_results_subset.NO.NA_sel[,indx.FC] > 0)
        indx.selection_DOWN<-which(edgeR_results_subset.NO.NA_sel$Minus_logpval >=1.3 &
                                     edgeR_results_subset.NO.NA_sel[,indx.FC] < 0)
        
        
        
        if(DEBUG == 1)
        {
          cat("indx.FC_0\n")
          cat(str(indx.FC))
          cat("\n")
          cat("indx.selection_UP_0\n")
          cat(str(indx.selection_UP))
          cat("\n")
          cat("indx.selection_DOWN_0\n")
          cat(str(indx.selection_DOWN))
          cat("\n")
         
        }
        
        
        edgeR_results_subset.NO.NA_sel$Volcano_label<-""
        
        indx.labelling<-which(edgeR_results_subset.NO.NA_sel$Symbol%in%overlap_Symbol)
        
        edgeR_results_subset.NO.NA_sel$Volcano_label[indx.labelling]<-edgeR_results_subset.NO.NA_sel$Symbol[indx.labelling]
        
        if(DEBUG == 1)
        {
          cat("edgeR_results_subset.NO.NA_sel_1\n")
          cat(str(edgeR_results_subset.NO.NA_sel))
          cat("\n")
          cat(sprintf(as.character(names(summary(as.factor(edgeR_results_subset.NO.NA_sel$Volcano_label))))))
          cat("\n")
          cat(sprintf(as.character(summary(as.factor(edgeR_results_subset.NO.NA_sel$Volcano_label)))))
          cat("\n")
        }
        
        A_MINUS_LOGPVAL<-summary(edgeR_results_subset.NO.NA_sel$Minus_logpval)
        
        
        if(DEBUG == 1)
        {
          cat("A_MINUS_LOGPVAL\n")
          cat(sprintf(as.character(names(A_MINUS_LOGPVAL))))
          cat("\n")
          cat(sprintf(as.character(A_MINUS_LOGPVAL)))
          cat("\n")
        }
        
        max_value<-A_MINUS_LOGPVAL[6]
        min_value<-A_MINUS_LOGPVAL[1]
        
        
        step<-round(abs(max_value-min_value)/10,0)
        
        if(step == 0)
        {
          
          step<-1
        }
        breaks.MINUS_LOGPVAL<-sort(unique(c(0,seq(min_value,max_value+step, by=step))))
        labels.MINUS_LOGPVAL<-as.character(round(breaks.MINUS_LOGPVAL),1)
        
        if(DEBUG == 1)
        {
          cat("step_MINUS_LOGPVAL\n")
          cat(sprintf(as.character(step)))
          cat("\n")
          cat("breaks.MINUS_LOGPVAL\n")
          cat(sprintf(as.character(breaks.MINUS_LOGPVAL)))
          cat("\n")
          cat("labels.MINUS_LOGPVAL\n")
          cat(sprintf(as.character(labels.MINUS_LOGPVAL)))
          cat("\n")
        }
        
        
        A_FC<-summary(edgeR_results_subset.NO.NA_sel[,indx.FC])
        
        
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
        breaks.FC<-sort(unique(c(0,seq(min_value,max_value+step, by=step))))
        labels.FC<-as.character(round(breaks.FC),1)
        
        if(DEBUG == 1)
        {
          cat("step_FC\n")
          cat(sprintf(as.character(step)))
          cat("\n")
          cat("labels.FC\n")
          cat(sprintf(as.character(labels.FC)))
          cat("\n")
        }
        
        if(DEBUG == 1)
        {
          cat("Graph_Part_START:\n")
          
        }
        
        volcano<-ggplot(data=edgeR_results_subset.NO.NA_sel,
                        aes(x=edgeR_results_subset.NO.NA_sel[,indx.FC], 
                            y=Minus_logpval,
                            label=Volcano_label)) +
          geom_point(size=4, color="gray")+
          geom_point(data=edgeR_results_subset.NO.NA_sel[indx.selection_UP,],
                     aes(x=edgeR_results_subset.NO.NA_sel[indx.selection_UP,indx.FC], 
                         y=Minus_logpval),color="red",size=4)+
          geom_point(data=edgeR_results_subset.NO.NA_sel[indx.selection_DOWN,],
                     aes(x=edgeR_results_subset.NO.NA_sel[indx.selection_DOWN,indx.FC], 
                         y=Minus_logpval), color="blue",size=4)+
          geom_text_repel(aes(label = Volcano_label), 
                          box.padding = 1,
                          max.overlaps = Inf,
                          show.legend = FALSE,
                          size=6)+
          scale_x_continuous(name=candidate_colnames_sel, 
                             breaks=breaks.FC,
                             labels=labels.FC, 
                             limits=c(breaks.FC[1]-0.1,breaks.FC[length(breaks.FC)]+0.1))+
          scale_y_continuous(name="-log10pval adjusted", 
                             breaks=breaks.MINUS_LOGPVAL,
                             labels=labels.MINUS_LOGPVAL, 
                             limits=c(breaks.MINUS_LOGPVAL[1],breaks.MINUS_LOGPVAL[length(breaks.MINUS_LOGPVAL)]+0.1))+
          theme_bw()+
          theme(axis.title.y=element_text(size=24, family="sans"),
                axis.title.x=element_text(size=24, family="sans"),
                axis.text.y=element_text(angle=0,size=24, color="black", family="sans"),
                axis.text.x=element_text(angle=0, size=24, color="black", family="sans"),
                legend.title=element_text(size=16,color="black", family="sans"),
                legend.text=element_text(size=12,color="black", family="sans"))+
          theme(legend.position="hidden",legend.title=element_blank(), legend.text = element_text(size=10))+
          geom_vline(xintercept=0,linetype="dashed")+
          ggtitle(paste(gene_set_ID, paste("adj.pval:" , round(unique(MySigDB_results_subset_MIN_sel$adjusted.p.val),6), sep=" "),sep="\n"))+
          theme(plot.title = element_text(size =18, face="bold"))+
          ggeasy::easy_center_title()
        
      
        
        path_graphs<-paste(out,'volcano_plots','/',LogFC_sel_2,'/','gene_sets','/',sep='')
        
        if (file.exists(path_graphs)){
          
          
        }else{
          
          dir.create(file.path(path_graphs))
          
        }#path_graphs
        
        setwd(path_graphs)
        
        svgname<-paste(paste("Volcano_gene_set",LogFC_sel_2,gene_set_ID, sep='_'),".svg",sep='')
       
        makesvg = TRUE
        
        if (makesvg == TRUE)
        {
          ggsave(svgname, plot= volcano,
                 device="svg",
                 height=10, width=12)
        }
        
        if(DEBUG == 1)
        {
          cat("Graph_Part_END:\n")
        }
      
    }#k in 1:dim(MySigDB_results_subset_MIN)[1]
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
    make_option(c("--MySigDB_results"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--edgeR_results"), type="character", default=NULL, 
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
  
  
  volcano_plotter(opt)
  
  
}


###########################################################################

system.time( main() )