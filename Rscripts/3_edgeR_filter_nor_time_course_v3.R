
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

filter_and_nor = function(option_list)
{
  suppressMessages(library("edgeR", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  
  
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
  
  #### Read manifest_FINAL_with_results file ----
  
  #"manifest_FINAL_with_results.rds"
  
  manifest_FINAL_with_results<-readRDS(file=opt$manifest_FINAL_with_results)
  

  cat("manifest_FINAL_with_results_0\n")
  cat(str(manifest_FINAL_with_results))
  cat("\n")
  
  #### Read gene_expression_results file ----
  
  gene_expression_results<-as.data.frame(fread(file=opt$gene_expression_results, sep="\t", header=T), stringsAsFactors=F)

  cat("gene_expression_results_0\n")
  cat(str(gene_expression_results))
  cat("\n")
  
  
  #### create the expression matrix ----
  
  gene_expression_matrix<-as.matrix(gene_expression_results[,-which(colnames(gene_expression_results) == 'ensembl_gene_id')])
  row.names(gene_expression_matrix)<-gene_expression_results$ensembl_gene_id
  
  cat("gene_expression_matrix_0\n")
  cat(str(gene_expression_matrix))
  cat("\n")
  
  #### Make the primary edgeR object ----
  
  edgeR_primary <- DGEList(counts=gene_expression_matrix, group=colnames(gene_expression_matrix))
  
  cat("edgeR_primary_0\n")
  cat(str(edgeR_primary))
  cat("\n")

  #### Filter out lowly expressed genes ----
  
  keep <- filterByExpr(edgeR_primary)
  
  cat("keep_0\n")
  cat(str(keep))
  cat("\n")
  
  
  edgeR_primary <- edgeR_primary[keep, , keep.lib.sizes=FALSE]
  
  cat("edgeR_primary_FILTERED\n")
  cat(str(edgeR_primary))
  cat("\n")
  
  #### normalization ----
  
  edgeR_primary<-calcNormFactors(edgeR_primary, method = "TMM")
  
  cat("edgeR_primary_Normalized\n")
  cat(str(edgeR_primary))
  cat("\n")
  
  
  Nor_factors<-edgeR_primary$samples
  
  cat("Nor_factors_Normalized\n")
  cat(str(Nor_factors))
  cat("\n")
  
  Nor_factors_subset<-Nor_factors[,c(which(colnames(Nor_factors) == 'group'),which(colnames(Nor_factors) == 'norm.factors'))]
  Nor_factors_subset[order(Nor_factors_subset$group),]
  
  cat("Nor_factors_subset_0\n")
  cat(str(Nor_factors_subset))
  cat("\n")
  
  
  Nor_factors_subset_matrix<-t(as.matrix(Nor_factors_subset[,-which(colnames(Nor_factors_subset) == 'group')]))
  colnames(Nor_factors_subset_matrix)<-Nor_factors_subset$group
  row.names(Nor_factors_subset_matrix)<-'TMM_factors'
  
  cat("Nor_factors_subset_matrix_0\n")
  cat(str(Nor_factors_subset_matrix))
  cat("\n")
  
  Nor_factors_subset_matrix_inverse<-as.vector(round(1/Nor_factors_subset_matrix,2))
  
  Nor_factors_subset_matrix_inverse<-setNames(Nor_factors_subset_matrix_inverse, colnames(Nor_factors_subset_matrix))
  
  cat("Nor_factors_subset_matrix_inverse_0\n")
  cat(str(Nor_factors_subset_matrix_inverse))
  cat("\n")
  
  #### Recover matrixes after filtering, divide by the TIMM scale factors and output pre-nor and post normalization matrixes for plotting-----
  
  PRE_NOR<-as.matrix(edgeR_primary)
  
  cat("PRE_NOR_0\n")
  cat(str(PRE_NOR))
  cat("\n")
  
  
  
  Nor_factors_subset_matrix_inverse_scaled<-t(matrix(Nor_factors_subset_matrix_inverse  , length(Nor_factors_subset_matrix_inverse) , nrow(PRE_NOR) ))

  colnames(Nor_factors_subset_matrix_inverse_scaled)<-names(Nor_factors_subset_matrix_inverse)
  
  cat("Nor_factors_subset_matrix_inverse_scaled_0\n")
  cat(str(Nor_factors_subset_matrix_inverse_scaled))
  cat("\n")
  
  POST_NOR<-PRE_NOR*Nor_factors_subset_matrix_inverse_scaled
  
  cat("POST_NOR_0\n")
  cat(str(POST_NOR))
  cat("\n")
  
  
  POST_NOR<-apply(POST_NOR, c(1,2), round)
  
  cat("POST_NOR_1\n")
  cat(str(POST_NOR))
  cat("\n")
  
  
  #### graphs ----
  
  path_graphs<-paste(out,'graphs','/',sep='')
  
  if (file.exists(path_graphs)){
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
  ##### Long matrix for PRE-NOR & POST_NOR----
  
  PRE_NOR_df<-as.data.frame(PRE_NOR)
  PRE_NOR_df$ensembl_gene_id<-row.names(PRE_NOR)
  
  cat("PRE_NOR_df_0\n")
  cat(str(PRE_NOR_df))
  cat("\n")
  
  PRE_NOR_LONG<-melt(PRE_NOR_df, id_variable='ensembl_gene_id', variable.name="Sample_label", value.name="counts")
  PRE_NOR_LONG$NOR<-'PRE'
  
  cat("PRE_NOR_LONG_0\n")
  cat(str(PRE_NOR_LONG))
  cat("\n")
  
  POST_NOR_df<-as.data.frame(POST_NOR)
  POST_NOR_df$ensembl_gene_id<-row.names(POST_NOR)
  
  cat("POST_NOR_df_0\n")
  cat(str(POST_NOR_df))
  cat("\n")
  
  POST_NOR_LONG<-melt(POST_NOR_df, id_variable='ensembl_gene_id', variable.name="Sample_label", value.name="counts")
  POST_NOR_LONG$NOR<-'POST'
  
  cat("POST_NOR_LONG_0\n")
  cat(str(POST_NOR_LONG))
  cat("\n")
  
  #### Representation ----
  
  REP<-rbind(PRE_NOR_LONG,POST_NOR_LONG)
  
  REP$NOR<-factor(REP$NOR,
                  levels=c('PRE','POST'),
                  ordered=T)
  
  cat("REP_0\n")
  cat(str(REP))
  cat("\n")
  cat(sprintf(as.character(names(summary(REP$NOR)))))
  cat("\n")
  cat(sprintf(as.character(summary(REP$NOR))))
  cat("\n")
  
  REP<-merge(REP,manifest_FINAL_with_results, by="Sample_label")
  
  cat("REP_1\n")
  cat(str(REP))
  cat("\n")
  
  REP$log2_counts<-log2(REP$counts + 1)
  
  #### violin plots ----
  
  DEBUG<-1
  
  
  A<-round(summary(REP$log2_counts[!is.na(REP$log2_counts)]),2)
  
  
  # cat("summary_GENE_EXP\n")
  # cat(sprintf(as.character(names(A))))
  # cat("\n")
  # cat(sprintf(as.character(A)))
  # cat("\n")
  
  step<-abs(A[6]-A[1])/10
  
  if(step == 0)
  {
    
    step<-1
  }
  
  breaks.Rank<-unique(seq(from= A[1], to=A[6]+step,by=step))
  labels.Rank<-as.character(round(breaks.Rank,3))
  
  levels_ensembl_gene_id<-unique(REP$ensembl_gene_id)
  
  if(DEBUG == 1)
  {
    cat("labels.Rank:\t")
    cat(sprintf(as.character(labels.Rank)))
    cat("\n")
    
    cat("levels_ensembl_gene_id:\t")
    cat(str(levels_ensembl_gene_id))
    cat("\n")
    
    # quit(status = 1)
  }
  
  
  
  
  graph_log2_counts<-ggplot(data=REP,aes(x=Sample_label, y=log2_counts, fill=sample)) +
    geom_violin()+
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", width = 0.5)+
    theme_bw()+
    scale_x_discrete(name=NULL, drop=F)+
    scale_y_continuous(name="Gene EXP normalised log2(counts+1)",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
    scale_fill_brewer(palette = "Set3", drop=F)+
    ggtitle(paste("Normalization_graph:",length(levels_ensembl_gene_id),"total genes filtered", sep=' '))+
    ggeasy::easy_center_title()
  
  if(DEBUG == 1)
  {
    cat("Graph_Part_I:\t")
    
  }
  
  graph_log2_counts<-graph_log2_counts+
    facet_grid(REP$NOR ~ ., scales='free_x', space='free_x') +
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
    scale_x_discrete(name=NULL, drop=T)

  if(DEBUG == 1)
  {
    cat("Graph_Part_II:\t")

  }
  
  setwd(path_graphs)
  
  svgname<-paste("Normalization_plot_",type,".svg",sep='')
  makesvg = TRUE
  
  if (makesvg == TRUE)
  {
    ggsave(svgname, plot= graph_log2_counts,
           device="svg",
           height=10, width=12)
  }
  
  if(DEBUG == 1)
  {
    cat("Graph_Part_END:\t")
    
  }
  
  #### SAVE ----
  
  
  setwd(out)
  
  write.table(PRE_NOR, file="PRE_NOR.tsv", sep="\t", quote=F, row.names = T)
  write.table(POST_NOR, file="POST_NOR.tsv", sep="\t", quote=F, row.names = T)

  write.table(Nor_factors_subset_matrix_inverse, file="Nor_factors_subset_matrix_inverse.tsv", sep="\t", quote=F, row.names = F)
  
  write.table(Nor_factors_subset_matrix_inverse_scaled, file="Nor_factors_subset_matrix_inverse_scaled.tsv", sep="\t", quote=F, row.names = F)
  
  
}


time_course_DE = function(option_list)
{
  suppressMessages(library("edgeR", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  
  
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
  
  #### READ and transform selection_samples ----
  
  selection_samples = unlist(strsplit(opt$selection_samples, split=','))
  
  cat("selection_samples_\n")
  cat(sprintf(as.character(selection_samples)))
  cat("\n")
  
  #### READ and transform selection_treatment ----
  
  selection_treatment = unlist(strsplit(opt$selection_treatment, split=','))
  
  cat("selection_treatment_\n")
  cat(sprintf(as.character(selection_treatment)))
  cat("\n")
  
  #### READ and transform selection_time_point ----
  
  selection_time_point = unlist(strsplit(opt$selection_time_point, split=','))
  
  cat("selection_time_point_\n")
  cat(sprintf(as.character(selection_time_point)))
  cat("\n")
  
  #### Read manifest_FINAL_with_results file ----
  
  #"manifest_FINAL_with_results.rds"
  
  manifest_FINAL_with_results<-readRDS(file=opt$manifest_FINAL_with_results)
  
  manifest_FINAL_with_results$Genotype<-'NA'
  
  manifest_FINAL_with_results$Genotype[which(manifest_FINAL_with_results$sample%in%c("WT_A","WT_B","WT_C"))]<-'wt'
  manifest_FINAL_with_results$Genotype[which(manifest_FINAL_with_results$sample%in%c("clone_13","clone_27","clone_29"))]<-'homALT'
  manifest_FINAL_with_results$Genotype[which(manifest_FINAL_with_results$sample%in%c("del_233","del_235","del_287"))]<-'Del80'
  
  manifest_FINAL_with_results$Genotype<-factor(manifest_FINAL_with_results$Genotype,
                                               levels=c('wt','homALT','Del80'),
                                               ordered=T)
  
  
  cat("manifest_FINAL_with_results_0\n")
  cat(str(manifest_FINAL_with_results))
  cat("\n")
  cat(sprintf(as.character(names(summary(manifest_FINAL_with_results$sample)))))
  cat("\n")
  cat(sprintf(as.character(summary(manifest_FINAL_with_results$sample))))
  cat("\n")
  cat(sprintf(as.character(names(summary(manifest_FINAL_with_results$treatment)))))
  cat("\n")
  cat(sprintf(as.character(summary(manifest_FINAL_with_results$treatment))))
  cat("\n")
  cat(sprintf(as.character(names(summary(manifest_FINAL_with_results$time_point)))))
  cat("\n")
  cat(sprintf(as.character(summary(manifest_FINAL_with_results$time_point))))
  cat("\n")
  
  manifest_after_sel<-droplevels(manifest_FINAL_with_results[which(manifest_FINAL_with_results$sample%in%selection_samples &
                                                                     manifest_FINAL_with_results$treatment%in%selection_treatment &
                                                                     manifest_FINAL_with_results$time_point%in%selection_time_point),])
  
  
  
  selection_sample_labels<-unique(manifest_after_sel$Sample_label)
  
  cat("manifest_after_sel_0\n")
  cat(str(manifest_after_sel))
  cat("\n")
  cat(sprintf(as.character(names(summary(manifest_after_sel$sample)))))
  cat("\n")
  cat(sprintf(as.character(summary(manifest_after_sel$sample))))
  cat("\n")
  cat(sprintf(as.character(names(summary(manifest_after_sel$treatment)))))
  cat("\n")
  cat(sprintf(as.character(summary(manifest_after_sel$treatment))))
  cat("\n")
  cat(sprintf(as.character(names(summary(manifest_after_sel$time_point)))))
  cat("\n")
  cat(sprintf(as.character(summary(manifest_after_sel$time_point))))
  cat("\n")
  cat(sprintf(as.character(names(summary(manifest_after_sel$Genotype)))))
  cat("\n")
  cat(sprintf(as.character(summary(manifest_after_sel$Genotype))))
  cat("\n")
  cat("\n")
  cat("\n")
  cat(sprintf(as.character(selection_sample_labels)))
  cat("\n")
  
  design_df<-unique(manifest_after_sel[,c(which(colnames(manifest_after_sel) == 'sample'),
                                          which(colnames(manifest_after_sel) == 'Genotype'),
                                          which(colnames(manifest_after_sel) == 'treatment'),
                                          which(colnames(manifest_after_sel) == 'time_point'),
                                          which(colnames(manifest_after_sel) == 'Sample_label'))])
  
  cat("design_df_0\n")
  cat(str(design_df))
  cat("\n")
  
  
  basal_df<-as.data.frame(rbind(c("WT_A","wt","5nM_PMA","Basal","WT_A__untreated__16_hrs__S76"),
                                c("WT_B","wt","5nM_PMA","Basal","WT_B__untreated__16_hrs__S78"),
                                c("WT_C","wt","5nM_PMA","Basal","WT_C__untreated__16_hrs__S80"),
                                c("clone_13","homALT","5nM_PMA","Basal","clone_13__untreated__16_hrs__S82"),
                                c("clone_27","homALT","5nM_PMA","Basal","clone_27__untreated__16_hrs__S84"),
                                c("clone_29","homALT","5nM_PMA","Basal","clone_29__untreated__16_hrs__S86"),
                                c("del_233","Del80","5nM_PMA","Basal","del_233__untreated__16_hrs__S88"),
                                c("del_235","Del80","5nM_PMA","Basal","del_235__untreated__16_hrs__S90"),
                                c("del_287","Del80","5nM_PMA","Basal","del_287__untreated__16_hrs__S92")
  ), stringsAsFactors = F)
  
  colnames(basal_df)<-c("sample","Genotype","treatment","time_point","Sample_label")
  
  
  
  cat("basal_df_0\n")
  cat(str(basal_df))
  cat("\n")
  
  design_df<-rbind(design_df,basal_df)
  
  cat("design_df_1\n")
  cat(str(design_df))
  cat("\n")
  
  design_df$time_point<-factor(design_df$time_point,
                               levels=c('Basal','16_hrs','24_hrs','48_hrs','72_hrs'),
                               ordered=T)
  
  design_df<-design_df[order(design_df$Genotype,design_df$sample,design_df$time_point),]
  
  
  cat("design_df_2\n")
  cat(str(design_df))
  cat("\n")
 
  #### Design matrix  ----
  

  Group<-factor(paste(design_df$Genotype,design_df$time_point, sep='.'))
  Group<-relevel(Group, ref="wt.Basal") ### edgeR doesn't seem to work well with ordered factors
  
  cat("Group_0\n")
  cat(sprintf(as.character(Group)))
  cat("\n")
  
  design_df$Group<-Group
  
  cat("design_df_3\n")
  cat(str(design_df))
  cat("\n")
  
  design_df<-design_df[order(design_df$Group),]
  
  
  Sample_label_vector<-as.character(design_df$Sample_label)
  
  cat("Sample_label_vector_0\n")
  cat(sprintf(as.character(Sample_label_vector)))
  cat("\n")
  
  
  setwd(Master_path_analysis)
  
  saveRDS(design_df, file="design_df.RDS")
  saveRDS(Sample_label_vector, file="Group.RDS")
  
  #### Read the filtered and normalised file ----
  
  setwd(out)
  
  POST_NOR<-as.matrix(read.table(file="POST_NOR.tsv", sep="\t"))
  
  cat("POST_NOR_0\n")
  cat(str(POST_NOR))
  cat("\n")
  
  POST_NOR_subset<-POST_NOR[,which(colnames(POST_NOR)%in%Sample_label_vector)]
  
  cat("POST_NOR_subset_0\n")
  cat(str(POST_NOR_subset))
  cat("\n")
  
  ### REORDER MATRIX TO MATCH GROUP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  POST_NOR_subset_df<-data.frame(POST_NOR_subset)
  POST_NOR_subset_df$ensembl_gene_id<-row.names(POST_NOR_subset)
  
  cat("POST_NOR_subset_1\n")
  cat(str(POST_NOR_subset))
  cat("\n")
  
  Sample_label_vector_2<-c("ensembl_gene_id",Sample_label_vector)
  
  cat("Sample_label_vector_2_0\n")
  cat(sprintf(as.character((Sample_label_vector_2))))
  cat("\n")
  
  POST_NOR_subset_df.dt<-setDT(POST_NOR_subset_df)
  
  cat("POST_NOR_subset_df.dt_0\n")
  cat(str(POST_NOR_subset_df.dt))
  cat("\n")
  
  POST_NOR_subset_df_reordered<-data.frame(setcolorder(POST_NOR_subset_df.dt, Sample_label_vector_2), stringsAsFactors = F)
  
  cat("POST_NOR_subset_df_reordered_0\n")
  cat(str(POST_NOR_subset_df_reordered))
  cat("\n")
  
  POST_NOR_subset_FINAL<-as.matrix(POST_NOR_subset_df_reordered[,-1])
  row.names(POST_NOR_subset_FINAL)<-POST_NOR_subset_df_reordered$ensembl_gene_id
  
  cat("POST_NOR_subset_FINAL_0\n")
  cat(str(POST_NOR_subset_FINAL))
  cat("\n")
  cat(sprintf(as.character(colnames(POST_NOR_subset_FINAL))))
  cat("\n")
  
  FLAG_order_columns<-sum(colnames(POST_NOR_subset_FINAL) == Sample_label_vector)
  
  cat("FLAG_order_columns_0\n")
  cat(str(FLAG_order_columns))
  cat("\n")
  
  if(FLAG_order_columns == dim(POST_NOR_subset_FINAL)[2])
  {
    #### Read into edgeR ----
    
    edgeR_primary <- DGEList(counts=POST_NOR_subset_FINAL, group=as.character(design_df$Group)) ########################## Make sure the order of the columns in the matrix matches the order of levels in Group !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    cat("edgeR_primary_0\n")
    cat(str(edgeR_primary))
    cat("\n")
    cat(sprintf(as.character(colnames(POST_NOR_subset_FINAL))))
    cat("\n")
    cat(sprintf(as.character(design_df$Group)))
    cat("\n")
    
   
    
    #### Design matrix ----
    
    
    design <- model.matrix(~ 0 + design_df$Group)
    
    colnames(design) <- levels(design_df$Group)
    
    cat("design_0\n")
    cat(str(design))
    cat("\n")
    cat(sprintf(as.character(colnames(design))))
    cat("\n")
    
    
    
    ### Limma test 
    
    Limma_Full_Rank_test<-is.fullrank(design)
    
    cat("Limma_Full_Rank_test_0\n")
    cat(sprintf(as.character(Limma_Full_Rank_test)))
    cat("\n")
    cat(sprintf(as.character(sum(Limma_Full_Rank_test))))
    cat("\n")
    
    #### Estimate dispersion ----
    
    edgeR_primary <- estimateDisp(edgeR_primary, design)
    
    cat("edgeR_primary_estimate dispersion\n")
    cat(str(edgeR_primary))
    cat("\n")
    
    #### fit a linear model ----
    
    edgeR_fit <- glmQLFit(edgeR_primary, design)
    
    cat("edgeR_fit_0\n")
    cat(str(edgeR_fit))
    cat("\n")
    
    
    
    setwd(Master_path_analysis)
    
    logCPM.obs <- cpm(edgeR_primary, log=TRUE, prior.count=edgeR_fit$prior.count)
    
    cat("logCPM.obs_0\n")
    cat(str(logCPM.obs))
    cat("\n")
    
    logCPM.fit <- cpm(edgeR_fit, log=TRUE)
    
    names(logCPM.fit)<-names(logCPM.obs)
    
    cat("logCPM.fit_0\n")
    cat(str(logCPM.fit))
    cat("\n")
    
    ## Save counts matrixes
    
    write.table(file='CPM_obs.tsv', logCPM.obs, quote=F, sep="\t")
    write.table(file='CPM_fit.tsv', logCPM.fit, quote=F, sep="\t")
    write.table(file='POST_NOR_subset.tsv', POST_NOR_subset, quote=F, sep="\t")
    
    
    
    
    ### contrast matrix ----
    
    my.contrasts <- makeContrasts(
      homALT.BasalVSwt.Basal = homALT.Basal-wt.Basal,
      homALT.16_hrsVSwt.16_hrs = homALT.16_hrs-wt.16_hrs,
      homALT.24_hrsVSwt.24_hrs = homALT.24_hrs-wt.24_hrs,
      homALT.48_hrsVSwt.48_hrs = homALT.48_hrs-wt.48_hrs,
      homALT.72_hrsVSwt.72_hrs = homALT.72_hrs-wt.72_hrs,
      Del80.BasalVSwt.Basal = Del80.Basal-wt.Basal,
      Del80.16_hrsVSwt.16_hrs = Del80.16_hrs-wt.16_hrs,
      Del80.24_hrsVSwt.24_hrs = Del80.24_hrs-wt.24_hrs,
      Del80.48_hrsVSwt.48_hrs = Del80.48_hrs-wt.48_hrs,
      Del80.72_hrsVSwt.72_hrs = Del80.72_hrs-wt.72_hrs,
      homALT.BasalVSDel80.Basal = homALT.Basal-Del80.Basal,
      homALT.16_hrsVSDel80.16_hrs = homALT.16_hrs-Del80.16_hrs,
      homALT.24_hrsVSDel80.24_hrs = homALT.24_hrs-Del80.24_hrs,
      homALT.48_hrsVSDel80.48_hrs = homALT.48_hrs-Del80.48_hrs,
      homALT.72_hrsVSDel80.72_hrs = homALT.72_hrs-Del80.72_hrs,
      wt.BasalVSwt.16_hrs = wt.16_hrs-wt.Basal,
      wt.16_hrsVSwt.24_hrs = wt.24_hrs-wt.16_hrs,
      wt.24_hrsVSwt.48_hrs = wt.48_hrs-wt.24_hrs,
      wt.48_hrsVSwt.72_hrs = wt.72_hrs-wt.48_hrs,
      levels=design)
    
    cat("my.contrasts_0\n")
    cat(str(my.contrasts))
    cat("\n")
    
    
    #### LOOP to calculate contrast -----
    
    
    Master_df<-as.data.frame(rbind(cbind(rep('homALT_vs_wt',5),c("homALT.BasalVSwt.Basal","homALT.16_hrsVSwt.16_hrs","homALT.24_hrsVSwt.24_hrs","homALT.48_hrsVSwt.48_hrs","homALT.72_hrsVSwt.72_hrs")),
                     cbind(rep('Del80_vs_wt',5),c("Del80.BasalVSwt.Basal","Del80.16_hrsVSwt.16_hrs","Del80.24_hrsVSwt.24_hrs","Del80.48_hrsVSwt.48_hrs","Del80.72_hrsVSwt.72_hrs")),
                     cbind(rep('homALT_vs_Del80',5),c("homALT.BasalVSDel80.Basal","homALT.16_hrsVSDel80.16_hrs","homALT.24_hrsVSDel80.24_hrs","homALT.48_hrsVSDel80.48_hrs","homALT.72_hrsVSDel80.72_hrs")),
                     cbind(rep('Time_course_wt',3),c("wt.BasalVSwt.16_hrs","wt.24_hrsVSwt.48_hrs","wt.48_hrsVSwt.72_hrs"))), stringsAsFactors = F)
    
    colnames(Master_df)<-c('analysis','contrasts')
    
    Master_df$analysis<-factor(Master_df$analysis,
                               levels=c("Time_course_wt","homALT_vs_Del80","Del80_vs_wt","homALT_vs_wt"),
                               ordered = T)
    
    cat("Master_df_0\n")
    cat(str(Master_df))
    cat("\n")
    
    array_levels_analysis<-levels(Master_df$analysis)
    
    cat("array_levels_analysis_0\n")
    cat(str(array_levels_analysis))
    cat("\n")
    
    for(i in 1:length(array_levels_analysis))
    {
      array_levels_analysis_sel<-array_levels_analysis[i]
      
      cat("--------------------------------->\t")
      cat(sprintf(as.character(array_levels_analysis_sel)))
      cat("\n")
      
      sub_path_analysis<-paste(Master_path_analysis,array_levels_analysis_sel,'/',sep='')
      
      if (file.exists(sub_path_analysis)){
        
        
      }else{
        
        dir.create(file.path(sub_path_analysis))
        
      }#sub_path_analysis
      
      Master_df_sel<-Master_df[which(Master_df$analysis == array_levels_analysis_sel),]
      
      cat("Master_df_sel_0\n")
      cat(str(Master_df_sel))
      cat("\n")
      
      qlf <- glmQLFTest(edgeR_fit, contrast=my.contrasts[,Master_df_sel$contrasts])
      
      results <- as.data.frame(topTags(qlf,n= Inf, adjust.method = "BH", sort.by = "PValue", p.value = 1))
      
      results$ensembl_gene_id<-row.names(results)
      row.names(results)<-NULL
      
      # cat("results_0\n")
      # cat(str(results))
      # cat("\n")
      
      multiVals <- function(x) paste(x,collapse=";")
      results$Symbol <- mapIds(org.Hs.eg.db, keys=results$ensembl_gene_id, keytype="ENSEMBL",
                               column="SYMBOL", multiVals=multiVals)
      
      # cat("results_1\n")
      # cat(str(results))
      # cat("\n")
      
      setwd(sub_path_analysis)
      
      write.table(file=paste('results_',array_levels_analysis_sel,'.tsv', sep=''), results, quote=F,row.names = F, sep="\t")
      
      
      # #################################################
      # quit(status = 1)
      
    }# i in 1:length(array_levels_analysis)
  }#FLAG_order_columns == dim(POST_NOR_subset_FINAL)[2]
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
  
  filter_and_nor(opt)
  time_course_DE(opt)


  
}


###########################################################################

system.time( main() )