#!/bin/bash

eval "$(conda shell.bash hook)"



Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
Experiment=$2
analysis=$3

output_dir=$(echo "$MASTER_ROUTE""$Experiment""/")

#rm -rf $output_dir
#mkdir -p $output_dir

Master_path_analysis=$(echo "$output_dir""$analysis""/")

#rm -rf $Master_path_analysis
#mkdir -p $Master_path_analysis



Log_files=$(echo "$Master_path_analysis""/""Log_files/")

#rm -rf $Log_files
#mkdir -p $Log_files



#### collect_results_and_normalize

module load R/4.1.0

Rscript_collect_results_and_normalize=$(echo "$Rscripts_path""2_Collect_normalize_edgeR_DeSeq2.R")

type=$(echo "collect_results_and_normalize")


outfile_collect_results_and_normalize=$(echo "$Log_files""outfile_""$type""_""$Experiment"".log")
touch $outfile_collect_results_and_normalize
echo -n "" > $outfile_collect_results_and_normalize
name_collect_results_and_normalize=$(echo "$type""_""$Experiment""_job")

path_to_results=$(echo "/scratch/manuel.tardaguila/""$Experiment""/star_rsem/")
manifest=$(echo "/home/manuel.tardaguila/RNA_seq_time_course/""manifest_FINAL.tsv")

myjobid_collect_results_and_normalize=$(sbatch --job-name=$name_collect_results_and_normalize --output=$outfile_collect_results_and_normalize --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_collect_results_and_normalize --path_to_results $path_to_results --manifest $manifest --type $type --out $output_dir")


#### edgeR part

module load R/4.1.0

scratch_path=$(echo "/scratch/manuel.tardaguila/""$Experiment""/")


type=$(echo "edgeR_suite")
outfile_edgeR_suite=$(echo "$Log_files""outfile_""$type""_""$analysis""_""$Experiment"".log")
touch $outfile_edgeR_suite
echo -n "" > $outfile_edgeR_suite
name_edgeR_suite=$(echo "$type""_""$Experiment""_job")


Rscript_edgeR_suite=$(echo "$Rscripts_path""3_edgeR_filter_nor_time_course_v3.R")

gene_expression_results=$(echo "$output_dir""Gene_expression_results.tsv")
manifest_FINAL_with_results=$(echo "$output_dir""manifest_FINAL_with_results.rds")
#selection_samples=$(echo "WT_A,WT_B,WT_C,clone_13,clone_27,clone_29,del_233,del_235,del_287,hESC_MK")
selection_samples=$(echo "WT_A,WT_B,WT_C,clone_13,clone_27,clone_29,del_233,del_235,del_287")
selection_treatment=$(echo "5nM_PMA")
selection_time_point=$(echo "Basal,16_hrs,24_hrs,48_hrs,72_hrs")


myjobid_edgeR_suite=$(sbatch --dependency=afterany:$myjobid_collect_results_and_normalize --job-name=$name_edgeR_suite --output=$outfile_edgeR_suite --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4096 --parsable --wrap="Rscript $Rscript_edgeR_suite --gene_expression_results $gene_expression_results --manifest_FINAL_with_results $manifest_FINAL_with_results --Master_path_analysis $Master_path_analysis --selection_samples $selection_samples --selection_treatment $selection_treatment --selection_time_point $selection_time_point --type $type --out $output_dir")
myjobid_seff_edgeR_suite=$(sbatch --dependency=afterany:$myjobid_edgeR_suite --open-mode=append --output=$outfile_edgeR_suite --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_edgeR_suite >> $outfile_edgeR_suite")

# myjobid_edgeR_suite=$(sbatch --job-name=$name_edgeR_suite --output=$outfile_edgeR_suite --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4096 --parsable --wrap="Rscript $Rscript_edgeR_suite --gene_expression_results $gene_expression_results --manifest_FINAL_with_results $manifest_FINAL_with_results --Master_path_analysis $Master_path_analysis --selection_samples $selection_samples --selection_treatment $selection_treatment --selection_time_point $selection_time_point --type $type --out $output_dir")
# myjobid_seff_edgeR_suite=$(sbatch --dependency=afterany:$myjobid_edgeR_suite --open-mode=append --output=$outfile_edgeR_suite --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_edgeR_suite >> $outfile_edgeR_suite")

#### edgeR graphs

module load R/4.1.0

scratch_path=$(echo "/scratch/manuel.tardaguila/""$Experiment""/")


type=$(echo "edgeR_graphs")
outfile_edgeR_graphs=$(echo "$Log_files""outfile_""$type""_""$analysis""_""$Experiment"".log")
touch $outfile_edgeR_graphs
echo -n "" > $outfile_edgeR_graphs
name_edgeR_graphs=$(echo "$type""_""$Experiment""_job")


Rscript_edgeR_graphs=$(echo "$Rscripts_path""3_5_edgeR_graphs.R")

gene_expression_results=$(echo "$output_dir""Gene_expression_results.tsv")
manifest_FINAL_with_results=$(echo "$output_dir""manifest_FINAL_with_results.rds")
#selection_samples=$(echo "WT_A,WT_B,WT_C,clone_13,clone_27,clone_29,del_233,del_235,del_287,hESC_MK")
selection_samples=$(echo "WT_A,WT_B,WT_C,clone_13,clone_27,clone_29,del_233,del_235,del_287")
selection_treatment=$(echo "5nM_PMA")
selection_time_point=$(echo "Basal,16_hrs,24_hrs,48_hrs,72_hrs")

myjobid_edgeR_graphs=$(sbatch --dependency=afterany:$myjobid_edgeR_suite --job-name=$name_edgeR_graphs --output=$outfile_edgeR_graphs --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024 --parsable --wrap="Rscript $Rscript_edgeR_graphs --gene_expression_results $gene_expression_results --manifest_FINAL_with_results $manifest_FINAL_with_results --Master_path_analysis $Master_path_analysis --selection_samples $selection_samples --selection_treatment $selection_treatment --selection_time_point $selection_time_point --type $type --out $output_dir")
 myjobid_seff_edgeR_graphs=$(sbatch --dependency=afterany:$myjobid_edgeR_graphs --open-mode=append --output=$outfile_edgeR_graphs --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_edgeR_graphs >> $outfile_edgeR_graphs")

#myjobid_edgeR_graphs=$(sbatch --job-name=$name_edgeR_graphs --output=$outfile_edgeR_graphs --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024 --parsable --wrap="Rscript $Rscript_edgeR_graphs --gene_expression_results $gene_expression_results --manifest_FINAL_with_results $manifest_FINAL_with_results --Master_path_analysis $Master_path_analysis --selection_samples $selection_samples --selection_treatment $selection_treatment --selection_time_point $selection_time_point --type $type --out $output_dir")
#myjobid_seff_edgeR_graphs=$(sbatch --dependency=afterany:$myjobid_edgeR_graphs --open-mode=append --output=$outfile_edgeR_graphs --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_edgeR_graphs >> $outfile_edgeR_graphs")


#### Separate the scrip to launch graphs from 3_edgeR_filter_nor_time_course_v3.R so I can run in parallel the plots downstream

declare -a arr

analysis_array=$(echo "homALT_vs_wt,Del80_vs_wt,homALT_vs_Del80,Time_course_wt")

a=($(echo "$analysis_array" | tr "," '\n'))

declare -a arr

for i  in "${a[@]}"
    do
        analysis_array_sel=${i}
        echo "$analysis_array_sel"

        analysis_route=$(echo "$Master_path_analysis""/""$analysis_array_sel""/")
	echo "$analysis_route"

	### MySigDB_GSEA
	
	type=$(echo "$analysis_array_sel""_""MySigDB_GSEA")
	outfile_MySigDB_GSEA=$(echo "$Log_files""outfile_""$type""_""$analysis""_""$Experiment"".log")
	touch $outfile_MySigDB_GSEA
	echo -n "" > $outfile_MySigDB_GSEA
	name_MySigDB_GSEA=$(echo "$type""_""$Experiment""_job")

	
	Rscript_MySigDB_GSEA=$(echo "$Rscripts_path""9_MySigDb_analysis.R")

	edgeR_results=$(echo "$analysis_route""results_""$analysis_array_sel"".tsv")
	path_to_GMT=$(echo "/home/manuel.tardaguila/GMT_files/msigdb_v2023.1.Hs_files_to_download_locally/msigdb_v2023.1.Hs_GMTs/")
	search_terms=$(echo "PLATELET,ERYTHROCYTE,CUX1,MEGAKARYOCYTE,GATA1,GATA2,TET2,RUNX1,RUNX2,PDGFA,HEPATOCYTE,NEURON,LIPID,SPHINGOSINE")   # ADD HSC TERMS
	background_genes=$(echo "/home/manuel.tardaguila/GMT_files/hsapiens.HPA.name.gmt")
	
	myjobid_MySigDB_GSEA=$(sbatch --dependency=afterany:$myjobid_edgeR_suite --job-name=$name_MySigDB_GSEA --output=$outfile_MySigDB_GSEA --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4096 --parsable --wrap="Rscript $Rscript_MySigDB_GSEA --edgeR_results $edgeR_results --path_to_GMT $path_to_GMT --search_terms $search_terms --background_genes $background_genes  --type $type --out $analysis_route")
	
#	myjobid_MySigDB_GSEA=$(sbatch --job-name=$name_MySigDB_GSEA --output=$outfile_MySigDB_GSEA --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4096 --parsable --wrap="Rscript $Rscript_MySigDB_GSEA --edgeR_results $edgeR_results --path_to_GMT $path_to_GMT --search_terms $search_terms --background_genes $background_genes  --type $type --out $analysis_route")


	#### volcano plots

	
	type=$(echo "$analysis_array_sel""_""volcano_plots")
	outfile_volcano_plots=$(echo "$Log_files""outfile_""$type""_""$analysis""_""$Experiment"".log")
	touch $outfile_volcano_plots
	echo -n "" > $outfile_volcano_plots
	name_volcano_plots=$(echo "$type""_""$Experiment""_job")

	
	Rscript_volcano_plots=$(echo "$Rscripts_path""8_volcano_plotter_v2.R")

	edgeR_results=$(echo "$analysis_route""results_""$analysis_array_sel"".tsv")
	MySigDB_results=$(echo "$analysis_route""MySigDb_analysis""/""result""$analysis_array_sel""_MySigDB_GSEA"".tsv")



 	myjobid_volcano_plots=$(sbatch --dependency=afterany:$myjobid_MySigDB_GSEA --job-name=$name_volcano_plots --output=$outfile_volcano_plots --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=1024 --parsable --wrap="Rscript $Rscript_volcano_plots --edgeR_results $edgeR_results --MySigDB_results $MySigDB_results  --type $type --out $analysis_route")
	
#	myjobid_volcano_plots=$(sbatch --job-name=$name_volcano_plots --output=$outfile_volcano_plots --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=1024 --parsable --wrap="Rscript $Rscript_volcano_plots --edgeR_results $edgeR_results --MySigDB_results $MySigDB_results  --type $type --out $analysis_route")

      echo "->>>$myjobid_volcano_plots"
      arr[${#arr[@]}]="$myjobid_volcano_plots"

done


done_string=$(echo "--dependency=afterany:"""""${arr[@]}"""")
echo "$done_string"

dependency_string=$(echo $done_string|sed -r 's/ /:/g')

echo "$dependency_string"

#### collect_DE_and_gene_sets

module load R/4.1.0

Rscript_collect_DE_and_gene_sets=$(echo "$Rscripts_path""10_Collect_DE_and_gene_sets.R")

type=$(echo "DE_collect_and_gene_sets")


outfile_collect_DE_and_gene_sets=$(echo "$Log_files""outfile_""$type""_""$analysis""_""$Experiment"".log")
touch $outfile_collect_DE_and_gene_sets
echo -n "" > $outfile_collect_DE_and_gene_sets
name_collect_DE_and_gene_sets=$(echo "$type""_""$Experiment""_job")

Master_path_analysis=$(echo "$output_dir""$analysis""/")
analysis_array=$(echo "homALT_vs_wt,Del80_vs_wt,homALT_vs_Del80,Time_course_wt")
#analysis_array=$(echo "Time_course_wt,homALT_vs_wt")
maxGS_size=$(echo "500")
minGS_size=$(echo "1")

# myjobid_collect_DE_and_gene_sets=$(sbatch --job-name=$name_collect_DE_and_gene_sets --output=$outfile_collect_DE_and_gene_sets --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=512M --parsable --wrap="Rscript $Rscript_collect_DE_and_gene_sets --analysis_array $analysis_array --minGS_size $minGS_size --maxGS_size $maxGS_size --Master_path_analysis $Master_path_analysis --type $type --out $output_dir")
# myjobid_seff_collect_DE_and_gene_sets=$(sbatch --dependency=afterany:$myjobid_collect_DE_and_gene_sets --open-mode=append --output=$outfile_collect_DE_and_gene_sets --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_collect_DE_and_gene_sets >> $outfile_collect_DE_and_gene_sets")

myjobid_collect_DE_and_gene_sets=$(sbatch $dependency_string --job-name=$name_collect_DE_and_gene_sets --output=$outfile_collect_DE_and_gene_sets --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=512M --parsable --wrap="Rscript $Rscript_collect_DE_and_gene_sets --analysis_array $analysis_array --minGS_size $minGS_size --maxGS_size $maxGS_size --Master_path_analysis $Master_path_analysis --type $type --out $output_dir")
myjobid_seff_collect_DE_and_gene_sets=$(sbatch --dependency=afterany:$myjobid_collect_DE_and_gene_sets --open-mode=append --output=$outfile_collect_DE_and_gene_sets --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_collect_DE_and_gene_sets >> $outfile_collect_DE_and_gene_sets")

#### global dot plots

module load R/4.1.0


Rscript_Global_dotplots=$(echo "$Rscripts_path""11_Dot_plot_DE_plus_GSEA_v2.R")

type=$(echo "Global_dotplots")


outfile_Global_dotplots=$(echo "$Log_files""outfile_""$type""_""$analysis""_""$Experiment"".log")
touch $outfile_Global_dotplots
echo -n "" > $outfile_Global_dotplots
name_Global_dotplots=$(echo "$type""_""$Experiment""_job")

Master_path_analysis=$(echo "$output_dir""$analysis""/")
MySigDB_results_FINAL=$(echo "$Master_path_analysis""GSEA_ActivePathways_results.rds")
edgeR_results_FINAL=$(echo "$Master_path_analysis""DE_edgeR_results.rds")
plot_together=$(echo "homALT_vs_wt,Del80_vs_wt,homALT_vs_Del80")
plot_ctrl=$(echo "Time_course_wt")
ctrl_gene=$(echo "ENSG00000257923,ENSG00000005961")
selected_genes=$(echo "CUX1,IMPDH1,ABCB8,CA1,IKZF2,FYB1,GFI1B,ITGA2B,WAS,IL16,SOCS1,PECAM1,PDGFA,PSMB8,PRF1")
RMV_pathways=$(echo "GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH,DESCARTES_")
KEEP_pathways=$(echo "GOMF_,HP_,REACTOME_,CUX1_")


 myjobid_Global_dotplots=$(sbatch --dependency=afterany:$myjobid_collect_DE_and_gene_set --job-name=$name_Global_dotplots --output=$outfile_Global_dotplots --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_Global_dotplots --MySigDB_results_FINAL $MySigDB_results_FINAL --edgeR_results_FINAL $edgeR_results_FINAL --plot_together $plot_together --plot_ctrl $plot_ctrl --ctrl_gene $ctrl_gene --selected_genes $selected_genes --KEEP_pathways $KEEP_pathways --RMV_pathways $RMV_pathways  --type $type --out $Master_path_analysis")
 myjobid_seff_Global_dotplots=$(sbatch --dependency=afterany:$myjobid_Global_dotplots --open-mode=append --output=$outfile_Global_dotplots --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Global_dotplots >> $outfile_Global_dotplots")
# myjobid_Global_dotplots=$(sbatch --job-name=$name_Global_dotplots --output=$outfile_Global_dotplots --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_Global_dotplots --MySigDB_results_FINAL $MySigDB_results_FINAL --edgeR_results_FINAL $edgeR_results_FINAL --plot_together $plot_together --plot_ctrl $plot_ctrl --ctrl_gene $ctrl_gene --selected_genes $selected_genes --KEEP_pathways $KEEP_pathways --RMV_pathways $RMV_pathways  --type $type --out $Master_path_analysis")
# myjobid_seff_Global_dotplots=$(sbatch --dependency=afterany:$myjobid_Global_dotplots --open-mode=append --output=$outfile_Global_dotplots --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Global_dotplots >> $outfile_Global_dotplots")





