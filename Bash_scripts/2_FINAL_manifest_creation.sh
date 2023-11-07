#!/bin/bash
module load R/4.1.0
Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
Experiment=$2


Log_files=$(echo "$MASTER_ROUTE""Log_files/")

#rm -rf $Log_files
#mkdir -p $Log_files


#### Merge_sample_sheet_and_manifest

Rscript_Merge_sample_sheet_and_manifest=$(echo "$Rscripts_path""1_Merge_sample_sheet_and_manifest.R")

type=$(echo "Merge_sample_sheet_and_manifest")


outfile_Merge_sample_sheet_and_manifest=$(echo "$Log_files""outfile_""$type""_""$Experiment"".log")
touch $outfile_Merge_sample_sheet_and_manifest
echo -n "" > $outfile_Merge_sample_sheet_and_manifest
name_Merge_sample_sheet_and_manifest=$(echo "$type""_""$Experiment""_job")

sample_sheet=$(echo "$MASTER_ROUTE""sample_sheet_corrected.tsv")
manifest=$(echo "$MASTER_ROUTE""manifest.tsv")
path_samples=$(echo "/processing_data/post_processing/genomics/soranzo/RITM0022181_2/fastq/")
path_scratch=$(echo "/scratch/manuel.tardaguila/""$Experiment""/")


myjobid_Merge_sample_sheet_and_manifest=$(sbatch --job-name=$name_Merge_sample_sheet_and_manifest --output=$outfile_Merge_sample_sheet_and_manifest --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=512M --parsable --wrap="Rscript $Rscript_Merge_sample_sheet_and_manifest --sample_sheet $sample_sheet --manifest $manifest --path_scratch $path_scratch --type $type --path_samples $path_samples")



#### next_flow rna-seq pipeline

type=$(echo "nf_run")


outfile_nf_run=$(echo "$Log_files""outfile_""$type""_""$Experiment"".log")
touch $outfile_nf_run
echo -n "" > $outfile_nf_run
name_nf_run=$(echo "$type""_""$Experiment""_job")

#setup
module load nextflow/22.10.1;
module load singularity/3.6.3;
export XDG_RUNTIME_DIR=$TMPDIR; # tmp folder to the system temporary folder
export NXF_OPTS='-Xms1g -Xmx4g'; # limiting Java VM memory

config_file=$(echo "/home/manuel.tardaguila/Scripts/Luigi/common.config")
work_dir=$(echo "$path_scratch""work""/")
input=$(echo "$path_scratch""nfs_input.csv")
fasta=$(echo "/scratch/manuel.tardaguila/RITM0022181/reference_files/Ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
gtf=$(echo "/scratch/manuel.tardaguila/RITM0022181/reference_files/Ensembl/Homo_sapiens.GRCh38.110.gtf")

myjobid_nf_run=$(sbatch --dependency=afterany:$myjobid_Merge_sample_sheet_and_manifest --job-name=$name_nf_run --output=$outfile_nf_run --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem=4G --parsable --wrap="nextflow run nf-core/rnaseq -r 3.10.1 -c $config_file -resume -work-dir $work_dir -profile manuel --input $input --outdir $path_scratch --fasta $fasta --gtf $gtf --save_reference TRUE --aligner star_rsem")
