#!/bin/bash
#
#SBATCH-N 1
#SBATCH-n 8
#SBATCH--mem 250G
#SBATCH-t 0-1:00:00
#SBATCH-q long

##  How to run this pipeline on SUMNER cluster
##: sbatch --export=ALL,var1="mytest.fasta",var2="myoutput.txt" /projects/weinstocklab/comp/pipelines/pacbio/Hanako/pipeline_ref_seq_error_profiles.sh

export PATH=/projects/weinstocklab/comp/local/phrap:$PATH
export PATH=/projects/weinstocklab/comp/local/anaconda2/bin:$PATH




query_fasta=$var1
outfile=$var2
reference_fasta=github_mock36reference_FINAL.fasta

cross_match  -discrep_lists -tags -masklevel 0 ${query_fasta} ${reference_fasta} > ${outfile}

python github_summarize_cross_match_mapping.py ${outfile} $reference_fasta ${outfile}.snpindel.table.txt

awk  '{ print $5"\t"($2 + $3 + $4)"\t"$2"\t"$3"\t"$4}' ${outfile}.snpindel.table.txt | awk '!/SNP/ {print}' |awk 'BEGIN {FS = "[:]+" } ; {print $3}' | awk 'BEGIN{printf "CPN\tTOTAL\tSNP\tINS\tDEL\n"} {print}' > ${outfile}.tsv

python github_hardcode_model_positive_control.py