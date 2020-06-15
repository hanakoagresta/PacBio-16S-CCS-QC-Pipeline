#!/bin/bash
#
#SBATCH-N 1
#SBATCH-n 8
#SBATCH--mem 250G
#SBATCH-t 0-1:00:00
#SBATCH-q long

##  How to run this pipeline on SUMNER cluster
##: sbatch /projects/weinstocklab/comp/pipelines/pacbio/Hanako/TRIAL/test_pipeline_ref_seq_error_profiles.sh 

source ReqScripts4PacBioPipeline/config.sh

export PATH=${phrap4cross_match}:$PATH
export PATH=${python2env}:$PATH


cross_match  -discrep_lists -tags -masklevel 0 ${query_fasta} ${reference_fasta} > OUTPUT_16SPacbioPipeline/${outfile}

python ReqScripts4PacBioPipeline/summarize_cross_match_mapping.py OUTPUT_16SPacbioPipeline/${outfile} $reference_fasta OUTPUT_16SPacbioPipeline/${outfile}.snpindel.table.txt

awk  '{ print $5"\t"($2 + $3 + $4)"\t"$2"\t"$3"\t"$4}' OUTPUT_16SPacbioPipeline/${outfile}.snpindel.table.txt | awk '!/SNP/ {print}' |awk 'BEGIN {FS = "[:]+" } ; {print $3}' | awk 'BEGIN{printf "CPN\tTOTAL\tSNP\tINS\tDEL\n"} {print}' > OUTPUT_16SPacbioPipeline/${outfile}.tsv

export PATH=${python3env}:$PATH

python ReqScripts4PacBioPipeline/model_positive_control.py OUTPUT_16SPacbioPipeline/${outfile}.tsv OUTPUT_16SPacbioPipeline/${finaloutputfasta} ${allseqs_fasta} -PEC ${PEC} -P ${P} -UL ${UL} -DUFM ${DUFM}
