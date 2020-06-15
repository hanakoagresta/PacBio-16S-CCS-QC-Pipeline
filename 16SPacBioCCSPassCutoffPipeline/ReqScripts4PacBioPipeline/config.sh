#/bin/bash

############
# INSTRUCTIONS FOR EACH VARIABLE
#
# query_fasta is the path to the .fasta file with positive control samples
# outfile is a user-selected name for the files generated in intermediate steps. It should be a .txt file.
# reference_fasta is the path to the .fasta file with all of the 16S full length sequences for the positive control.
# python3env is the path to the bin where the python 3.7.3 environment is stored.
# python2env is the path to the bin where the python 2.7.15 environment is stored.
# phrap4cross_match is the path to the phrap module (which contains cross_match). Directions to download are found here: http://www.phrap.org/consed/distributions/README.29.0.txt
# allseqs_fasta is the path to the .fasta file with all of the sample sequences.
# finaloutputfasta is the user-selected name for the final output .fasta file with only sequences that fit the upper and lower CCS pass thresholds.
############

query_fasta="/projects/weinstocklab/comp/pipelines/pacbio/Hanako/InputFastas/subset_A01.fasta"
outfile="testingtesting3.txt"
reference_fasta="/projects/weinstocklab/comp/pipelines/pacbio/Hanako/16SPacBio_Pipeline/mock36reference_FINAL.fasta"
python3env="/projects/weinstocklab/students/Hanako/devel/conda/bin"
python2env="/projects/weinstocklab/comp/local/anaconda2/bin"
phrap4cross_match="/projects/weinstocklab/comp/local/phrap"
allseqs_fasta="/projects/weinstocklab/comp/pipelines/pacbio/Hanako/InputFastas/subset_A01.fasta"
finaloutputfasta="FINALOUTPUT.fasta"

###########
# THE FOLLOWING ARE CODED AS OPTIONAL VARIABLES IN THE STANDALONE model_positive_congrol.py SCRIPT. HOWEVER, WITHIN THIS PIPELINE, IT IS NECESSARY TO KEEP THESE AS DEFAULT SETTINGS OR SET THEM MANUALLY RATHER THAN DELETE THEM FROM THE config.sh FILE
#
# PEC is the percent error cutoff, and indicates the user-selected percent error of the model (based off of the user-selected percentile or mean) over which the corresponding CCS pass number is selected. The default is 1% and can be set to any number. An error will arise if the percent error cutoff is below the asymptote of the exponential decay model.
# P is the percentile, and indicates the percentile of the data that is used to create the model. The default is set to the 75th percentile and can be set to any integer.
# UL is the upper CCS pass limit, and indicates user selected upper limit. The default is set to 75 CCS passes and can be set to any integer.
# DUFM is the data used for the model, and can either be the user-selected percentile determined by P (denoted by "nthpercent") or the mean (denoted by "mean"). The default is "nthpercent" and the only two options for this selection are 1. nthpercent or 2. mean.
###########

PEC="1"
P="75"
UL="75"
DUFM="nthpercent"
