Before running the 16S PacBio pipeline, the following must be configured:
1. Ensure that you have the following scripts and files:
	- pipeline_ref_seq_error_profiles.sh
	- summarize_cross_match_mapping.py
	- config.sh
	- model_positive_congrol.py
2. Ensure that the following scripts are located in a directory labeled ReqScripts4PacBioPipeline
	- summarize_cross_match_mapping.py
        - config.sh
        - model_positive_congrol.py
3  Ensure that you have an empty directory labeled OUTPUT_16SPacbioPipeline. Create one if necessary. This will be where all of the outputs will be stored.
4. Update the config.sh file to include the paths to each of the variables listed
	- query_fasta is the path to the .fasta file with positive control samples
	- outfile is a user-selected name for the files generated in intermediate steps. It should be a .txt file.
	- reference_fasta is the path to the .fasta file with all of the 16S full length sequences for the positive control.
	- python3env is the path to the bin where the python 3.7.3 environment is stored.
	- python2env is the path to the bin where the python 2.7.15 environment is stored.
	- phrap4cross_match is the path to the phrap module (which contains cross_match). Directions to download are found here: http://www.phrap.org/consed/distributions/README.29.0.txt
	- finaloutputfasta is a user-selected name for the final .fasta file that has removed all sequences with a CCS pass number outside of the pipeline calculated range
   
5. Update the config.sh file to select the following OPTIONAL preferences in the pipeline:
	- percent cutoff 
	- percentile
	- upper CCS pass limit
	- data used for model
	Note that if you do not wish to vary from the default, leave the settings in place from config.sh. Do not delete the variables.
6. The pipeline must be run from the same directory as the test_pipeline_ref_seq_error_profiles.sh script
7. The final output of sequences with a new CCS pass cutoff will be located in a new directory called OUTPUT_16SPacbioPipeline


To run the 16S PacBio pipeline use the command:
sbatch test_pipeline_ref_seq_error_profiles.sh
