#!/bin/bash

# This script runs pyradiomics on the masks and images according to 'write_csv.sh' for the...
# manual and semi-automatic methods, then dilates and erodes the radiologist masks only then extracts the ICCs...
# based on the Radiologist Original Mask, assessing variation across: i)semi-automatic delineation, ii)eroded and iii) dilated manual delineations.

# Inputs: name of parameter file $1 (e.g. parameters.yml), dilation/erosion pixel number $2/$3,
# CSV name with features and ICC scores $4 (e.g. CSV_radio='csvs/bw35_icc_score_list_radio.csv')

param_filename=$1
dilation_amount=$2
erosion_amount=$3
CSV_radio=$4

################################################
######## 1-Run standard pyradiomics ##############
################################################
# export the path containing the pyradiomics module
cd /home/zm297/rds/hpc-work/Radiomics
export PATH="/home/zm297/.local/bin:$PATH"

echo Extracting Radiomics features using settings in YML file $1

# the below function writes a csv file with all the masks and images of the patients contained in the 'rds/Radiomics/Patients/' directory. 
# then, run pyradiomics to produce csv_output csv file

source write_csv.sh 'csv_input_radio.csv' '/home/zm297/rds/hpc-work/Radiomics/Patients/RADIOLOGIST/'

pyradiomics csv_input_radio.csv -o csv_output_radio.csv -f csv --param $1


##########################################
########### 2-do dilation. #################
##########################################

python dilation.py $dilation_amount '/home/zm297/rds/hpc-work/Radiomics/Patients/RADIOLOGIST/'

source write_csv_dilated.sh '/home/zm297/rds/hpc-work/Radiomics/csv_input_dilated_radio.csv' '/home/zm297/rds/hpc-work/Radiomics/Patients/RADIOLOGIST/'

pyradiomics csv_input_dilated_radio.csv -o csv_output_dilated_radio.csv -f csv --param $1

##########################################
########### 3-erosion... #############
##########################################

python erosion.py $erosion_amount '/home/zm297/rds/hpc-work/Radiomics/Patients/RADIOLOGIST/'

source write_csv_eroded.sh '/home/zm297/rds/hpc-work/Radiomics/csv_input_eroded_radio.csv' '/home/zm297/rds/hpc-work/Radiomics/Patients/RADIOLOGIST/'

pyradiomics csv_input_eroded_radio.csv -o csv_output_eroded_radio.csv -f csv --param $1

#######################################################
###### 4-get the ICC scores ######
#######################################################


echo Saving ICC scores to: $CSV_radio
python icc_scores_replica.py 'csv_output_radio.csv' 'csv_output_dilated_radio.csv' 'csv_output_eroded_radio.csv' $CSV_radio
echo Done ICC scores.



# further analysis...
#source binwidth_count.sh
#
