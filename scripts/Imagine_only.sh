#!/bin/bash

# Define the relative path to the output directory
file_dir="/afs/desy.de/user/t/tslau/belle2/Tak/B_to_KPiPiGamma/output/combined_output"

# Define the input name
charge="1" # Only 0 or 1
sqrt_s_min_str="0p778"
sqrt_s_max_str="1p6"
mode="1" # 0 for constant sqrt_s, 1 for 0.778 < sqrt_s < sqrt_s_max.
str1270="1p"
str1400="0p"
str1410="0p"
str1430="0p"
str1680="0p"
lambda_str="m1p0" # lambda should be between -1.0 and 1.0 (1.0 for fully right hand polarised, and -1.0 for left hand.)
date_str="0417"

# Loop over the strings "CosTheta", "Phi", and other variables  
if [ $mode -eq 0 ]; then
    jobname="Charge_${charge}_Res_${str1270}_${str1400}_${str1410}_${str1430}_${str1680}_Lambda_${lambda_str}_ConstSqS_${sqrt_s_max_str}"
else
    jobname="Charge_${charge}_Res_${str1270}_${str1400}_${str1410}_${str1430}_${str1680}_Lambda_${lambda_str}_VarSqS_${sqrt_s_min_str}_to_${sqrt_s_max_str}"
fi

# Specify input and output files
input_file="${file_dir}/${jobname}_nJXJc_${date_str}.dat"
output_file="${file_dir}/${jobname}_nJXJc_value_${date_str}.dat"

# Extract second column and save to output file
sed 's/([^,]*,\([^)]*\)).*/\1/' $input_file > $output_file

