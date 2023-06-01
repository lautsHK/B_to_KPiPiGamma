#!/bin/bash

# Define the relative path to the output directory
input_dir="./../output"
output_dir="./../output/combined_output"

# Define the input name
charge="1" # Only 0 or 1
sqrt_s_min_str="0p778"
sqrt_s_max_str="1p8"
mode="1" # 0 for constant sqrt_s, 1 for 0.778 < sqrt_s < sqrt_s_max.
str1270="1p36603"
str1400="m0p366025"
str1410="0p"
str1430="0p8"
str1680="0p"
lambda_str="m1p0" # lambda should be between -1.0 and 1.0 (1.0 for fully right hand polarised, and -1.0 for left hand.)


# Get the current date in mmdd format
date_str=$(date +"%m%d")

# Loop over the strings "CosTheta", "Phi", and other variables  
if [ $mode -eq 0 ]; then
    jobname="Charge_${charge}_Res_${str1270}_${str1400}_${str1410}_${str1430}_${str1680}_Lambda_${lambda_str}_ConstSqS_${sqrt_s_max_str}"
else
    jobname="Charge_${charge}_Res_${str1270}_${str1400}_${str1410}_${str1430}_${str1680}_Lambda_${lambda_str}_VarSqS_${sqrt_s_min_str}_to_${sqrt_s_max_str}"
fi

echo $jobname

for observables in CosTheta phi s s13 s23 s13s23 ME2 SqrtS J2 nJXJc P1 P2 P3; do

    # Define the input files and the output file
    input_files="${input_dir}/${jobname}_[0-9]*_${observables}.dat"
    output_file="${output_dir}/${jobname}_${observables}_${date_str}.dat"

    # Remove the output file if it already exists
    if [ -f "$output_file" ]; then
        rm "$output_file"
    fi

    # Loop over the input files and append their contents to the output file
    for file in $input_files; do
        cat "$file" >> "$output_file"
        #echo $file
    done

done

