import os

# Define the relative path to the output directory
file_dir = "/afs/desy.de/user/t/tslau/belle2/Tak/B_to_KPiPiGamma/output/combined_output"

# Define the input name
charge = "1" # Only 0 or 1
sqrt_s_min_str = "0p778"
sqrt_s_max_str = "1p8"
mode = "1" # 0 for constant sqrt_s, 1 for 0.778 < sqrt_s < sqrt_s_max.
str1270 = "0p"
str1400 = "1p"
str1410 = "0p"
str1430 = "0p"
str1680 = "0p"
lambda_str = "1p0" # lambda should be between -1.0 and 1.0 (1.0 for fully right hand polarised, and -1.0 for left hand.)
date_str = "0417"

# Define the jobname
if mode == "0":
    #jobname = "Charge" + str(charge) + "Res_" + str1270 + "_" + str1400 + ...
    #Below using a f-string, python3 is needed. For python2, please use the above one.
    jobname = f"Charge_{charge}_Res_{str1270}_{str1400}_{str1410}_{str1430}_{str1680}_Lambda_{lambda_str}_ConstSqS_{sqrt_s_max_str}"
else:
    jobname = f"Charge_{charge}_Res_{str1270}_{str1400}_{str1410}_{str1430}_{str1680}_Lambda_{lambda_str}_VarSqS_{sqrt_s_min_str}_to_{sqrt_s_max_str}"

# Specify input and output files
input_file1 = os.path.join(file_dir, f"{jobname}_J2_{date_str}.dat")
input_file2 = os.path.join(file_dir, f"{jobname}_nJXJc_{date_str}.dat")
input_file3 = os.path.join(file_dir, f"{jobname}_CosTheta_{date_str}.dat")
output_file = os.path.join(file_dir, f"{jobname}_omega_{date_str}.dat")

# Open the input files
with open(input_file1, 'r') as file_a, open(input_file2, 'r') as file_b, open(input_file3, 'r') as file_c:
    # Open the output file
    with open(output_file, 'w') as file_output:
        # Iterate over each line in the input files
        for line_a, line_b, line_c in zip(file_a, file_b, file_c):
            # Read the contents of each line
            a = float(line_a.strip().split(",")[0].strip("()"))
            b = float(line_b.strip().split(",")[1].strip("()"))
            c = float(line_c.strip())

            # Calculate the result
            result = (b * c) / (a * (1 + c**2))

            # Write the result to the output file
            #file_output.write(str(result) + "\n")
            file_output.write(f'{result:.6f}\n')

