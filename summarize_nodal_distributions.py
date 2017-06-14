#!/usr/bin/env python3

### This is the main script for processing MCMC output of Ambitus for plotting nodal reconstructions on maps. There are also subsidiary R modules for processing dataframes.
# Intermediate files from Ambitus must be saved (decline to delete at the end of the script, or ignore the question and kill the program).
# Copy these scripts and all relevant "reconstructionclean" files from the temporary folder of a finished ambitus run to an empty folder before running these. For instance, if PNOs 1-4 are needed, only copy the "reconstructionclean" files numbered 1-4.
# This main script can be run in one go, but it may be helpful to copy and paste each step into a terminal and examine output.

import subprocess # To call on bash
import sys # To process arguments
import argparse # Parse arguments
import os

print("Ambitus projector -- software for Bayesian ancestral niche reconstruction.\n")

# Parse command line arguments
parser = argparse.ArgumentParser(description='Script to delete unnecessary fields in locality documets that are tab-delimited and follow DarwinCore header terms.')
parser.add_argument('variable', action='store', help='Number of niche variables.')
parser.add_argument('ambitus_path', action='store', help='Path to desired ambitus run.')
parser.add_argument('node_number', action='store', help='Number of nodes -- can get from the output or from the nodelist file that ambitus produces.')
parser.add_argument('--list', nargs='+', help='List of climate variables.', required=True)
args = parser.parse_args()

# Define variables
process = int(args.variable)
ambitus_path = args.ambitus_path
node_number = int(args.node_number)

def shell_call(command):
	p = subprocess.Popen(command, shell=True)
	p.wait()
	sys.stdout.flush() 

# Concatenate individual MCMC chains across samples
# This step could fail for very large numbers of PNO samples -- this will vary depending on ARG_MAX -- in this case the following commands should be rewritten using find
for i in range(1, process + 1):
	shell_call("rm ./reconstruction_pnoclean{1}.tsv".format(ambitus_path, i))
	shell_call("for f in {0}/output/character{1}_boot*.RECONSTRUCTION.log; do sed -n '/^Iteration\t/,$p' {2} | tail -n+2 >> ./reconstruction_pnoclean{1}.tsv; done".format(ambitus_path, i, "${f}")) # Sed is to print header line and after, tail is to remove header

# The output can be extremely large, a problem that is handled by thinning the MCMC.
# This command will subsample every 50th sample.
# With the BayesTraits settings for the Heuchera paper, this results in thinning to 100 samples per MCMC
for i in range(1, process + 1):
	shell_call("awk 'NR == 1 || NR % 50 == 0' reconstruction_pnoclean{0}.tsv > reconstruction_pnoclean{0}_subsampled.tsv".format(i))

# Need to create correctheader.txt since the files used here are headerless.
# The correct_header is made automatically RECONSTRUCTION log file from the output directory of ambitus to obtain a header for the reconstruction files, as the only line of correctheader.txt.
# The original header is the line immediately before MCMC output begins in the RECONSTRUCTION log.
# This version of the script builds maps for all possible nodes for all environmental variables
# Rename the relevant columns of the header with "hybrid," "source" in correctheader.txt

for i in range(1, process + 1):
	shell_call("grep -P '^Iteration\t' {0}/output/character1_boot1_runID*.RECONSTRUCTION.log > ./correctheader.txt".format(ambitus_path))
	shell_call("sed -i 's/ - 1\t/_1\t/g' correctheader.txt")
	shell_call("sed -i 's/\tEst node/\tEst_node/g' correctheader.txt")
	shell_call("cat correctheader.txt reconstruction_pnoclean{0}_subsampled.tsv > PNO_distribution_{0}.tsv".format(i))

# The tsv files are the input for either the binary map or the binned probabilities map.

##################
# BINARY MAP CODE #
##################

# Create shell scripts for GDAL to calculate binary maps for each variable

for j in range(1, node_number + 1):
	for i in range(1, process + 1):
		for scenario in args.list:
			shell_call("echo '#!/usr/bin/env Rscript' > pno_binary_calculator_node{0}_pno{1}_{2}.r".format(j, i, scenario))
			shell_call("echo 'write(\"#!/usr/bin/env bash\", file = \"BINARY_classify_pixels_node{0}_pno{1}_{2}.sh\")' >> pno_binary_calculator_node{0}_pno{1}_{2}.r".format(j, i, scenario))
			shell_call("echo 'csv = read.csv(\"PNO_distribution_{0}.tsv\", header = TRUE, sep = \"\t\")' >> pno_binary_calculator_node{1}_pno{0}_{2}.r".format(i, j, scenario))
			shell_call("echo 'source_lower <- mean(csv$Est_node{0}_1) - (sd(csv$Est_node{0}_1) * 2)' >> pno_binary_calculator_node{0}_pno{1}_{2}.r".format(j, i, scenario))
			shell_call("echo 'source_upper <- mean(csv$Est_node{0}_1) + (sd(csv$Est_node{0}_1) * 2)' >> pno_binary_calculator_node{0}_pno{1}_{2}.r".format(j, i, scenario))
			# The following line was difficult to escape properly
			shell_call('''echo "output <- paste('gdal_calc.py -A ./scenarios/{0}{1}.tif --outfile=binary_node{2}_pno{1}_{0}.tif --calc \\"logical_and((A >= \', source_lower, \'), (A <=  \', source_upper, \'))\\"  --type=Float32 --NoDataValue=0 --overwrite', sep = \\"\\")" >> pno_binary_calculator_node{2}_pno{1}_{0}.r'''.format(scenario, i, j))
			shell_call("echo 'write(output, file = \"BINARY_classify_pixels_node{0}_pno{1}_{2}.sh\", append = TRUE)' >> pno_binary_calculator_node{0}_pno{1}_{2}.r".format(j, i, scenario))
			shell_call("chmod a+x pno_binary_calculator_node{0}_pno{1}_{2}.r".format(j, i, scenario))
			shell_call("./pno_binary_calculator_node{0}_pno{1}_{2}.r".format(j, i, scenario))
			shell_call("chmod a+x BINARY_classify_pixels_node{0}_pno{1}_{2}.sh".format(j, i, scenario))
			shell_call("./BINARY_classify_pixels_node{0}_pno{1}_{2}.sh".format(j, i, scenario))
	
	# Build gdal command strings			
	for scenario in args.list:
		letter_list = []
		file_list = []
		for i in range(1, process + 1):
			k = str(chr(i + 96)) # Since a is ascii 97
			letter_list.append(k.upper())
			l = "binary_node{0}_pno{1}_{2}.tif".format(j, i, scenario)
			file_list.append(l)
				
		letter_string = "*".join(letter_list)
		print(letter_string)
					
		raster_list = []
		for y, z in zip(letter_list, file_list):
				raster_list.append("".join(["-", y, " ", z]))
			
		raster_string = " ".join(raster_list)
		print(raster_string)
		
		# Combine binary maps for each variable			
		print("gdal_calc.py {0} --outfile BINARY_final_range_product_node{1}_{2}.tif --calc=\"{3}\"".format(raster_string, j, scenario, letter_string))
	
		shell_call("gdal_calc.py {0} --outfile BINARY_final_range_product_node{1}_{2}.tif --calc=\"{3}\"".format(raster_string, j, scenario, letter_string))
	
	# Clean up intermediate files
	for i in range(1, process + 1):
		for scenario in args.list:
			shell_call("rm binary_node{0}_pno{1}_{2}.tif".format(j, i, scenario))
			shell_call("rm BINARY_classify_pixels_node{0}_pno{1}_{2}.sh".format(j, i, scenario))
			shell_call("rm pno_binary_calculator_node{0}_pno{1}_{2}.r".format(j, i, scenario))

##############################
## BINNED PROBABILITIES CODE #
##############################

for j in range(8, node_number + 1):
	for i in range(1, process + 1):
		for scenario in args.list:
			shell_call("echo '#!/usr/bin/env Rscript' > pno_binned_probabilities_calculator_pno{0}_node{1}.r".format(i, j))
			shell_call("echo 'csv = read.csv(\"PNO_distribution_{0}.tsv\", header = TRUE, sep = \"\t\")' >> pno_binned_probabilities_calculator_pno{0}_node{1}.r".format(i, j))
			shell_call("echo 'breaks <- seq(from = min(csv$Est_node{0}_1), to = max(csv$Est_node{0}_1), length.out = 51)' >> pno_binned_probabilities_calculator_pno{1}_node{0}.r".format(j, i))
			shell_call("echo 'h <- hist(csv$Est_node{0}_1, breaks = breaks, plot=FALSE)' >> pno_binned_probabilities_calculator_pno{1}_node{0}.r".format(j, i))
			shell_call("echo 'h$counts=h$counts/sum(h$counts)' >> pno_binned_probabilities_calculator_pno{0}_node{1}.r".format(i, j))
			shell_call("echo 'breaklimits <- embed(h$breaks, 2)[, 2:1]' >> pno_binned_probabilities_calculator_pno{0}_node{1}.r".format(i, j))
			shell_call("echo 'data <- data.frame(breaks = breaklimits, probability = h$counts)' >> pno_binned_probabilities_calculator_pno{0}_node{1}.r".format(i, j))
			shell_call("echo 'write.csv(data, file = \"bins_pno{0}_node{1}.csv\")' >> pno_binned_probabilities_calculator_pno{0}_node{1}.r".format(i, j))
			shell_call("chmod a+x pno_binned_probabilities_calculator_pno{0}_node{1}.r".format(i, j))
			shell_call("./pno_binned_probabilities_calculator_pno{0}_node{1}.r".format(i, j))

			shell_call("echo '#!/usr/bin/env Rscript' > pno_binned_probabilities_calculator_step2_pno{0}_node{1}_{2}.r".format(i, j, scenario))
			shell_call('''echo "write('#!/usr/bin/env bash', file = 'GDAL_classify_pixels_pno{0}_node{1}_{2}.sh')" >> pno_binned_probabilities_calculator_step2_pno{0}_node{1}_{2}.r'''.format(i, j, scenario))
			shell_call('''echo "csv = read.csv(\\"bins_pno{0}_node{1}.csv\\", header = TRUE, sep = \\",\\")" >> pno_binned_probabilities_calculator_step2_pno{0}_node{1}_{2}.r'''.format(i, j, scenario))
			shell_call('''echo 'for (i in 1:nrow(csv)) {{' >> pno_binned_probabilities_calculator_step2_pno{0}_node{1}_{2}.r'''.format(i, j, scenario)) # Double {{ is to escape the character
			shell_call('''echo "output <- paste('gdal_calc.py -A ./scenarios/{0}{1}.tif --outfile=binned_pno{1}_node{2}_{0}\', i, \'.tif --calc \\"logical_and((A >= \', csv\$breaks.1[i], \'), (A <  \', csv\$breaks.2[i], \'))*\', csv\$probability[i], \'\\" --type=Float32 --NoDataValue=0 --overwrite', sep = \\"\\")" >> pno_binned_probabilities_calculator_step2_pno{1}_node{2}_{0}.r'''.format(scenario, i, j))		
			shell_call('''echo "write(output, file = 'GDAL_classify_pixels_pno{0}_node{1}_{2}.sh', append = TRUE)" >> pno_binned_probabilities_calculator_step2_pno{0}_node{1}_{2}.r'''.format(i, j, scenario))
			shell_call('''echo "output = ''" >> pno_binned_probabilities_calculator_step2_pno{0}_node{1}_{2}.r'''.format(i, j, scenario))
			shell_call("echo '}}' >> pno_binned_probabilities_calculator_step2_pno{0}_node{1}_{2}.r".format(i, j, scenario)) # Double }} is to escape the character
			shell_call("chmod a+x pno_binned_probabilities_calculator_step2_pno{0}_node{1}_{2}.r".format(i, j, scenario))
			shell_call("./pno_binned_probabilities_calculator_step2_pno{0}_node{1}_{2}.r".format(i, j, scenario))
			shell_call("chmod a+x GDAL_classify_pixels_pno{0}_node{1}_{2}.sh".format(i, j, scenario))
			shell_call("./GDAL_classify_pixels_pno{0}_node{1}_{2}.sh".format(i, j, scenario))

	for i in range(1, process + 1):
		for scenario in args.list:
			shell_call("mkdir binned_pno{0}_node{1}_{2}/".format(i, j, scenario))
			shell_call("mv binned_pno{0}_node{1}_{2}*.tif binned_pno{0}_node{1}_{2}/".format(i, j, scenario))

	for i in range(1, process + 1):
		for scenario in args.list:
			shell_call("gdalbuildvrt virtual_pno{0}_node{1}_{2}.vrt binned_pno{0}_node{1}_{2}/*.tif".format(i, j, scenario))
			shell_call("gdal_translate -of GTiff virtual_pno{0}_node{1}_{2}.vrt combined_pno{0}_node{1}_{2}.tif".format(i, j, scenario))
	
	for scenario in args.list:
		file_list = []
		for i in range(1, process + 1):
			letter_list = []
			file_list = []
			for i in range(1, process + 1):
				k = str(chr(i + 96)) # Since a is ascii 97
				letter_list.append(k.upper())
				l = "combined_pno{0}_node{1}_{2}.tif".format(i, j, scenario)
				file_list.append(l)
				
			letter_string = "*".join(letter_list)
			print(letter_string)
						
			raster_list = []
			for y, z in zip(letter_list, file_list):
					raster_list.append("".join(["-", y, " ", z]))
				
			raster_string = " ".join(raster_list)
			print(raster_string)
			
			# Multiply across variables for combined probabilities			
			shell_call("gdal_calc.py {0} --outfile final_range_product_node{1}_{2}.tif --calc=\"{3}\"".format(raster_string, j, scenario, letter_string))
		
	for i in range(1, process + 1):
		for scenario in args.list:
			shell_call("rm binned_pno{0}_node{1}_{2}/*; rmdir binned_pno{0}_node{1}_{2}/".format(i, j, scenario))

