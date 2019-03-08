#!/bin/bash

## Define usage function
usage () {
	echo '
     Description:
	 This script extracts the normalization factors from a hic file for all chromosomes and resolutions.

     Dependencies:
	 java
	 juicer tools (aiden lab)

     Usage:
	-f : path/to/.hic
	-j : path/to/juicer_tools.jar
	-n : normalization method
	-o : file/out/path/.txt

     Example:
	./extractNorm.sh -f MN_inter.hic -j juicer_tools.1.9.9_jcuda.0.8.jar -n KR -o normFactors.txt
	
' #| less

}

## Define default values
file_DEFAULT="inter_30.hic"
nrm_DEFAULT="KR"
juicer_tools_DEFAULT="juicer_tools.1.9.9_jcuda.0.8.jar"
output_file_DEFAULT="."

## Parse command-line arguments with getopts
while getopts f:j:n:o: ARGS;
do
	case "${ARGS}" in
		f)
		   file=${OPTARG};;
		j)
		   juicer_tools=${OPTARG};;
		n)
		   nrm=${OPTARG};;
		o)
		   output_file=${OPTARG};;
		*)
		   usage
		   exit 1
		;;
	esac
done

# VARIABLE when it is returned
: ${file=$file_DEFAULT}
: ${juicer_tools=$juicer_tools_DEFAULT}
: ${nrm=$nrm_DEFAULT}
: ${output_file=$output_file_DEFAULT}

shift $((OPTIND -1))

# Delcare all chromosomes and resolutions
declare -a chrs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 "X")
declare -a res=(5000 10000 25000 50000 100000 250000 500000 1000000 2500000)
filebase=$output_file

# Loop through all resolutions and chromomosomes with a new file for each resolution (in parallel)
for i in "${res[@]}"
do
	(filename=$i$filebase
	for j in "${chrs[@]}"
	do
		awk -v OFS="\t" -v a="$i" -v b="$j" 'BEGIN {print "vector", "normFactors", "chr"b, a, "BP"}' >> $filename
	   	java -jar $juicer_tools dump norm $nrm $file "$j" BP "$i" >> $filename
	done)&
	wait
done

echo "Done!"


