#!/bin/bash

## Define usage function
usage () {
	echo '
     Description:
	 This script extracts reads from a .hic file corresponding to a bedpe file using the straw api.

     Usage:
	-f : path/to/.hic
	-b : path/to/bedpe
	-n : normalization
	-r : resolution
	-o : file/out/path/.txt
	
' #| less

}

## Define default values
file_DEFAULT="inter_30.hic"
bedpe_DEFAULT="loops.txt"
norm_DEFAULT="NONE"
res_DEFAULT="10000"
output_file_DEFAULT="./output_extractCounts.txt"

## Parse command-line arguments with getopts
while getopts f:b:n:r:o: ARGS;
do
	case "${ARGS}" in
		f)
		   file=${OPTARG};;
	   	b)
   		   bedpe=${OPTARG};;
		n)
		   norm=${OPTARG};;
		r)
		   res=${OPTARG};;
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
: ${bedpe=$bedpe_DEFAULT}
: ${norm=$norm_DEFAULT}
: ${res=$res_DEFAULT}
: ${output_file=$output_file_DEFAULT}

shift $((OPTIND -1))

## Clean bedpe file to remove chr prefix
sed 's/chr//g' $bedpe > cleaned.bedpe
bedpe0=$bedpe
bedpe='cleaned.bedpe'

## Get bedpe file into arrays
chr1=($(awk '{print $1}' $bedpe))
x1=($(awk '{print $2}' $bedpe))
x2=($(awk '{print $3}' $bedpe))
chr2=($(awk '{print $4}' $bedpe))
y1=($(awk '{print $5}' $bedpe))
y2=($(awk '{print $6}' $bedpe))

## Loop through each line of the cleaned bedpe file
for i in "${!chr1[@]}"
do
	a=${x1["$i"]}
	b=${y1["$i"]}
	./straw $norm $file ${chr1["$i"]}:${x1["$i"]}:${x2["$i"]} ${chr2["$i"]}:${y1["$i"]}:${y2["$i"]} BP $res | awk -v a=$a -v b=$b '{if($1 == a && $2 == b) print $3}' >> temp
done

## Add results to original bedpe file and send to output file
paste $bedpe0 temp > $output_file

## Remove temporary files
rm cleaned.bedpe temp
