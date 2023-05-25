#!/bin/bash
# Script for generating STAR genome references

USAGE="	USAGE:
$0 [-h] [-g name -d dir] [-l] 

	-g path to the reference genome file
	-d path to the output directory for genome reference
"

# get command line arguments
while getopts g:d:l:h option
do
case "${option}"
in
g) GENOME_FILE=${OPTARG};;
d) GENOME_DIR=${OPTARG};;
h) echo "$USAGE"
       exit;;
esac
done


# test if input paths are not empty
if [ "$GENOME_FILE" == "" ]; then
	echo "Empty name of the genome file!"
	echo "$USAGE"
	exit 1
fi

if [ "$GENOME_DIR" == "" ]; then
	echo "Empty path to the genome directory!"
	echo "$USAGE"
	exit 1
fi

if [ ! -d "$GENOME_DIR" ]; then
	mkdir -p "$GENOME_DIR"
fi

if [ ! -f "$GENOME_FILE" ]; then 
	echo "Input file '$GENOME_FILE' doesn't exist"
	exit 1
fi


# count genomeSAindexNbases parameter for STAR
LENGTH=$(tail -n +2 $GENOME_FILE | wc -m)
echo "$LENGTH"
BASE=$(python -c "import math;print (round(min(14, math.log2($LENGTH)/2 - 1)) - 1)")
echo "$BASE"

STAR --runMode genomeGenerate --genomeDir $GENOME_DIR --genomeFastaFiles $GENOME_FILE --genomeSAindexNbases $BASE



