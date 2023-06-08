#!/bin/bash
#This script will run macs for YEAST S. cerevisiae samples using the specified 'test' and 'control' BAM
# files in the working directory. $1=test file, $2=control file, $3=output file prefix
echo -e What is the name of the experimental dataset?
read name1
echo -e What is the name of the control dataset?
read name2
echo -e What prefix would you like added to the output files?
read name3


#change --extsize (requires nomodel)
macs2 callpeak -t $name1 -c $name2 -f BAM -g 1.2e7 -n $name3 --nomodel --extsize 147 -B -q 0.01