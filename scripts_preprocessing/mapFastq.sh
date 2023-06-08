#!/bin/bash

#!/bin/bash

# the script maps all fastq files within a folder to hg19, removes unmapped reads
# outputs are sorted and indexed indexed BAM files


# ! set up pathway to bowtie2 sacCer3 index
sacCer3index="path/bowtie2-2.2.6/indexes/sacCer3.bt2/sacCer3"


#unzip files:
for zippedFastq in *.fastq.gz 
do 
	# get the fastq file name and clean file name without any extensions
	fastqFile="${zippedFastq%.*}"
	fileName="${zippedFastq%%.*}"
	
	# first unzip the file
	gunzip $zippedFastq
	
	# Map to hg19 with bowtie2
	echo "Mapping following file:"
	echo $fastqFile
	
	unsortedBAM="${fileName}_unsorted.bam"
	bowtie2 -x $sacCer3index -p 6 --time -U $fastqFile | samtools view -S -b - > $unsortedBAM
	
	#sort and index the bam file
	sorted="${fileName}_sorted"
	sortedBAM="${sorted}.bam"
	sortedBAI="${sortedBAM}.bai"
	samtools sort $unsortedBAM $sorted
	samtools index $sortedBAM
	
	# remove the unsorted BAM file
	rm $unsortedBAM

	#filter unmapped reads
	filteredBAM="${fileName}_filtered.bam"
	echo -e The numbers of mapped and unmapped reads, respectively:
	samtools view -c -F 4 $sortedBAM
	samtools view -c -f 4 $sortedBAM
	samtools view -h -F 4 -b $sortedBAM > $filteredBAM
	
	# remove unfiltered files
	rm $sortedBAM
	rm $sortedBAI

	#sort and index
	sorted="${fileName}_filtered_sorted"
	sortedBAM="${sorted}.bam"
	sortedBAI="${sortedBAM}.bai"
	samtools sort $filteredBAM $sorted
	samtools index $sortedBAM
	

	echo -e Final read count for this dataset:
	samtools view -c -F 4 $sortedBAM
	
	# zip the fastq file back
	gzip $fastqFile
done
