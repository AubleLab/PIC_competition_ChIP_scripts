 #!/bin/bash
 
# The script goes through subdirectories in the current directory and maps the 
# fastq files (ending 1.fq.gz / 2.fq.gz) - paired-end to S. pombe genome

# get path to index files
indexFiles="path/hisat2_index/genome"

# get the name of current directory
currentDir=$(pwd)

# create folder in current directory, where final bam files will be stored
mkdir BAM_Spombe
bamDirectory="${currentDir}/BAM_Spombe"


# go through all subdirectiories
for d in */  
do
	# do not run fastqc in the subdirectory with the final BAM files
	if [ "$d" != "BAM_Spombe/" ]
    then
        echo "$d"
    
    	# change path to subdirectory
    	cd $d
    	
    	# find mate 1 - file ending with 1.fq.gz -> remove path extension
    	mate1="$(find . -type f -name *1.fq.gz)"
    	mate1File="${mate1##*/}"
    	#echo "$mate1File"
    	
    	# find mate 2 - file ending with 2.fq.gz -> remove path extension
    	mate2="$(find . -type f -name *2.fq.gz)"
    	mate2File="${mate2##*/}"
    	#echo "$mate2File"
    	
    	# create base name for processed files
    	baseName="${mate1File%%_*}"
    	#echo "$baseName"
    	
    	# create name for unsorted and sorted bam file
    	unsortedBAM="${baseName}_Spombe_unsorted.bam"
    	#echo "$unsortedBAM"
    	
    	sortedBAM_prefix="${baseName}_Spombe_sorted"
    	sortedBAM="${baseName}_Spombe_sorted.bam"
    	#echo "$sortedBAM"
    	
    	BAI="${sortedBAM}.bai"
    	#echo "$BAI"
    	
    	
    	# now let's map the fastq files with hisat2 - paired reads and NEBNext Ultra kit forward stranded
    	# don't make special alignment that will get passed to stringtie - for yeast we will just use featureCount
    	# - use 7 threads
    	hisat2 -p 7 --rna-strandness RF -x $indexFiles -1 $mate1File -2 $mate2File | samtools view -S -b - > $unsortedBAM
    	
    	echo -e Sorting and indexing BAM files...
    	#sort the bam file - use 7 threads
    	samtools sort -@ 7 $unsortedBAM $sortedBAM_prefix
    	
    	#generate index file
    	samtools index $sortedBAM $BAI
    	
    	
    	# remove unnecessary files and move BAM+BAI files to final destination
    	rm $unsortedBAM
    	mv $sortedBAM $bamDirectory
    	mv $BAI $bamDirectory
    
    	# go back to the original directory
    	cd ..
    fi
    
done


echo "  "
echo "Done!"
