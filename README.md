# Genome-scale chromatin interaction dynamic measurements for key components of the RNA Pol II general transcription machinery
Scripts accompanying manuscript


Scripts for mapping raw reads to *sacCer3* index and peak calling are placed in the */scripts_preprocessing* folder.
### 1. Map fastq files
From the directory with the fastq.gz files run `mapFastq.sh` script. !!! Line 8 pointing to the sacCer
index must be changed to the local path to the sacCer3 index!!!

### 2. Call peaks
After running `peakCalling.sh` script you will be asked following questions:

*What is the name of the experimental dataset?*
Provide full path to the BAM file from which you want to call peaks.

*What is the name of the control dataset?*
Provide full path to the BAM file from control input.

*What prefix would you like added to the output files?*
Provide name you want for your peak files.
