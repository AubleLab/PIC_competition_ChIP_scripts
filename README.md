# Genome-scale chromatin interaction dynamic measurements for key components of the RNA Pol II general transcription machinery

 Scripts accompanying manuscript

## ___________________ COMPETITION CHIP ___________________

## Few preprocessing steps

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

## Estimate residence times

The scripts for estimating residence times are placed in the folder */scripts_residenceTimeEstimates*. The scripts are using normalized count tables that are placed in the folder */data/normalized_count_tables*, normalized Western blotting measurements placed from */data/westerns.csv* or alternatively already fitted data from Western blotting: */data/westerns_fitParam.csv*.

**Stefan Scripts**

### 2. Additional reliably fast sites


## Scripts and data to create figures

All R scripts for analysis of preprocessed files are in the */scripts_figures* folder. These rely on data placed in the */data* folder. Each script name starts with a prefix indicating the figure and panel name in the manuscript. Before running the scripts set the working directory path to the global *PIC_competition_ChIP_scripts* directory.


## ___________________ NASCENT RNA-SEQ ___________________
