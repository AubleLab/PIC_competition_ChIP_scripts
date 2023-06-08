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

### 1. Residence time estimated using differential equations

The final results from the fitting can be found in */data/time_estimates_first_iter* folder.

### 2. Additional reliably fast sites

All sites are first fitted with Hill equation using `fit_allGTF.R` script. Results created with this script match the results placed in the */data/Hill_fits* folder. The script `additioanlFastSites.R` takes all the Hill fits (skips all the genes/sites that were reliably fitted in the first iteration) and identifies sites with delta Km (difference in estimated half-lives between fitted Western blots and competition ChIP) smaller than 2. The list of these sites, which we all classify as < 1min is then placed to / can be found in */data/time_estimates_second_iter/additionalFastSites.csv*.

**All final residence time estimates can be found in /data/residence_times_all.csv**


## Scripts and data to create figures

All R scripts for analysis of preprocessed files are in the */scripts_figures* folder. These rely on data placed in the */data* folder. Each script name starts with a prefix indicating the figure and panel name in the manuscript. Before running the scripts set the working directory path to the global *PIC_competition_ChIP_scripts* directory.


## ___________________ NASCENT RNA-SEQ ___________________

## Preprocessing

The scripts for mapping the paired-end reads from RNA-seq are here: */scripts_preprocessing_RNA*.

## 1. Map to *S. cerevisiae*

Change line 7 in `map_RNA_Scerevisiae.sh` to path pointing to the local sacCer3 hisat2 index file. After that run the script within folder containing folders with fastq.gz files - for each sample have a separate folder. BAM alignment files are created.  

## 2. Map to *S. pombe* - used for normalization

Change line 7 in `map_RNA_Spombe.sh` to path pointing to the local sacCer3 hisat2 index file. After that run the script within folder containing folders with fastq.gz files - for each sample have a separate folder. BAM alignment files are created.

## Count tables
