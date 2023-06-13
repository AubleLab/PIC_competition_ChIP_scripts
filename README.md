# Genome-scale chromatin interaction dynamic measurements for key components of the RNA Pol II general transcription machinery

 Scripts accompanying manuscript

## ___________________ RECREATE FIGURES ___________________

The R scripts to recreate all the figures from the manuscript can be found in */scripts_figures* folder. The prefix in the name of the script indicates figure and panel number. All the necessary preprocessed files can be found in the */data* folder (all scripts point to appropriate file within the folder and do not need path changing, when working directory is set to the main folder).

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

The scripts for TBP fitting: `TBP_Fit_Turnover_Model_RD_Western_Norm_Exact_Init.nb`, TFIIA: `TFIIA_Fit_Turnover_Model_RD_Western_Norm_Exact_Init.nb`, TFIIB: `TFIIB_Fit_Turnover_Model_RD_Western_Norm_Exact_Init.nb`, TFIIE: `TFIIE_Fit_Turnover_Model_RD_Western_Norm_Exact_Init.nb`, and TFIIF: `TFIIF_Fit_Turnover_Model_RD_Western_Norm_Exact_Init.nb`.

These are custom Mathamatica scripts developed by Stefan Bekiranov that fit normalized TBP, TFIIA, TFIIB, TFIIF and TFIIE competition ChIP sequencing data to a mass action turnover model as detailed in the Methods section of Kupkova et al.  These scripts are provided as is.  The transcription factor Western data is entered manually.  The header of the normalized competition ChIP ratio data is stripped before reading it into the Mathematica environment.  The input and output directories have to be modified in a user specific manner.  Determining the best integer value for the Hill model fit to the transcription factor Western data is the first step in the overall fitting procedure.  The best integer value for the Hill model fit to the Western data will be used to fit the Hill model to competition ChIP sequencing data.  It will also be used in the mass action turnover model via the function cB[t].  The initial guess parameters t1 for the Hill model encoded in the function inductionModel[X, t1, t] (4th line of For loop) and resTimeInit (15th line of For loop) for the ratio of induced over endogenous transcription factor occupancies derived from the turnover model and encoded in modelb[ka, kd]/modela[ka, kd] are iteratively determined/optimized after at least one initial round of fitting as detailed in Kupkova et al.  The integer value and t1 initial guess of the Hill model and resTimeInit are the only parameters that are modified/tuned for fitting the mass action turnover model to competition ChIP sequencing data.  While ka (the on-rate) is fit to competition ChIP data, its value is unreliable.  The fitted value of kd (the off-rate) and associated statistics including estimates of the error and significance of kd as well as the overall fit of competition ChIP sequencing data to the ratio of induced over endogenous transcription factor occupancies derived from the turnover model are written to a table in .tsv format.    


The final results from the fitting can be found in */data/time_estimates_first_iter* folder.

### 2. Additional reliably fast sites

All sites are first fitted with Hill equation using `fit_allGTF.R` script. Results created with this script match the results placed in the */data/Hill_fits* folder. The script `additioanlFastSites.R` takes all the Hill fits (skips all the genes/sites that were reliably fitted in the first iteration) and identifies sites with delta Km (difference in estimated half-lives between fitted Western blots and competition ChIP) smaller than 2. The list of these sites, which we all classify as < 1min is then placed to / can be found in */data/time_estimates_second_iter/additionalFastSites.csv*.

**All final residence time estimates can be found in /data/residence_times_all.csv**


## Scripts and data to create figures

All R scripts for analysis of preprocessed files are in the */scripts_figures* folder. These rely on data placed in the */data* folder. Each script name starts with a prefix indicating the figure and panel name in the manuscript. Before running the scripts set the working directory path to the global *PIC_competition_ChIP_scripts* directory.


## ___________________ NASCENT RNA-SEQ ___________________

## Preprocessing

The scripts for mapping the paired-end reads from RNA-seq and creating count tables are here: */scripts_preprocessing_RNA*.

### 1. Map to *S. cerevisiae*

Change line 7 in `map_RNA_Scerevisiae.sh` to path pointing to the local sacCer3 hisat2 index file. After that run the script within folder containing folders with fastq.gz files - for each sample have a separate folder. BAM alignment files are created.  

### 2. Map to *S. pombe* - used for normalization

Change line 7 in `map_RNA_Spombe.sh` to path pointing to the local sacCer3 hisat2 index file. After that run the script within folder containing folders with fastq.gz files - for each sample have a separate folder. BAM alignment files are created.

### 3. Count tables

To create a count table with raw counts, use script `create_count_table.R`. In the script based on the location in your local settings, change path pointing to *S. cerevisiae* BAM files (line 7), and *S. cerevisiae* FASTA file (line 27). The final raw counts are also available here: *data/RNA_count_tables/countCerevisiae_genes_minusStranded.csv*.

Normalized count table used for DTA analysis (only samples with thiouracil added) can be created with `count_table_for_DTA.R` script. The final count table is available here: */data/RNA_count_tables/countCerevisiae_genes_minusStranded__ThiouracilYESonly_pombeNormalized.csv*

## Estimating synthesis rates

The scripts for RNA synthesis rates estimates from count tables are here: */scripts_RNA_synthesisRatesEstimates*.

To calculate synthesis rates from both 20 min and 60 min post galactose time points together, use `DTA_analysis.R`. To calculate synthesis rates for each time point, use `DTA_analysis_dynamic.R`.
