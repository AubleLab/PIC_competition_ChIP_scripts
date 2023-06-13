rm(list = ls())

library(Rsubread)

# pathToBAMFiles = string giving path to folder with BAM files
# path to GTF = string giving fill path to GTF file
# countedFeature = coverage over what should be considered? exons? whole genes? etc. - 
#                   give string that corresponds to ingo in GTF file, e.g. "exon"
# pathToChrAliases = path to file that has chromosoma names in GTF file and hista index
# groupFeatures = should e.g. exons be grouped into genes or reported individually?
# fractionMultimapping = if a read maps to more than one location, should it be reported as 1 
#                       at both locations (FALSE) or 1/n, where n = number of sites it mapped to (TRUE)
# stranded = strandess: 0-no, 1-yes, 2-reverse strand
# fastaFile = path to fast file

# If a feature overlaps with two different exon, it's counted once
# multimapping rads are allowed - this can cause numbers not to be integers, because
# the read is then  reported as 1/n, where n = number of sites it mapped to
createCountTable = function(pathToBAMFiles, 
                            pathtoGTFfile="path/Saccharomyces_cerevisiae.R64-1-1.106.gtf", 
                            countedFeature="gene", 
                            attributeType="gene_id",
                            pathToChrAliases="data/RNA_support_files/chrAliasesForFeatureCounts.txt",
                            groupFeatures=TRUE,
                            fractionMultimapping=FALSE,
                            stranded=1,
                            fastaFile="path/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa", 
                            multiOlap=FALSE){
  # list all bam files in the folder
  BAMfiles = list.files(pathToBAMFiles, pattern = "*.bam$",full.names = T)
  
  countTable = featureCounts(files=BAMfiles,
                             annot.ext=pathtoGTFfile,
                             isGTFAnnotationFile=TRUE,
                             GTF.featureType=countedFeature,
                             GTF.attrType=attributeType,
                             chrAliases=pathToChrAliases,
                             useMetaFeatures=groupFeatures,
                             allowMultiOverlap=multiOlap,
                             countMultiMappingReads=TRUE,
                             fraction=fractionMultimapping,
                             strandSpecific=stranded,
                             genome=fastaFile,
                             isPairedEnd=TRUE)
  
}


#### !!!! -----Reversely stranded-----!!!!!
# 1)  count cerevisiae - whole genes vs. exons only
countCerevisiae = createCountTable(pathToBAMFiles="path/yeast_RNA/BAM",
                                   countedFeature="gene", stranded = 2)
countCerevisiae_final = countCerevisiae$counts
write.csv(countCerevisiae_final, "data/RNA_count_tables/countCerevisiae_genes_minusStranded.csv", quote = F)


