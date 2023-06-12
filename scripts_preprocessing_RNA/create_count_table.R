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
                            pathtoGTFfile="/Users/lilcrusher/annotations/yeast/Saccharomyces_cerevisiae.R64-1-1.106.gtf", 
                            countedFeature="gene", 
                            attributeType="gene_id",
                            pathToChrAliases="/Users/lilcrusher/annotations/yeast/chrAliasesForFeatureCounts.txt",
                            groupFeatures=TRUE,
                            fractionMultimapping=FALSE,
                            stranded=1,
                            fastaFile="/Users/lilcrusher/annotations/yeast/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa", 
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
                             isPairedEnd=TRUE,
                             nthreads=12)
  
}

# # 1)  count cerevisiae - whole genes vs. exons only --wrong strandness for our data
# countCerevisiae = createCountTable(pathToBAMFiles="/Users/lilcrusher/yeast_RNA/BAM",
#                                    countedFeature="gene")
# countCerevisiae_final = countCerevisiae$counts
# write.csv(countCerevisiae_final, "countTables/countCerevisiae_genes.csv", quote = F)
# 
# 
# 
# countCerevisiae_exonsOnly = createCountTable(pathToBAMFiles="/Users/lilcrusher/yeast_RNA/BAM",
#                                    countedFeature="exon")
# countCerevisiae_exonsOnly_final = countCerevisiae_exonsOnly$counts
# write.csv(countCerevisiae_exonsOnly_final, "countTables/countCerevisiae_mergedExons.csv", quote = F)


# #2)  count cerevisiae DEDUPLICATED - whole genes vs. exons only
# countCerevisiae_deduplicated = createCountTable(pathToBAMFiles="/Users/lilcrusher/yeast_RNA/BAM_deduplicated",
#                                                 countedFeature="gene")
# countCerevisiae_deduplicated_final = countCerevisiae_deduplicated$counts
# write.csv(countCerevisiae_deduplicated_final, "countTables/countCerevisiae_deduplicated_genes.csv", quote = F)
# 
# countCerevisiae_exonsOnly_deduplicated = createCountTable(pathToBAMFiles="/Users/lilcrusher/yeast_RNA/BAM_deduplicated",
#                                                           countedFeature="exon")
# countCerevisiae_exonsOnly_deduplicated_final = countCerevisiae_exonsOnly_deduplicated$counts
# write.csv(countCerevisiae_exonsOnly_deduplicated_final, "countTables/countCerevisiae_deduplicated_mergedExons.csv", quote = F)


# # 3) count pombe: genes
# 
# countPombe = createCountTable(pathToBAMFiles="/Users/lilcrusher/yeast_RNA/BAM_Spombe/",
#                               countedFeature="gene",
#                               pathtoGTFfile = "/Users/lilcrusher/annotations/yeast/S_pombe/Schizosaccharomyces_pombe.ASM294v2.53.gtf", 
#                               pathToChrAliases = "/Users/lilcrusher/annotations/yeast/S_pombe/chrAliasesForFeatureCounts.txt",
#                               attributeType = "ID",
#                               fastaFile = "/Users/lilcrusher/annotations/yeast/S_pombe/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa")
# countPombe_final = countPombe$counts
# write.csv(countPombe_final, "countTables/countPombe_genes2222222.csv", quote = F)
# 
# pombeSamples = (colnames(countPombe_final))


#pombeCounts = data.frame(colSums(countPombe_final)) 

#### !!!! -----Reversely stranded-----!!!!!
# 1)  count cerevisiae - whole genes vs. exons only
countCerevisiae = createCountTable(pathToBAMFiles="/Users/lilcrusher/yeast_RNA/BAM",
                                   countedFeature="gene", stranded = 2)
countCerevisiae_final = countCerevisiae$counts
write.csv(countCerevisiae_final, "countTables/countCerevisiae_genes_minusStranded.csv", quote = F)



countCerevisiae_exonsOnly = createCountTable(pathToBAMFiles="/Users/lilcrusher/yeast_RNA/BAM",
                                             countedFeature="exon", stranded = 2)
countCerevisiae_exonsOnly_final = countCerevisiae_exonsOnly$counts
write.csv(countCerevisiae_exonsOnly_final, "countTables/countCerevisiae_mergedExons_minusStranded.csv", quote = F)


# 2) count pombe: genes
countPombe = createCountTable(pathToBAMFiles="/Users/lilcrusher/yeast_RNA/BAM_Spombe/",
                              countedFeature="gene",
                              pathtoGTFfile = "/Users/lilcrusher/annotations/yeast/S_pombe/Schizosaccharomyces_pombe.ASM294v2.53.gtf",
                              pathToChrAliases = "/Users/lilcrusher/annotations/yeast/S_pombe/chrAliasesForFeatureCounts.txt",
                              attributeType = "ID",
                              fastaFile = "/Users/lilcrusher/annotations/yeast/S_pombe/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa",
                              stranded = 2)
countPombe_final = countPombe$counts
write.csv(countPombe_final, "countTables/countPombe_genes_minusStranded.csv", quote = F)
# 
# pombeSamples = (colnames(countPombe_final))
# 
# # 3) count cerevisiae DEDUPLICATED
# countCerevisiae_deduplicated = createCountTable(pathToBAMFiles="/Users/lilcrusher/yeast_RNA/BAM_Spombe_deduplicated",
#                                                 countedFeature="gene", 
#                                                 stranded = 2)
# countCerevisiae_deduplicated_final = countCerevisiae_deduplicated$counts
# write.csv(countCerevisiae_deduplicated_final, "countTables/countCerevisiae_deduplicated_genes_minusStranded.csv", quote = F)
# 
# countCerevisiae_exonsOnly_deduplicated = createCountTable(pathToBAMFiles="/Users/lilcrusher/yeast_RNA/BAM_Spombe_deduplicated",
#                                                           countedFeature="exon", 
#                                                           stranded = 2)
# countCerevisiae_exonsOnly_deduplicated_final = countCerevisiae_exonsOnly_deduplicated$counts
# write.csv(countCerevisiae_exonsOnly_deduplicated_final, "countTables/countCerevisiae_deduplicated_mergedExons_minusStranded.csv", quote = F)


#pombeCounts = data.frame(colSums(countPombe_final)) 

# #### ---- unstranded ------
# countCerevisiae = createCountTable(pathToBAMFiles="/Users/lilcrusher/yeast_RNA/BAM",
#                                    countedFeature="gene", stranded = 0)
# countCerevisiae_final = countCerevisiae$counts
# write.csv(countCerevisiae_final, "countTables/countCerevisiae_genes_unstranded.csv", quote = F)
# 
# 
# 
# # 2) count pombe: genes
# countPombe = createCountTable(pathToBAMFiles="/Users/lilcrusher/yeast_RNA/BAM_Spombe/",
#                               countedFeature="gene",
#                               pathtoGTFfile = "/Users/lilcrusher/annotations/yeast/S_pombe/Schizosaccharomyces_pombe.ASM294v2.53.gtf", 
#                               pathToChrAliases = "/Users/lilcrusher/annotations/yeast/S_pombe/chrAliasesForFeatureCounts.txt",
#                               attributeType = "ID",
#                               fastaFile = "/Users/lilcrusher/annotations/yeast/S_pombe/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa", 
#                               stranded = 0)
# countPombe_final = countPombe$counts
# write.csv(countPombe_final, "countTables/countPombe_genes_unstranded.csv", quote = F)
# 
# pombeSamples = (colnames(countPombe_final))


# #### !!!! -----Reversely stranded-----!!!!! allow multiple overlaps
# # 1)  count cerevisiae - whole genes (exons are usually the same)
# countCerevisiae = createCountTable(pathToBAMFiles="/Users/lilcrusher/yeast_RNA/BAM",
#                                    countedFeature="gene", stranded = 2, multiOlap = TRUE)
# countCerevisiae_final = countCerevisiae$counts
# write.csv(countCerevisiae_final, "countTables/countCerevisiae_genes_minusStranded_multiOlap.csv", quote = F)
# 
# 
# # 2) count pombe: genes
# countPombe = createCountTable(pathToBAMFiles="/Users/lilcrusher/yeast_RNA/BAM_Spombe/",
#                               countedFeature="gene",
#                               pathtoGTFfile = "/Users/lilcrusher/annotations/yeast/S_pombe/Schizosaccharomyces_pombe.ASM294v2.53.gtf",
#                               pathToChrAliases = "/Users/lilcrusher/annotations/yeast/S_pombe/chrAliasesForFeatureCounts.txt",
#                               attributeType = "ID",
#                               fastaFile = "/Users/lilcrusher/annotations/yeast/S_pombe/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa",
#                               stranded = 2, 
#                               multiOlap = TRUE)
# countPombe_final = countPombe$counts
# write.csv(countPombe_final, "countTables/countPombe_genes_minusStranded_multiOlap.csv", quote = F)
# 
# pombeSamples = (colnames(countPombe_final))
# 
# # 3) count cerevisiae DEDUPLICATED
# countCerevisiae_deduplicated = createCountTable(pathToBAMFiles="/Users/lilcrusher/yeast_RNA/BAM_Spombe_deduplicated",
#                                                 countedFeature="gene",
#                                                 stranded = 2, 
#                                                 multiOlap = TRUE)
# countCerevisiae_deduplicated_final = countCerevisiae_deduplicated$counts
# write.csv(countCerevisiae_deduplicated_final, "countTables/countCerevisiae_deduplicated_genes_minusStranded_multiOlap.csv", quote = F)






