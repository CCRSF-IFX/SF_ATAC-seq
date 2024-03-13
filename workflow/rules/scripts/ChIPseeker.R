
library(ChIPseeker)
#library(args[1])
library("ggplot2")
library(clusterProfiler)
library(ReactomePA)
args <- commandArgs(trailingOnly = TRUE)


#args <- list("hg38", "Genrich/GH3564.narrowPeak", "Genrich/GH3565.narrowPeak", "Genrich/GH3566.narrowPeak", "Genrich/GH3568.narrowPeak", "Genrich/GH3562.narrowPeak", "Genrich/GH3560.narrowPeak", "Genrich/GH3563.narrowPeak", "Genrich/GH3570.narrowPeak", "Genrich/GH3561.narrowPeak", "Genrich/GH3571.narrowPeak", "Genrich/GH3567.narrowPeak", "Genrich/GH3569.narrowPeak")
#genome <- args[1]
genome <- snakemake@params['ref']

#length <-length(args[2:length(args)])
print (genome)
#print (length)

if (genome == "mm10") {
  library(org.Mm.eg.db)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  orgdb = "org.Mm.eg.db"
  organism = "mmu"
}else if(genome == "hg38" ){
  library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  orgdb = "org.Hs.eg.db"
  organism = "hsa"
}else if(genome == "hg19" ){
    library(org.Hs.eg.db)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    orgdb = "org.Hs.eg.db"
    organism = "hsa"
}else if(genome == "rheMac10" ){
  library(org.Mmu.eg.db)
  library(TxDb.Mmulatta.UCSC.rheMac10.refGene)
  txdb <- TxDb.Mmulatta.UCSC.rheMac10.refGene
  orgdb = "org.Mmu.eg.db"
  organism = "Mmu"
}

input <- snakemake@input
print (input)
#input <- args[2:length(args)]
input1 <- input[sapply(input, file.size) > 0]
#input <- lapply(strsplit(args[2:length(args)], split='/', fixed=TRUE), function(x) (x[2]))
nm <- lapply(gsub('.narrowPeak', '' , lapply(gsub('Genrich/', '', input1), function(x) (x))), function(y) (y))

print (input1)
names(input1) <- nm

#new_input <- paste0("Genrich",sep="/", list.files('Genrich', 'narrowPeak$'))
#new_nm <- lapply(gsub('.narrowPeak', '' , lapply(gsub('Genrich/', '', new_input), function(x) (x))), function(y) (y))
#names(new_input) <- new_nm

#input = list('N1_N2.narrowPeak', 'T1_T2.narrowPeak')
#names(input) <- c('N1_N2', 'T1_T2')

#peak <- readPeakFile('N1_N2.narrowPeak', header = FALSE)

#head(peak)
#annotatePeak('peak')
#covplot(peak, weightCol="V7")

#########################################################
## Profile of ChIP peaks binding to TSS regions
#########################################################
#pdf('ChIPseeker.pdf', onefile = TRUE)
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
promoter
#tagMatrix <- getTagMatrix(peak, windows=promoter)
#tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
tagMatrixList_all <- lapply(input1, getTagMatrix, windows=promoter)
#print (tagMatrixList)
tagMatrixList <- tagMatrixList_all[sapply(tagMatrixList_all, nrow) > 0]
#########################################################
## Average Profile of ChIP peaks binding to TSS region
#########################################################

if (is.list(tagMatrixList) & length(tagMatrixList) != 0) {
  png('ChIPseeker/Average_Profile_of_ATAC_peaks_binding_to_TSS_region.png',res = 200, width = 6, height = 4, units = "in")
  #par(mar=c(1,1,1,1))
  print(plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Frequency"))
  dev.off()
}
  #png('ChIPseeker/Average_Profile_of_ATAC_peaks_binding_to_TSS_region_CI_mqc.png',res = 200, width = 6, height = strwidth(nm[1], units='inches') * 1.5 * length, units = "in")
  #print(plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row"))
  #dev.off()

#########################################################
## Heatmap
#########################################################

  #png('ChIPseeker/Heatmap_of_ATAC_binding_to_TSS_regions_mqc.png',res = 200, width = strwidth(nm[1], units = 'inches')* 1.5 * length, height =  strwidth(nm[1], units = 'inches')* 1.5 * length * 2, units = "in")
  #par(mar=c(1,1,1,1))
  #par("mar")
  #print(tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL))
  #dev.off()

#########################################################
## Peak Annotation
#########################################################
png('ChIPseeker/Genomic_Annotation_among_different_ATACseq_data_mqc.png',res = 200, width = 6, height = 4, units = "in")
peakAnnoList <- lapply(input1, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList)
dev.off()

#########################################################
## Visualize distribution of TF-binding loci relative to TSS
#########################################################

png('ChIPseeker/Distribution_of_Binding_Sites_among_different_ATACseq_data_mqc.png', res = 200, width = 8, height = 4, units = "in")
plotDistToTSS(peakAnnoList)
dev.off()


##############################################################
## Functional profiles comparison
##############################################################


#genes1 = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
#names(genes1) = sub("_", "\n", names(genes1))
#compKEGG <- compareCluster(geneCluster = genes1, fun = "enrichKEGG", pvalueCutoff  = 0.05, pAdjustMethod = "BH")
#ht = max(7, min(ceiling(sum(compKEGG@compareClusterResult$Count > 0)/4)))
#wt = max(ceiling(max(nchar(compKEGG@compareClusterResult$Description))/8) + (length/2), 7)
#png('ChIPseeker/Functional_profiles_KEGG_comparison_mqc.png',res = 600, width = wt, height = 10, units = "in")
#dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")
#dev.off()

#compGO <- compareCluster(geneCluster = genes1, fun= "enrichGO", pvalueCutoff  = 0.05, pAdjustMethod = "BH", OrgDb = orgdb)
#ht = max(7, min(ceiling(sum(compGO@compareClusterResult$Count > 0)/4)))
#wt = max(ceiling(max(nchar(compGO@compareClusterResult$Description))/8) + (length/2), 7)
#png('ChIPseeker/Functional_profiles_GO_comparison.png',res = 200, width = wt, height = 6, units = "in")
#p <- dotplot(compGO, showCategory = 15, title = "GO Pathway Enrichment Analysis", font.size = 8)
#p + theme(axis.text.x=element_text(angle=90, hjust=1))
#dev.off()


##############################################################
## Venn Diagram
##############################################################

#png('ChIPseeker/Overlap_of_peaks_and_annotated_genes_mqc.png', res = 300, width = 10, height = 6, units = "in")
#genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
#vennplot(genes)

#dev.off()

##############################################################
## Upset plot
##############################################################

#png("Visualize_Genomic_Annotation_Upset_Pie_mqc.png" , res = 300, width = 10, height = 6, units = "in")
#upsetplot(peakAnnoList,vennpie=TRUE)
#dev.off()

#########################################################
## Average Profile of ChIP peaks binding to TSS region
#########################################################

#plotAvgProf(tagMatrix, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

#plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)

#plotDistToTSS(peakAnno,title="Distribution of transcription factor-binding loci\nrelative to TSS")

#########################################################
## Peak Annotation
#########################################################

#peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db")
#plotAnnoBar(peakAnno)

##############################################################
## Visualize distribution of TF-binding loci relative to TSS
##############################################################

#plotDistToTSS(peakAnno,title="Distribution of transcription factor-binding loci\nrelative to TSS")

##############################################################
## Functional profiles comparison
##############################################################

#genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
#names(genes) = sub("_", "\n", names(genes))
#compKEGG <- compareCluster(geneCluster   = genes, fun = "enrichKEGG", pvalueCutoff  = 0.05, pAdjustMethod = "BH")
#dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")
