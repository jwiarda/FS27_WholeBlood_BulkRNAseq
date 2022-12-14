library(DESeq2)
library(tidyverse)
library(readxl)
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(EnhancedVolcano)
library(gplots)
library(writexl)

### Upload count data and target file:

# Load target file into R:
target <-read_excel("/Users/jayne.wiarda/Desktop/Experiments/FS27_Salmonella/RNAseq_Bulk_WholeBlood/FS27_TargetFile.xlsx")
#View(target) # view full target file
target$FileName <- sub(".*-FS27-", "", target$FileName)  # change file names to match colnames(counts)
target$FileName <- sub("_R1_.*", "", target$FileName)  
head(target)

# Load counts data back into R:
counts <- read.table("/Users/jayne.wiarda/Desktop/Experiments/FS27_Salmonella/RNAseq_Bulk_WholeBlood/featureCounts/FS27_GeneCountTable_LanesSeparate.txt", header = TRUE)
colnames(counts) # make sure columns 2-37 are in same order as entries in 'label' column of target
rownames(counts) <- counts$Geneid
counts <- counts[,-c(1:6)] # remove unneccesary columns (not containing counts data)
colnames(counts) <- sub(".*.FS27.", "", colnames(counts))  # change column names to match target$FileName
colnames(counts) <- sub("_R1_.*", "", colnames(counts))  

identical(colnames(counts), target$FileName) # make sure the sample IDs match (identical = TRUE)

# Assess variability of technical replicates (biological replicates sequenced across two lanes):
sums <- data.frame(colSums(counts))
barplot(sums$colSums.counts.) # check that every two samples have similar reads detected (indicating biological replicate sequenced across two lanes)
cor <- cor(counts)
cor <- melt(cor)
head(cor)
ggplot(data = cor, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = .95, limit = c(.9, 1)) # should see extremely high correlations for every two samples (indicating biological replicate sequenced across two lanes)
# All looks well across technical replicates

# Merge technical replicates:
colnames(counts) <- substr(colnames(counts), 1, 3) # convert column names to only animal ID
counts <- t(rowsum(t(counts), group = colnames(counts), na.rm = T)) # sum counts for columns with identical names

fDat <- counts
sums <- data.frame(colSums(fDat))
barplot(sums$colSums.fDat.) # see total library size per sample

# Further refine target file:
target$treatment <- paste(target$TrtGroup, target$DaysPostInfection, sep = '')
#colnames(target)
target <- target[,c(3, 4, 7, 19)] # subset to only experimental variables of animal ID, treatment, and dpi
target <- distinct(target) # remove duplicate rows
target <- as.data.frame(target)
rownames(target) <- target$AnimalID
target <- target[,c(2:4)]
head(target)

identical(colnames(counts), rownames(target)) # make sure the sample IDs match (identical = TRUE)

# Define treatment factor:
treatment <- factor(target$treatment)
treatment

# Create DESeq2 matrix:
dds<-DESeqDataSetFromMatrix(counts, colData=target, design = ~treatment)
dds

# filter genes
keep<- rowSums(counts(dds))>=5
dds<-dds[keep,]
dds

# apply DESeq2
d.dds<- DESeq(dds)
d.dds

# Create PCA
vsdb<- varianceStabilizingTransformation(d.dds)
plotPCA(vsdb, intgroup=c("TrtGroup"))
plotPCA(vsdb, intgroup=c("DaysPostInfection"))
plotPCA(vsdb, intgroup=c("treatment"))
#plotPCA(vsdb, intgroup=c("treatment")) + stat_ellipse(aes(fill = treatment), level = 0.95, geom = 'polygon', alpha = 0.1)
plotPCA(vsdb, intgroup=c("treatment")) + stat_ellipse(aes(fill = treatment), level = 0.95)

# DGE results
res1<-results(d.dds, contrast=c("treatment","Sal07", "Crl07"),alpha=0.05, lfcThreshold = 0)
summary(res1)

res2<-results(d.dds, contrast=c("treatment","Crl07","Crl02"),alpha=0.05, lfcThreshold = 0)
summary(res2)

res3<-results(d.dds, contrast=c("treatment","Sal02","Crl02"),alpha=0.05, lfcThreshold = 0)
summary(res3)

res4<-results(d.dds, contrast=c("treatment","Sal07","Sal02"),alpha=0.05, lfcThreshold = 0)
summary(res4)
# we see only res1 comparison yielded any DE genes

# volcano plots for each comp:
EnhancedVolcano(toptable = res1,
                title = "Crl07 (left) v Sal07 (right)",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 1,
                ylim = c(0, 11),
                xlim = c(-9, 9)) #+ coord_flip()

EnhancedVolcano(toptable = res2,
                title = "Crl02 (left) v Crl07 (right)",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 1,
                ylim = c(0, 11),
                xlim = c(-9, 9)) #+ coord_flip()

EnhancedVolcano(toptable = res3,
                title = "Crl02 (left) v Sal02 (right)",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 1,
                ylim = c(0, 11),
                xlim = c(-9, 9)) #+ coord_flip()

EnhancedVolcano(toptable = res4,
                title = "Sal02 (left) v Sal07 (right)",
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 0,
                cutoffLineWidth = 0,
                col = c('grey40', 'grey40', 'grey40', 'red3'),
                colAlpha = 1,
                ylim = c(0, 11),
                xlim = c(-9, 9)) #+ coord_flip()

# Make heatmaps
genes <- subset(res2, padj < 0.05)
vsdb_table<- as.data.frame(assay(vsdb))
a<-vsdb_table[rownames(genes),]
my_palette <- colorRampPalette(c("navy", "grey80", "red3"))(n = 299)
heatmap.2(as.matrix(a), Rowv=T, dendrogram= "both",
          scale="row", key=T, keysize=1.5, density.info="none", trace="none", cexCol=0.9, col = my_palette) # heatmap with all samples
a <- a[,-c(11:20)] # remove Sal samples
heatmap.2(as.matrix(a), Rowv=T, dendrogram= "both",
          scale="row", key=T, keysize=1.5, density.info="none", trace="none", cexCol=0.9, col = my_palette) # heatmap with only Crl samples

genes <- subset(res3, padj < 0.05)
vsdb_table<- as.data.frame(assay(vsdb))
a<-vsdb_table[rownames(genes),]
my_palette <- colorRampPalette(c("navy", "grey80", "red3"))(n = 299)
heatmap.2(as.matrix(a), Rowv=T, dendrogram= "both",
          scale="row", key=T, keysize=1.5, density.info="none", trace="none", cexCol=0.9, col = my_palette) # heatmap with all samples
a <- a[,-c(1, 7:11, 17:20)] # remove Crl07 & Sal07 samples
heatmap.2(as.matrix(a), Rowv=T, dendrogram= "both",
          scale="row", key=T, keysize=1.5, density.info="none", trace="none", cexCol=0.9, col = my_palette) # heatmap with only 2dpi samples

genes <- subset(res4, padj < 0.05)
vsdb_table<- as.data.frame(assay(vsdb))
a<-vsdb_table[rownames(genes),]
my_palette <- colorRampPalette(c("navy", "grey80", "red3"))(n = 299)
heatmap.2(as.matrix(a), Rowv=T, dendrogram= "both",
          scale="row", key=T, keysize=1.5, density.info="none", trace="none", cexCol=0.9, col = my_palette) # heatmap with all samples
a <- a[,-c(1:10)] # remove Crl samples
heatmap.2(as.matrix(a), Rowv=T, dendrogram= "both",
          scale="row", key=T, keysize=1.5, density.info="none", trace="none", cexCol=0.9, col = my_palette) # heatmap with only Sal samples

# Append Ensembl IDs:
annot <- read_excel("/Users/jayne.wiarda/Desktop/Experiments/FS27_Salmonella/RNAseq_Bulk_WholeBlood/UpdatedGeneNameListForSus97GTF_06302021_JEW_SKS.xlsx")
annot <- filter(annot, FinalList %in% rownames(dds))# reduce new gene sybmol annotation to include only genes in the filtered gene list for our single-cell dataset
colnames(annot)
annot <- annot[,-c(2:8)]
colnames(annot) <- c('EnsemblID', 'GeneName')

res1$GeneName <- rownames(res1)
res1 <- as.data.frame(res1)
res1 <- merge(res1, annot, by = 'GeneName')

res2$GeneName <- rownames(res2)
res2 <- as.data.frame(res2)
res2 <- merge(res2, annot, by = 'GeneName')

res3$GeneName <- rownames(res3)
res3 <- as.data.frame(res3)
res3 <- merge(res3, annot, by = 'GeneName')

res4$GeneName <- rownames(res4)
res4 <- as.data.frame(res4)
res4 <- merge(res4, annot, by = 'GeneName')

# save DGE results
write_xlsx(res1, '/Users/jayne.wiarda/Desktop/Experiments/FS27_Salmonella/RNAseq_Bulk_WholeBlood/Results_DGE_DESeq2/Results_DGE_Crl07vSal07.xlsx')
write_xlsx(res2, '/Users/jayne.wiarda/Desktop/Experiments/FS27_Salmonella/RNAseq_Bulk_WholeBlood/Results_DGE_DESeq2/Results_DGE_Crl02vCrl07.xlsx')
write_xlsx(res3, '/Users/jayne.wiarda/Desktop/Experiments/FS27_Salmonella/RNAseq_Bulk_WholeBlood/Results_DGE_DESeq2/Results_DGE_Crl02vSal02.xlsx')
write_xlsx(res4, '/Users/jayne.wiarda/Desktop/Experiments/FS27_Salmonella/RNAseq_Bulk_WholeBlood/Results_DGE_DESeq2/Results_DGE_Sal02vSal07.xlsx')

sessionInfo()
#R version 4.2.1 (2022-06-23)
#Platform: aarch64-apple-darwin20 (64-bit)
#Running under: macOS Monterey 12.6

#Matrix products: default
#LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#attached base packages:
#  [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] writexl_1.4.1               gplots_3.1.3                EnhancedVolcano_1.16.0      ggrepel_0.9.2               reshape2_1.4.4              readxl_1.4.1                forcats_0.5.2               stringr_1.5.0              
#[9] dplyr_1.0.10                purrr_0.3.5                 readr_2.1.3                 tidyr_1.2.1                 tibble_3.1.8                ggplot2_3.4.0               tidyverse_1.3.2             DESeq2_1.38.1              
#[17] SummarizedExperiment_1.28.0 Biobase_2.58.0              MatrixGenerics_1.10.0       matrixStats_0.63.0          GenomicRanges_1.50.1        GenomeInfoDb_1.34.4         IRanges_2.32.0              S4Vectors_0.36.1           
#[25] BiocGenerics_0.44.0        

#loaded via a namespace (and not attached):
#  [1] bitops_1.0-7           fs_1.5.2               lubridate_1.9.0        bit64_4.0.5            RColorBrewer_1.1-3     httr_1.4.4             tools_4.2.1            backports_1.4.1        utf8_1.2.2             R6_2.5.1              
#[11] KernSmooth_2.23-20     DBI_1.1.3              colorspace_2.0-3       withr_2.5.0            tidyselect_1.2.0       bit_4.0.5              compiler_4.2.1         cli_3.4.1              rvest_1.0.3            xml2_1.3.3            
#[21] DelayedArray_0.24.0    labeling_0.4.2         caTools_1.18.2         scales_1.2.1           XVector_0.38.0         pkgconfig_2.0.3        dbplyr_2.2.1           fastmap_1.1.0          rlang_1.0.6            rstudioapi_0.14       
#[31] RSQLite_2.2.19         farver_2.1.1           generics_0.1.3         jsonlite_1.8.4         gtools_3.9.4           BiocParallel_1.32.4    googlesheets4_1.0.1    RCurl_1.98-1.9         magrittr_2.0.3         GenomeInfoDbData_1.2.9
#[41] Matrix_1.5-3           Rcpp_1.0.9             munsell_0.5.0          fansi_1.0.3            lifecycle_1.0.3        stringi_1.7.8          MASS_7.3-58.1          zlibbioc_1.44.0        plyr_1.8.8             grid_4.2.1            
#[51] blob_1.2.3             parallel_4.2.1         crayon_1.5.2           lattice_0.20-45        Biostrings_2.66.0      haven_2.5.1            annotate_1.76.0        hms_1.1.2              KEGGREST_1.38.0        locfit_1.5-9.6        
#[61] pillar_1.8.1           geneplotter_1.76.0     codetools_0.2-18       reprex_2.0.2           XML_3.99-0.13          glue_1.6.2             modelr_0.1.10          png_0.1-8              vctrs_0.5.1            tzdb_0.3.0            
#[71] cellranger_1.1.0       gtable_0.3.1           assertthat_0.2.1       cachem_1.0.6           xtable_1.8-4           broom_1.0.1            googledrive_2.0.0      gargle_1.2.1           AnnotationDbi_1.60.0   memoise_2.0.1         
#[81] timechange_0.1.1       ellipsis_0.3.2   