library(ggplot2)
library(readxl)

# Create dot plots of selected GO processes 

## Crl02vCrl07
GO <- read_excel('/Users/jayne.wiarda/Desktop/Experiments/FS27_Salmonella/RNAseq_Bulk_WholeBlood/Results_GO/ResultsGO_Crl02vCrl07_NegFC.xlsx')
colnames(GO) <- c('GO_biological_process', 'Number_BackgroundGenes', 'Number_AnalyzedGenes', 
                  'Expected_AnalyzedGenes', 'OverUnder', 'Fold_Enrichment', 'Raw_p_value', 'FDR')
GO <- subset(GO, OverUnder == "+" & FDR < 0.1)
dim(GO) # see no significant biological processes to plot

GO <- read_excel('/Users/jayne.wiarda/Desktop/Experiments/FS27_Salmonella/RNAseq_Bulk_WholeBlood/Results_GO/ResultsGO_Crl02vCrl07_PosFC.xlsx')
colnames(GO) <- c('GO_biological_process', 'Number_BackgroundGenes', 'Number_AnalyzedGenes', 
                  'Expected_AnalyzedGenes', 'OverUnder', 'Fold_Enrichment', 'Raw_p_value', 'FDR')
GO <- subset(GO, OverUnder == "+" & FDR < 0.1)
dim(GO) # see no significant biological processes to plot

## Sal02vSal07
GO1 <- read_excel('/Users/jayne.wiarda/Desktop/Experiments/FS27_Salmonella/RNAseq_Bulk_WholeBlood/Results_GO/ResultsGO_Sal02vSal07_NegFC.xlsx')
colnames(GO1) <- c('GO_biological_process', 'Number_BackgroundGenes', 'Number_AnalyzedGenes', 
                   'Expected_AnalyzedGenes', 'OverUnder', 'Fold_Enrichment', 'Raw_p_value', 'FDR')
GO1 <- subset(GO1, OverUnder == "+" & FDR < 0.1)
dim(GO1) # see we have some significant processes to plot
GO1$Treatment <- rep('Sal02', nrow(GO1))

GO2 <- read_excel('/Users/jayne.wiarda/Desktop/Experiments/FS27_Salmonella/RNAseq_Bulk_WholeBlood/Results_GO/ResultsGO_Sal02vSal07_PosFC.xlsx')
colnames(GO2) <- c('GO_biological_process', 'Number_BackgroundGenes', 'Number_AnalyzedGenes', 
                   'Expected_AnalyzedGenes', 'OverUnder', 'Fold_Enrichment', 'Raw_p_value', 'FDR')
GO2 <- subset(GO2, OverUnder == "+" & FDR < 0.1)
dim(GO2) # see we have some significant processes to plot
GO2$Treatment <- rep('Sal07', nrow(GO2))

GO <- rbind(GO1, GO2) # merge GO results
# View(GO) # Identify up to 10 processes of interest in the top 30 (based on lowest FDRs) for both the PosFC & negFC lists...since NegFC list only had 9 processes, we will take all of those

terms <- c('cell surface receptor signaling pathway (GO:0007166)', 'positive regulation of immune response (GO:0050778)',
           'positive regulation of immune system process (GO:0002684)', 'response to biotic stimulus (GO:0009607)', 
           'immune response-regulating signaling pathway (GO:0002764)', 'innate immune response (GO:0045087)',
           'response to cytokine (GO:0034097)', 'cellular response to cytokine stimulus (GO:0071345)',
           'cytokine-mediated signaling pathway (GO:0019221)', 'cell communication (GO:0007154)',
           'positive regulation of DNA-templated transcription (GO:0045893)', 'positive regulation of biosynthetic process (GO:0009891)',
           'positive regulation of RNA metabolic process (GO:0051254)', 'adaptive immune response (GO:0002250)',
           'positive regulation of cytosolic calcium ion concentration (GO:0007204)', 'positive regulation of cellular biosynthetic process (GO:0031328)',
           'positive regulation of macromolecule biosynthetic process (GO:0010557)', 'lymphocyte activation (GO:0046649)',
           'positive regulation of RNA biosynthetic process (GO:1902680)', 'leukocyte activation (GO:0045321)')
GO <- GO[GO$GO_biological_process %in% terms, ]
GO$GO_biological_process <- factor(GO$GO_biological_process, levels = terms)
ggplot(GO) +
  geom_point(aes(x = GO_biological_process, 
                 y = Treatment,
                 size = Fold_Enrichment,
                 color = FDR)) +
  theme_bw() +
  scale_color_gradient(low = 'red3', high = 'gold') +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))

## Crl02vSal02
GO1 <- read_excel('/Users/jayne.wiarda/Desktop/Experiments/FS27_Salmonella/RNAseq_Bulk_WholeBlood/Results_GO/ResultsGO_Crl02vSal02_NegFC.xlsx')
colnames(GO1) <- c('GO_biological_process', 'Number_BackgroundGenes', 'Number_AnalyzedGenes', 
                  'Expected_AnalyzedGenes', 'OverUnder', 'Fold_Enrichment', 'Raw_p_value', 'FDR')
GO1 <- subset(GO1, OverUnder == "+" & FDR < 0.1)
dim(GO1) # see we have some significant processes to plot
GO1$Treatment <- rep('Crl02', nrow(GO1))

GO2 <- read_excel('/Users/jayne.wiarda/Desktop/Experiments/FS27_Salmonella/RNAseq_Bulk_WholeBlood/Results_GO/ResultsGO_Crl02vSal02_PosFC.xlsx')
colnames(GO2) <- c('GO_biological_process', 'Number_BackgroundGenes', 'Number_AnalyzedGenes', 
                   'Expected_AnalyzedGenes', 'OverUnder', 'Fold_Enrichment', 'Raw_p_value', 'FDR')
GO2 <- subset(GO2, OverUnder == "+" & FDR < 0.1)
dim(GO2) # see we have some significant processes to plot
GO2$Treatment <- rep('Sal02', nrow(GO2))

GO <- rbind(GO1, GO2) # merge GO results
# View(GO) # Identify up to 10 processes of interest in the top 30 (based on lowest FDRs) for both the PosFC & negFC lists...since NegFC list only had 9 processes, we will take all of those

terms <- c(GO1$GO_biological_process, 'defense response (GO:0006952)', 'response to biotic stimulus (GO:0009607)',
           'defense response to other organism (GO:0098542)', 'response to bacterium (GO:0009617)', 
           'response to external stimulus (GO:0009605)', 'innate immune response (GO:0045087)',
           'macrophage activation (GO:0042116)', 'response to molecule of bacterial origin (GO:0002237)',
           'immune effector process (GO:0002252)', 'myeloid leukocyte activation (GO:0002274)')
GO <- GO[GO$GO_biological_process %in% terms, ]
GO$GO_biological_process <- factor(GO$GO_biological_process, levels = terms)
ggplot(GO) +
  geom_point(aes(x = GO_biological_process, 
                 y = Treatment,
                 size = Fold_Enrichment,
                 color = FDR)) +
  theme_bw() +
  scale_color_gradient(low = 'red3', high = 'gold') +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))

## Crl07vSal07
GO1 <- read_excel('/Users/jayne.wiarda/Desktop/Experiments/FS27_Salmonella/RNAseq_Bulk_WholeBlood/Results_GO/ResultsGO_Crl07vSal07_NegFC.xlsx')
colnames(GO1) <- c('GO_biological_process', 'Number_BackgroundGenes', 'Number_AnalyzedGenes', 
                   'Expected_AnalyzedGenes', 'OverUnder', 'Fold_Enrichment', 'Raw_p_value', 'FDR')
GO1 <- subset(GO1, OverUnder == "+" & FDR < 0.1)
dim(GO1) # see no significant biological processes to plot

#GO2 <- read_excel('/Users/jayne.wiarda/Desktop/Experiments/FS27_Salmonella/RNAseq_Bulk_WholeBlood/Results_GO/ResultsGO_Crl07vSal07_PosFC.xlsx')
# this gene list didn't exist, so we are all done here

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
#  [1] readxl_1.4.1  ggplot2_3.4.0

#loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.9             locfit_1.5-9.6         lattice_0.20-45        png_0.1-8              Biostrings_2.66.0      assertthat_0.2.1       utf8_1.2.2             cellranger_1.1.0       R6_2.5.1              
#[10] GenomeInfoDb_1.34.4    RSQLite_2.2.19         httr_1.4.4             pillar_1.8.1           zlibbioc_1.44.0        rlang_1.0.6            rstudioapi_0.14        annotate_1.76.0        blob_1.2.3            
#[19] S4Vectors_0.36.1       Matrix_1.5-3           labeling_0.4.2         BiocParallel_1.32.4    geneplotter_1.76.0     RCurl_1.98-1.9         bit_4.0.5              munsell_0.5.0          DelayedArray_0.24.0   
#[28] compiler_4.2.1         pkgconfig_2.0.3        BiocGenerics_0.44.0    tidyselect_1.2.0       KEGGREST_1.38.0        tibble_3.1.8           GenomeInfoDbData_1.2.9 IRanges_2.32.0         codetools_0.2-18      
#[37] matrixStats_0.63.0     XML_3.99-0.13          fansi_1.0.3            crayon_1.5.2           dplyr_1.0.10           withr_2.5.0            bitops_1.0-7           grid_4.2.1             xtable_1.8-4          
#[46] gtable_0.3.1           lifecycle_1.0.3        DBI_1.1.3              magrittr_2.0.3         scales_1.2.1           cli_3.4.1              cachem_1.0.6           farver_2.1.1           XVector_0.38.0        
#[55] vctrs_0.5.1            generics_0.1.3         RColorBrewer_1.1-3     tools_4.2.1            bit64_4.0.5            Biobase_2.58.0         glue_1.6.2             MatrixGenerics_1.10.0  parallel_4.2.1 
#[64] fastmap_1.1.0          AnnotationDbi_1.60.0   colorspace_2.0-3       GenomicRanges_1.50.1   memoise_2.0.1         