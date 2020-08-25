############# DLBC_maftools ##################
# MAF files from Difuse Large B-cell lymphomas

# maftools: Summarize, Analyze and Visualize Mutation Anotated Files (MAF) Files
# URL: https://www.bioconductor.org/packages/release/bioc/html/maftools.html
# version 2.4.05


## Installing packages
packages_bioconductor <- c("TCGAbiolinks","maftools","BSgenome.Hsapiens.UCSC.hg38")
packages_cran <- c("tidyverse", "pheatmap")

#use this function to check if each package is on the local machine
#if a package is installed, it will be loaded
#if any are not, the missing package(s) will be installed from Bioconductor and loaded
package.check <- lapply(packages_bioconductor, FUN = function(x) {
        if (!require(x, character.only = TRUE)) {
                BiocManager::install(x, dependencies = TRUE)
                library(x, character.only = TRUE)
        }
})
package.check <- lapply(packages_cran, FUN = function(x) {
        if (!require(x, character.only = TRUE)) {
                install.packages(x, dependencies = TRUE)
                library(x, character.only = TRUE)
        }
})

rm(packages_cran, packages_bioconductor, package.check)

#setwd()

## Reading DLBC Maf files ---------------------------

# download MAF aligned against hg19
# it saves data inside a GDCdata and project name directory 
maf <- GDCquery_Maf("DLBC", pipelines = "muse", directory = "GDCdata")
sort(colnames(maf))

# MAF object contains main maf file, summarized data and any associated sample annotations
dlbc.maf <- read.maf(maf = maf, useAll = T) 

# checking
getSampleSummary(dlbc.maf) # 37 samples (Tumor_Sample_Barcode) -> not 29??
getGeneSummary(dlbc.maf) # 2632 genes (hugo)
getClinicalData(dlbc.maf) # 37 samples, no clinical data  
getFields(dlbc.maf) # 120 variables 

# writes an output file
write.mafSummary(maf = dlbc.maf, basename = 'dlbc.maf')


## Reading clinical indexed data ------------------

# same as in data portal
clinical <- GDCquery_clinic(project = "TCGA-DLBC", type = "clinical", save.csv = FALSE)
sort(colnames(clinical))

colnames(clinical)[1] <- "Tumor_Sample_Barcode"
# clinical$Overall_Survival_Status <- 1
# clinical$Overall_Survival_Status[which(clinical$vital_status == "Dead")] <- 0
clinical$time <- clinical$days_to_death
clinical$time[is.na(clinical$days_to_death)] <- clinical$days_to_last_follow_up[is.na(clinical$days_to_death)]

# create object for survival analysis 
dlbc.mafclin <- read.maf(maf = maf, clinicalData = clinical, isTCGA = T)


## Vizualizing -----------------------------------

# displays variants in each sample and variant types summarized by Variant_Classification
plotmafSummary(maf = dlbc.maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

# oncoplot for top ten mutated genes (costumize oncoplots!)
oncoplot(maf = dlbc.maf, top = 20)

# transition and transversions
mafDlbc.titv <- titv(maf = dlbc.maf, plot = FALSE, useSyn = TRUE)
plotTiTv(res = mafDlbc.titv)

# lollipop plot for IGHV2-70, which is one of the most frequent mutated gene in DLBC
lollipopPlot(maf = dlbc.maf, gene = 'CARD11', AACol = 'Amino_acids', showMutationRate = TRUE)
lollipopPlot(maf = dlbc.mafclin, gene = 'MUC16', AACol = 'Amino_acids', showDomainLabel = FALSE)

# rainfall plots highlights hyper-mutated genomic regions
# Kataegis: genomic segments containing six or more consecutive mutations with an average inter-mutation distance of less than or equal to 1,00 bp
rainfallPlot(maf = dlbc.maf, detectChangePoints = TRUE, pointSize = 0.6)

# mutation load
dlbcl.mutload <- tcgaCompare(maf = dlbc.maf)

# plots Variant Allele Frequencies
# clonal genes usually have mean allele frequency around ~50% assuming pure sample
# it looks for column t_vaf containing vaf information, if the field name is different from t_vaf, we can manually specify it using argument vafCol
# plotVaf(maf = luad.maf, vafCol = 'i_TumorVAF_WU')

## Comparing two cohorts (MAFs)

# Considering only genes mutated in at-least in 5 samples in one of the cohort to avoid bias
# lusc.vs.luad <- mafCompare(m1 = lusc.maf, m2 = luad.maf, m1Name = 'Lusc', m2Name = 'Luad', minMut = 5)
# forestPlot(mafCompareRes = lusc.vs.luad, pVal = 0.1, color = c('royalblue', 'maroon'), geneFontSize = 0.8)
# genes = c("TTTN", "TP53", "MIC16", "RYR2", "CSMD5")
# coOncoplot(m1 = lusc.maf, m2 = luad.maf, m1Name = 'Lusc', m2Name = 'Luad', genes = genes, removeNonMutated = TRUE)
# coBarplot(m1 = lusc.maf, m2 = luad.maf, m1Name = 'Lusc', m2Name = 'Luad')
# lollipopPlot2(m1 = lusc.maf, m2 = luad.maf, gene = 'TP53', AACol1 = 'AA_MAF', AACol2 = 'AA_MAF', showMutationRate = TRUE)


## Analysis -------------------------------------

# exclusive/co-occurance event analysis on top 10 mutated genes (pair-wise Fisher’s Exact test)
somaticInteractions(maf = dlbc.maf, top = 10, pvalue = c(0.05, 0.1))

# detecting cancer driver genes based on positional clustering
# most of the variants in cancer causing genes are enriched at few specific loci (aka hot-spots)
dlbc.sig <- oncodrive(maf = dlbc.maf, AACol = 'AA_MAF', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = dlbcl.sig, fdrCutOff = 0.1, useFraction = TRUE) # AACol = 'Amino_acids'?

# pfam domain
dlbc.pfam = pfamDomains(maf = dlbc.maf, AACol = 'AA_MAF', top = 10) # idem


## Reading gistic or CNV -------------------------

all.lesions <- system.file("extdata", "all_lesions.conf_99.txt", package = "maftools")
amp.genes <- system.file("extdata", "amp_genes.conf_99.txt", package = "maftools")
del.genes <- system.file("extdata", "del_genes.conf_99.txt", package = "maftools")
scores.gis <- system.file("extdata", "scores.gistic", package = "maftools")

dlbc.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes,
                         gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = TRUE)

# GISTIC object
dlbc.gistic

# checking 
getSampleSummary(dlbc.gistic) # 191 samples (Tumor_Sample_Barcode) -> not 29??
getGeneSummary(dlbc.gistic) # 2632 genes (hugo)
getCytobandSummary(dlbc.gistic) # 16 cytobands

write.GisticSummary(gistic = dlbc.gistic, basename = 'dlbc.gistic')


## Vizualizing Gistic ---------------------------------

# genome plot
gisticChromPlot(gistic = dlbc.gistic, markBands = "all")

# bubble plot
gisticBubblePlot(gistic = dlbc.gistic)

# oncoplot 
gisticOncoPlot(gistic = dlbc.gistic, clinicalData = getClinicalData(x = dlbc.mafclin), clinicalFeatures = 'ann_arbor_clinical_stage', sortByAnnotation = FALSE, top = 10) # error 

gisticOncoPlot(gistic = dlbc.gistic, top = 20)


## Survival Analysis ------------------------------

# it requires input data with Tumor_Sample_Barcode, binary event (1/0) and time to event.
mafSurvival(maf = dlbc.mafclin, genes = 'IGHV2-70', time = 'time', Status = 'vital_status', isTCGA = TRUE)

# identify a set of genes (of size 2) to predict poor prognostic groups
prog_geneset <- survGroup(maf = dlbc.mafclin, top = 20, geneSetSize = 2, time = "days_to_last_follow_up", Status = "vital_status", verbose = FALSE)


## Clinical enrichment analysis ------------------

aacs.enrich <- clinicalEnrichment(maf = dlbc.mafclin, clinicalFeature = 'ann_arbor_clinical_stage')

# identify mutations associated with clinicalFeature
# results are returned as a list. Significant associations p-value < 0.05
aacs.enrich$groupwise_comparision[p_value < 0.05] # it takes too long!!

plotEnrichmentResults(enrich_res = aacs.enrich, pVal = 0.05)


## Oncogenic Signaling Pathways --------------------

OncogenicPathways(maf = dlbc.maf)

PlotOncogenicPathways(maf = dlbc.maf, pathways = "NOTCH")
# tumor suppressor genes (red) and oncogenes (blue)


## Tumor heterogeneity and MATH scores -------------
# heterogeneity can be inferred by clustering variant allele frequencies
# < code here >


## Mutational signatures --------------------------
# signatures can be extracted by decomposing matrix of nucleotide substitutions and compared to validated signatures
# < code here >


-----------------------------------------
-----------------------------------------
     
## References
     
# citation("maftools")
# Mayakonda A, Lin DC, Assenov Y, Plass C, Koeffler HP. 2018. Maftools: efficient and comprehensive analysis of somatic variants in cancer. Genome Resarch PMID: 30341162
     
# about MAF files: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
# about download data: https://www.bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/query.html
# about costumize oncoplots: http://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/oncoplots.html
# about validated signatures: https://cancer.sanger.ac.uk/cosmic/signatures
     
# sessionInfo()
