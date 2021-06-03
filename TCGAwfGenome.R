#### TCGAWorkflow - Genome Analysis #####

# This workflow provides a series of biologically focused integrative analyses of different molecular data.
# Bioconductor version: Release (3.12)
# Vignette: https://www.bioconductor.org/packages/release/workflows/vignettes/TCGAWorkflow/inst/doc/TCGAWorkflow.html
# Publication: https://f1000research.com/articles/5-1542/v2.


### 1. Install Packages ------------------------------
if (!"BiocManager" %in% rownames(installed.packages()))
     install.packages("BiocManager")
BiocManager::install("TCGAWorkflow")
BiocManager::install("TCGAWorkflowData")

library(TCGAWorkflow)
library(TCGAWorkflowData) # subset of LGG and GMB TCGA data
library(DT) # to visualize the results
library(tidyverse)

## 2. Access and Downloda Data ----------------------------
# Data from GDC data portal and GDC Legacy Archive can be accessed by TCGAbiolinks
# Data from Firehose can be accessed by RTCGAToolbox

library(TCGAbiolinks)
# Obs: The data in the legacy database has been aligned to hg19
 query.met.gbm <- GDCquery(
   project = "TCGA-GBM",
   legacy = TRUE,
   data.category = "DNA methylation",
   platform = "Illumina Human Methylation 450",
   barcode = c("TCGA-76-4926-01B-01D-1481-05", "TCGA-28-5211-01C-11D-1844-05")
 )
 GDCdownload(query.met.gbm)
 met.gbm.450 <- GDCprepare(
   query = query.met.gbm,
   save = TRUE,
   save.filename = "gbmDNAmet450k.rda",
   summarizedExperiment = TRUE
 )

 query.met.lgg <- GDCquery(
   project = "TCGA-LGG",
   legacy = TRUE,
   data.category = "DNA methylation",
   platform = "Illumina Human Methylation 450",
   barcode = c("TCGA-HT-7879-01A-11D-2399-05", "TCGA-HT-8113-01A-11D-2399-05")
 )
 GDCdownload(query.met.lgg)
 met.lgg.450 <- GDCprepare(query = query.met.lgg,
                           save = TRUE,
                           save.filename = "lggDNAmet450k.rda",
                           summarizedExperiment = TRUE)
 met.gbm.lgg <- SummarizedExperiment::cbind(met.lgg.450, met.gbm.450)

 query.exp.lgg <- GDCquery(
   project = "TCGA-LGG",
   legacy = TRUE,
   data.category = "Gene expression",
   data.type = "Gene expression quantification",
   platform = "Illumina HiSeq",
   file.type = "results",
   sample.type = "Primary solid Tumor"
 )
 GDCdownload(query.exp.lgg)
 exp.lgg <- GDCprepare(query = query.exp.lgg, save = TRUE, save.filename = "lggExp.rda")

 query.exp.gbm <- GDCquery(
   project = "TCGA-GBM",
   legacy = TRUE,
   data.category = "Gene expression",
   data.type = "Gene expression quantification",
   platform = "Illumina HiSeq",
   file.type = "results",
   sample.type = "Primary solid Tumor"
 )
 GDCdownload(query.exp.gbm)
 exp.gbm <- GDCprepare(query = query.exp.gbm, save = TRUE, save.filename = "gbmExp.rda")
 exp.gbm.lgg <- SummarizedExperiment::cbind(exp.lgg, exp.gbm)


# Copy number variation aligned to hg38 ----------------------
query.cnv.acc <- GDCquery(
   project = "TCGA-ACC",
   data.category = "Copy Number Variation",
   data.type = "Copy Number Segment",
   barcode = c("TCGA-OR-A5KU-01A-11D-A29H-01", "TCGA-OR-A5JK-01A-11D-A29H-01")
 )
GDCdownload(query.cnv.acc)
data <- GDCprepare(query.cnv.acc)
cnv.acc <- GDCprepare(query = query.cnv.acc, save = TRUE, save.filename = "accCNV.rda")

query.cnvmask.acc <- GDCquery(
   project = "TCGA-ACC",
   data.category = "Copy Number Variation",
   data.type = "Masked Copy Number Segment",
   sample.type = c("Primary solid Tumor")
 ) # see the barcodes with getResults(query)$cases
GDCdownload(query.cnvmask.acc)
cnvMask.acc <- GDCprepare(query = query.cnv.acc, save = TRUE, save.filename = "accCNVmask.rda")


## Summarized Expreriment -------------------------------------------------
# This object stores rectangular matrices of experimental results and metadata, produced by sequencing and microarray experiments.
library(SummarizedExperiment)

# Load object from TCGAWorkflowData package
# this object will be created in the further sections
data(GBMIllumina_HiSeq) 

# get expression matrix
data <- assay(gbm.exp)
datatable(
     data = data[1:10,], 
     options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
     rownames = TRUE
)

# get genes information
genes.info <- rowRanges(gbm.exp)
genes.info

# get sample information
sample.info <- colData(gbm.exp)
datatable(
     data = as.data.frame(sample.info), 
     options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
     rownames = FALSE
)


### Clinical Data -------------------------------------------------

## 1.Download only the indexed GDC clinical data 
gbm_clin <- GDCquery_clinic(project = "TCGA-GBM", type = "Clinical")
lgg_clin <- GDCquery_clinic(project = "TCGA-LGG", type = "Clinical")

# Bind the results, as the columns might not be the same, apply
# plyr rbind.fill, to have all columns from both files
clinical <- plyr::rbind.fill(gbm_clin,lgg_clin)
datatable(clinical[1:10,], options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

## 2.Fetch clinical data directly from the clinical XML files
# if barcode is not set, it will consider all samples
query.clinical <- GDCquery(
     project = "TCGA-GBM",
     file.type = "xml",
     data.category = "Clinical",
     barcode = c("TCGA-08-0516","TCGA-02-0317")
) 
GDCdownload(query.clinical)

clinic.patient <- GDCprepare_clinic(query.clinical, clinical.info = "patient")
datatable(clinic.patient, options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

clinic.drug <- GDCprepare_clinic(query.clinical, clinical.info = "drug")
datatable(clinic.drug, options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

clinic.radiation <- GDCprepare_clinic(query.clinical, clinical.info = "radiation")
datatable(clinic.radiation, options = list(scrollX = TRUE,  keys = TRUE), rownames = FALSE)

clinic.admin <- GDCprepare_clinic(query.clinical, clinical.info = "admin")
datatable(clinic.admin, options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)


### Mutation Annotation Format (MAF) --------------------------
# Somatic (or Public) MAF files, are derived from the GDC annotated VCF files, which are further processed to remove low quality and potential germline variants.

maf.lgg <- GDCquery_Maf(tumor = "LGG", pipelines = "mutect2")
data(mafMutect2LGGGBM) 
datatable(maf.lgg[1:10,], options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

subtypes.gmb <- TCGAquery_subtype(tumor = "gbm")
datatable(subtypes.gmb[1:10,], options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)


## Downloading data from Broad TCGA GDAC -------------
# RTCGAToolbox provides access to Firehose Level 3 and 4 data through the function getFirehoseData. 
library(RTCGAToolbox)

# Get the last run dates
lastRunDate <- getFirehoseRunningDates()[1]

# get DNA methylation data, RNAseq2 and clinical data for GBM
firehose.gbm <- getFirehoseData(
   dataset = "GBM",
   runDate = lastRunDate,
   gistic2Date = getFirehoseAnalyzeDates(1),
   Methylation = FALSE,
   clinical = TRUE,
   RNASeq2GeneNorm  = FALSE,
   Mutation = TRUE,
   fileSizeLimit = 10000
 )

firehose.mut <- getData(firehose.gbm,"Mutation")
firehose.clin <- getData(firehose.gbm,"clinical")

# Download GISTIC results
# GISTIC identify genes targeted by somatic copy-number alterations (SCNAs)
lastanalyzedate <- getFirehoseAnalyzeDates(1)
gistic.gmb <- getFirehoseData("GBM",GISTIC = TRUE, gistic2Date = lastanalyzedate)

# get GISTIC results
gistic.allbygene <- getData(gistic.gmb, type = "GISTIC", platform = "AllByGene")
gistic.thresholedbygene <- getData(gistic.gmb, type = "GISTIC", platform = "ThresholdedByGene")

data(GBMGistic)
gistic.allbygene %>% head() %>% gt::gt() # no package called 'dt'
gistic.thresholedbygene %>% head() %>% gt::gt()


### 3. Genomic Analysis ----------------------------

## CNV data pre-processing ----------------------
# CNVs are genomic regions greater than 1 kb with an alteration of copy number between two conditions (e.g., Tumor versus Normal).
query.nocnv.gmb <- GDCquery(
   project = "TCGA-GBM",
   data.category = "Copy number variation",
   legacy = TRUE,
   file.type = "nocnv_hg19.seg",
   sample.type = c("Primary solid Tumor")
 )

# to reduce time we will select only 20 samples
query.gbm.nocnv$results[[1]] <- query.gbm.nocnv$results[[1]][1:20,]
GDCdownload(query.nocnv.gmb, files.per.chunk = 100)
gbm.nocnv <- GDCprepare(query.nocnv.gmb, save = TRUE, save.filename = "GBMnocnvhg19.rda")

## CNV analysis from Broad Institute ------
# For hg38 analysis: https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files
# File: SNP6 GRCh38 Liftover Probeset File for CNV Analysis

# Retrieve probes meta file from broadinstitute website for hg19
gdac.root <- "ftp://ftp.broadinstitute.org/pub/GISTIC2.0/hg19_support/"
file <- paste0(gdac.root, "genome.info.6.0_hg19.na31_minus_frequent_nan_probes_sorted_2.1.txt")
if(!file.exists(basename(file))) downloader::download(file, basename(file))
markersMatrix <-  readr::read_tsv(basename(file), col_names = FALSE, col_types = "ccn", progress = FALSE)
save(markersMatrix, file = "markersMatrix.rda", compress = "xz")

## CNV identification (GAIA) ---------------
#  GAIA is based on a conservative permutation test allowing the estimation of the probability distribution of the contemporary mutations expected for non-driver markers.
cancer <- "GBM"
message(paste0("Starting ", cancer))

# get objects created above
data(GBMnocnvhg19)
data(markersMatrix)

head(cnvMatrix)
head(markersMatrix)

# Add label (0 for loss, 1 for gain)
cnvMatrix <- cnvMatrix %>% dplyr::filter(abs(cnvMatrix$Segment_Mean) > 0.3)
cnvMatrix$Aberration <- ifelse(cnvMatrix$Segment_Mean > 0.3, 1, 0)
# Remove "Segment_Mean" and change col.names
cnvMatrix$Segment_Mean <- NULL
colnames(cnvMatrix) <- c("Sample.Name", "Chromosome", "Start", "End", "Num.of.Markers", "Aberration")

# Substitute Chromosomes "X" and "Y" with "23" and "24"
cnvMatrix$Chromosome[cnvMatrix$Chromosome == "X"] <- "23"
cnvMatrix$Chromosome[cnvMatrix$Chromosome == "Y"] <- "24"
cnvMatrix$Chromosome <- as.integer(cnvMatrix$Chromosome)
cnvMatrix <- cnvMatrix %>% as.data.frame()

# Recurrent CNV identification with GAIA
colnames(markersMatrix) <- c("Probe.Name", "Chromosome", "Start")
unique(markersMatrix$Chromosome)

markersMatrix$Chromosome[markersMatrix$Chromosome == "X"] <- "23"
markersMatrix$Chromosome[markersMatrix$Chromosome == "Y"] <- "24"
markersMatrix$Chromosome <- as.integer(markersMatrix$Chromosome)
markersMatrix <- markersMatrix %>% as.data.frame()

markerID <- paste(markersMatrix$Chromosome,markersMatrix$Start, sep = ":")
# Removed duplicates
markersMatrix <- markersMatrix[!duplicated(markerID),]
# Filter markersMatrix for common CNV
markerID <- paste(markersMatrix$Chromosome,markersMatrix$Start, sep = ":")

file <- "ftp://ftp.broadinstitute.org/pub/GISTIC2.0/hg19_support/CNV.hg19.bypos.111213.txt"
if(!file.exists(basename(file))) downloader::download(file, basename(file))
commonCNV <- readr::read_tsv(basename(file), progress = FALSE)
#data(CNV.hg19.bypos.111213)
commonID <- paste(commonCNV$Chromosome,commonCNV$Start, sep = ":")
markersMatrix_fil <- markersMatrix[!markerID %in% commonID,]

library(gaia)
set.seed(200)
markers_obj <- load_markers(
     marker_matrix = markersMatrix_fil %>% as.data.frame()
)

nbsamples <- length(unique(cnvMatrix$Sample.Name))
cnv_obj <- load_cnv(
     segmentation_matrix = cnvMatrix  %>% as.data.frame(), 
     markers_list = markers_obj, 
     num_of_samples = nbsamples
)

suppressWarnings({
     results <- runGAIA(
          cnv_obj = cnv_obj,
          markers_obj = markers_obj,
          output_file_name = paste0("GAIA_",cancer,"_flt.txt"),
          aberrations = -1,  # -1 to all aberrations
          chromosomes = 9, # -1 to all chromosomes
          approximation = TRUE, # Set to TRUE to speed up the time requirements
          num_iterations = 5000, # Reduced to speed up the time requirements
          threshold = 0.25)
})
# Set q-value threshold
# Use a smalled value for your analysis. We set this as high values
# due to the small number of samples which did not reproduced
# results with smaller q-values
threshold <- 0.3

# Plot the results
mode(results) <- "numeric" # transform character to numeric
RecCNV <- results %>% as.data.frame()
RecCNV$score <- 0

# replace 0 q-value to min non-zero
minval <- format(min(RecCNV[RecCNV[,"q-value"] != 0,"q-value"]), scientific = FALSE)
minval <- substring(minval,1, nchar(minval) - 1) %>% as.numeric
RecCNV$`q-value`[RecCNV$`q-value` == 0] <- minval

RecCNV$score <- -log10(RecCNV$`q-value`)
RecCNV[RecCNV[,"q-value"] == as.numeric(minval),]

gaiaCNVplot(RecCNV,threshold)
save(results, RecCNV, threshold, file = paste0(cancer,"_CNV_results.rda"))

## Recurrent CNV gene annotation (GenomicRanges) -----------------
library(GenomicRanges)

# Get gene information from GENCODE using biomart
genes <- TCGAbiolinks:::get.GRCh.bioMart(genome = "hg19")
genes <- genes[genes$external_gene_name != "" & genes$chromosome_name %in% c(1:22,"X","Y"),]
genes[genes$chromosome_name == "X", "chromosome_name"] <- 23
genes[genes$chromosome_name == "Y", "chromosome_name"] <- 24
genes$chromosome_name <- sapply(genes$chromosome_name,as.integer)
genes <- genes[order(genes$start_position),]
genes <- genes[order(genes$chromosome_name),]
genes <- genes[,c("external_gene_name", "chromosome_name", "start_position","end_position")]
colnames(genes) <- c("GeneSymbol","Chr","Start","End")
genes_GR <- makeGRangesFromDataFrame(genes,keep.extra.columns = TRUE)
save(genes_GR,genes,file = "genes_GR.rda", compress = "xz")

# Get gene information from GENCODE using biomart
data(genes_GR) # downloaded in the previous step (available in TCGAWorkflowData)

load(paste0(cancer,"_CNV_results.rda"))
sCNV <- RecCNV[RecCNV[,"q-value"] <= threshold,c(1:4,6)]
sCNV <- sCNV[order(sCNV$Chromosome,sCNV$`Region Start [bp]`),]
colnames(sCNV) <- c("Chr","Aberration","Start","End","q-value")
sCNV_GR <- makeGRangesFromDataFrame(sCNV,keep.extra.columns = TRUE)

hits <- findOverlaps(genes_GR, sCNV_GR, type = "within")
sCNV_ann <- cbind(sCNV[subjectHits(hits),],genes[queryHits(hits),])
AberrantRegion <- paste0(sCNV_ann[,1],":",sCNV_ann[,3],"-",sCNV_ann[,4])
GeneRegion <- paste0(sCNV_ann[,7],":",sCNV_ann[,8],"-",sCNV_ann[,9])
AmpDel_genes <- cbind(sCNV_ann[,c(6,2,5)],AberrantRegion,GeneRegion)
AmpDel_genes[AmpDel_genes[,2] == 0,2] <- "Del"
AmpDel_genes[AmpDel_genes[,2] == 1,2] <- "Amp"
rownames(AmpDel_genes) <- NULL

save(RecCNV, AmpDel_genes, file = paste0(cancer,"_CNV_results.rda"))

knitr::kable(head(AmpDel_genes), caption = "Chromosome 9 recurrent deleted genes in LGG")

## Multiple Genomic Alteration Events (MAFTools) -------
library(maftools)

# Download Mutation Annotation Format (MAF) files
LGGmut <- GDCquery_Maf(tumor = "LGG", pipelines = "mutect2")
GBMmut <- GDCquery_Maf(tumor = "GBM", pipelines = "mutect2")

# Merge them
mut <- plyr::rbind.fill(LGGmut, GBMmut)
save(mut,file ="mafMutect2LGGGBM.rda", compress = "xz") # change 1st parameter: maf to mut

# recovering data from TCGAWorkflowData package.
data(mafMutect2LGGGBM)

# To prepare for maftools we will also include clinical data
# For a mutant vs WT survival analysis 
# get indexed clinical patient data for GBM samples
gbm_clin <- GDCquery_clinic(project = "TCGA-GBM", type = "Clinical")
# get indexed clinical patient data for LGG samples
lgg_clin <- GDCquery_clinic(project = "TCGA-LGG", type = "Clinical")
# Bind the results, as the columns might not be the same,
# we will will plyr rbind.fill, to have all columns from both files
clinical <- plyr::rbind.fill(gbm_clin,lgg_clin)
colnames(clinical)[1] <- "Tumor_Sample_Barcode"

# we need to create a binary variable 1 is dead 0 is not dead
plyr::count(clinical$vital_status)
clinical$Overall_Survival_Status <- 1 # dead
clinical$Overall_Survival_Status[which(clinical$vital_status != "Dead")] <- 0

# If patient is not dead we don't have days_to_death (NA)
# we will set it as the last day we know the patient is still alive
clinical$time <- clinical$days_to_death
clinical$time[is.na(clinical$days_to_death)] <- clinical$days_to_last_follow_up[is.na(clinical$days_to_death)]

# Create object to use in maftools
maf <- read.maf(maf = mut, clinicalData = clinical, isTCGA = TRUE)

# plotmafSummary 
plotmafSummary(maf = maf,rmOutlier = TRUE,addStat = 'median',dashboard = TRUE)

# oncoplot 
oncoplot(maf = maf,top = 20,legendFontSize = 8,clinicalFeatures = c("tissue_or_organ_of_origin") )

# mafSurvival
plot <- mafSurvival(maf = maf,genes = "TP53",time = 'time',Status = 'Overall_Survival_Status', isTCGA = TRUE)

# LGGmut <- GDCquery_Maf(tumor = "LGG", pipelines = "mutect2")


### Genomic aberration overview - Circos plot ----------------- 

# Retrieve curated mutations for selected cancer (e.g. "LGG") 
data(mafMutect2LGGGBM)
# Select only potentially damaging mutations
LGGmut %>% 
     dplyr::filter(Variant_Classification %in% 
                        c(  "Missense_Mutation",
                             "Nonsense_Mutation",
                             "Nonstop_Mutation",
                             "Frame_Shift_Del",
                             "Frame_Shift_Ins" )
     )
# Select recurrent mutations (identified in at least two samples)
mut.id <- paste0(LGGmut$Chromosome,":",LGGmut$Start_Position,"-",LGGmut$End_Position,"|",LGGmut$Reference_Allele)
mut <- cbind(mut.id, LGGmut)
# Prepare selected mutations data for circos plot
s.mut <- mut[mut$mut.id %in% unique(mut.id[duplicated(mut.id)]),]
s.mut <- s.mut[,c("Chromosome","Start_Position","End_Position","Variant_Classification","Hugo_Symbol")]
s.mut <- unique(s.mut)
typeNames <- unique(s.mut[,4])
type <- c(4:1)
names(type) <- typeNames[1:4]
Type <- type[s.mut[,4]]
s.mut <- cbind(s.mut,Type)
s.mut <- s.mut[,c(1:3,6,4,5)]

# Load recurrent CNV data for selected cancer (e.g. "LGG")
load("GBM_CNV_results.rda")
# Prepare selected sample CNV data for circos plot
s.cnv <- as.data.frame(RecCNV[RecCNV[,"q-value"] <= threshold,c(1:4,6)])
s.cnv <- s.cnv[,c(1,3,4,2)]
s.cnv[s.cnv$Chromosome == 23,"Chromosome"] <- "X"
s.cnv[s.cnv$Chromosome == 24,"Chromosome"] <- "Y"
Chromosome <- paste0("chr",s.cnv[,1])
s.cnv <- cbind(Chromosome, s.cnv[,-1])
s.cnv[,1] <- as.character(s.cnv[,1])
s.cnv[,4] <- as.character(s.cnv[,4])
s.cnv <- cbind(s.cnv,CNV=1)
colnames(s.cnv) <- c("Chromosome","Start_position","End_position","Aberration_Kind","CNV")

library(circlize)
# Draw genomic circos plot
par(mar=c(1,1,1,1), cex=1)
circos.initializeWithIdeogram()
# Add CNV results
colors <- c("forestgreen","firebrick")
names(colors)  <- c(0,1)
circos.genomicTrackPlotRegion(
     s.cnv,  ylim = c(0,1.2),
     panel.fun = function(region, value, ...) {
          circos.genomicRect(region, value, ytop.column = 2, ybottom = 0,
                             col = colors[value[[1]]], 
                             border="white")
          cell.xlim = get.cell.meta.data("cell.xlim")
          circos.lines(cell.xlim, c(0, 0), lty = 2, col = "#00000040")
     })
# Add mutation results
colors <- c("blue","green","red","gold")
names(colors)  <- typeNames[1:4]
circos.genomicTrackPlotRegion(
     s.mut, ylim = c(1.2,4.2),
     panel.fun = function(region, value, ...) {
          circos.genomicPoints(region, value, cex = 0.8, pch = 16, col = colors[value[[2]]], ...)
     })

circos.clear()

legend(
     -0.2,
     0.2,
     bty = "n",
     y.intersp = 1,
     c("Amp", "Del"),
     pch = 15,
     col = c("firebrick", "forestgreen"),
     title = "CNVs",
     text.font = 1,
     cex = 0.4,
     title.adj = 0
)
legend(
     -0.2,
     0,
     bty = "n",
     y.intersp = 1,
     names(colors),
     pch = 16,
     col = colors,
     title = "Mutations",
     text.font = 1,
     cex = 0.4,
     title.adj = 0
)

# Zoom in one chromossome
par(mar=c(1,1,1,1),cex=1.5)
circos.par(
     "start.degree" = 90, 
     canvas.xlim = c(0, 1), 
     canvas.ylim = c(0, 1), 
     gap.degree = 270,
     cell.padding = c(0, 0, 0, 0), 
     track.margin = c(0.005, 0.005)
)
circos.initializeWithIdeogram(chromosome.index = "chr17")
circos.par(cell.padding = c(0, 0, 0, 0))
# Add CNV results
colors <- c("forestgreen","firebrick")
names(colors)  <- c(0,1)
circos.genomicTrackPlotRegion(
     data = s.cnv,  ylim = c(0,1.2),
     panel.fun = function(region, value, ...) {
          circos.genomicRect(region, value, ytop.column = 2, ybottom = 0,
                             col = colors[value[[1]]], 
                             border="white")
          cell.xlim = get.cell.meta.data("cell.xlim")
          circos.lines(cell.xlim, c(0, 0), lty = 2, col = "#00000040")
     }
)

# Add mutation results representing single genes
genes.mut <- paste0(s.mut$Hugo_Symbol,"-",s.mut$Type)
s.mutt <- cbind(s.mut,genes.mut)
n.mut <- table(genes.mut)
idx <- !duplicated(s.mutt$genes.mut)
s.mutt <- s.mutt[idx,]
s.mutt <- cbind(s.mutt,num=n.mut[s.mutt$genes.mut])
genes.num <- paste0(s.mutt$Hugo_Symbol," (",s.mutt$num.Freq,")")
s.mutt <- cbind(s.mutt[,-c(6:8)],genes.num)
s.mutt[,6] <- as.character(s.mutt[,6])
s.mutt[,4] <- s.mutt[,4]/2
s.mutt$num.Freq <- NULL
colors <- c("blue","green","red","gold")
names(colors)  <- typeNames[1:4]
circos.genomicTrackPlotRegion(
     data = s.mutt, ylim = c(0.3,2.2), track.height = 0.05,
     panel.fun = function(region, value, ...) {
          circos.genomicPoints(region, value, cex = 0.4, pch = 16, col = colors[value[[2]]], ...)
     }
)

circos.genomicTrackPlotRegion(s.mutt, ylim = c(0, 1), track.height = 0.1, bg.border = NA)
i_track = get.cell.meta.data("track.index")

circos.genomicTrackPlotRegion(
     data = s.mutt, ylim = c(0,1),
     panel.fun = function(region, value, ...) {
          circos.genomicText(
               region, value, 
               y = 1, 
               labels.column = 3,
               col = colors[value[[2]]],
               facing = "clockwise", adj = c(1, 0.5),
               posTransform = posTransform.text, cex = 0.4, niceFacing = TRUE)
     }, track.height = 0.1, bg.border = NA
)

circos.genomicPosTransformLines(
     data = s.mutt,
     posTransform = function(region, value)
          posTransform.text(region, 
                            y = 1, 
                            labels = value[[3]],
                            cex = 0.4, track.index = i_track+1),
     direction = "inside", track.index = i_track
)

circos.clear()

legend(
     x = 0.25,
     y =  0.2, 
     bty = "n", 
     y.intersp = 1, 
     c("Amp","Del"),
     pch = 15, 
     col = c("firebrick","forestgreen"),
     title = "CNVs", 
     text.font = 1, 
     cex = 0.4, 
     title.adj = 0
)
legend(
     x = 0, 
     y = 0.2, 
     bty = "n", 
     y.intersp = 1, 
     names(colors), 
     pch = 16, 
     col = colors, 
     title = "Mutations", 
     text.font = 1, cex = 0.4, 
     title.adj = 0
)


### 4. DNA Methylation Analysis  ----------------------------
# Samples from TCGAbiolinks
lgg.samples <- matchedMetExp("TCGA-LGG", n = 10)
gbm.samples <- matchedMetExp("TCGA-GBM", n = 10)
samples <- c(lgg.samples,gbm.samples)

# For methylation it is quicker in this case to download the tar.gz file
# and get the samples we want instead of downloading files by files
query <- GDCquery(
   project = c("TCGA-LGG","TCGA-GBM"),
   data.category = "DNA methylation",
   platform = "Illumina Human Methylation 450",
   legacy = TRUE,
   barcode = samples
 )
GDCdownload(query)
met <- GDCprepare(query, save = FALSE)

# We will use only chr9 to make the example faster
met <- subset(met,subset = as.character(seqnames(met)) %in% c("chr9"))
# This data is avaliable in the package (object elmerExample)
data(elmerExample)

# Mean methylation ----------------------------
# Plot a barplot for the groups in the disease column in the
# summarizedExperiment object
install.packages('SummarizedExperiment')
library(SummarizedExperiment)

# remove probes with NA (similar to na.omit)
met <- met[rowSums(is.na(assay(met))) == 0,]

df <- data.frame(
     "Sample.mean" = colMeans(assay(met), na.rm = TRUE),
     "groups" = met$project_id
)

library(ggpubr)
ggpubr::ggboxplot(
     data = df,
     y = "Sample.mean",
     x = "groups",
     color = "groups",
     add = "jitter",
     ylab = expression(paste("Mean DNA methylation (", beta, "-values)")),
     xlab = ""
) + stat_compare_means() 


# Differentially methylated CpG sites ---------- 
# change TCGAanalyze_DMC to TCGAanalyze_DMR function
# https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/414
dmc <- TCGAanalyze_DMR(
     data = met,
     groupCol = "project_id", # a column in the colData matrix
     group1 = "TCGA-GBM", # a type of the disease type column
     group2 = "TCGA-LGG", # a type of the disease column
     p.cut = 0.05,
     diffmean.cut = 0.15,
     save = FALSE,
     legend = "State",
     plot.filename = "LGG_GBM_metvolcano.png",
     cores = 1 # if set to 1 there will be a progress bar
)


# DNA Methylation heatmap -------------------------
library(ComplexHeatmap)
clinical <- plyr::rbind.fill(gbm_clin,lgg_clin)

# get the probes that are Hypermethylated or Hypomethylated
# met is the same object of the section 'DNA methylation analysis'
status.col <- "status"
probes <- rownames(dmc)[grep("hypo|hyper",dmc$status,ignore.case = TRUE)]
sig.met <- met[probes,]


# top annotation, which sampples are LGG and GBM
# We will add clinical data as annotation of the samples
# we will sort the clinical data to have the same order of the DNA methylation matrix
clinical.order <- clinical[match(substr(colnames(sig.met),1,12),clinical$bcr_patient_barcode),]

ta <- HeatmapAnnotation(
     df = clinical.order[, c("disease", "gender", "vital_status", "race")],
     col = list(
          disease = c("LGG" = "grey", "GBM" = "black"),
          gender = c("male" = "blue", "female" = "pink")
     )
)

# row annotation: add the status for LGG in relation to GBM
# For exmaple: status.gbm.lgg Hypomethyated means that the
# mean DNA methylation of probes for lgg are hypomethylated
# compared to GBM ones.
ra = rowAnnotation(
     df = dmc[probes, status.col],
     col = list(
          "status.TCGA.GBM.TCGA.LGG" =
               c("Hypomethylated" = "orange",
                 "Hypermethylated" = "darkgreen")
     ),
     width = unit(1, "cm")
)

heatmap  <- Heatmap(
     matrix = assay(sig.met),
     name = "DNA methylation",
     col = matlab::jet.colors(200),
     show_row_names = FALSE,
     cluster_rows = TRUE,
     cluster_columns = FALSE,
     show_column_names = FALSE,
     bottom_annotation = ta,
     column_title = "DNA Methylation"
) 
# Save to pdf
png("heatmap.png",width = 600, height = 400)
draw(heatmap, annotation_legend_side =  "bottom")
dev.off()

## 5. Motif Analysis -----------------
library(rGADEM)
library(BSgenome.Hsapiens.UCSC.hg19)
library(motifStack)
library(SummarizedExperiment)
library(dplyr)

probes <- rowRanges(met)[rownames(dmc)[grep("hypo|hyper",dmc$status,ignore.case = TRUE)],]

# Get hypo/hyper methylated probes and make a 200bp window 
# surrounding each probe.
install.packages('GenomicRanges')
library(GenomicRanges)

# ERROR FROM HERE
sequence <- GRanges(
     seqnames = as.character(seqnames(probes)),
     IRanges(
          start = ranges(probes) %>% as.data.frame() %>% dplyr::pull("start") - 100,
          end = ranges(probes) %>% as.data.frame() %>% dplyr::pull("end") + 100), 
     strand = "*"
) 

#look for motifs
gadem <- GADEM(sequence, verbose = FALSE, genome = Hsapiens)

# How many motifs were found?
nMotifs(gadem)

# get the number of occurrences
nOccurrences(gadem)

# view all sequences consensus
consensus(gadem)

# Print motif
pwm <- getPWM(gadem)
pfm  <- new("pfm",mat = pwm[[1]],name = "Novel Site 1")
plotMotifLogo(pfm)

# Number of instances of motif 1?
length(gadem@motifList[[1]]@alignList)


tabl <- "
|                  Histone marks                 |                                                   Role                                                  |
|:----------------------------------------------:|:-------------------------------------------------------------------------------------------------------:|
| Histone H3 lysine 4 trimethylation (H3K4me3)   | Promoter regions [@heintzman2007distinct,@bernstein2005genomic]                                      |
| Histone H3 lysine 4 monomethylation (H3K4me1)  | Enhancer regions [@heintzman2007distinct]                                                           |
| Histone H3 lysine 36 trimethylation (H3K36me3) | Transcribed regions                                                                                     |
| Histone H3 lysine 27 trimethylation (H3K27me3) | Polycomb repression [@bonasio2010molecular]                                                         |
| Histone H3 lysine 9 trimethylation (H3K9me3)   | Heterochromatin regions  [@peters2003partitioning]                                                  |
| Histone H3 acetylated at lysine 27 (H3K27ac)   | Increase activation of genomic elements [@heintzman2009histone,@rada2011unique,@creyghton2010histone] |
| Histone H3 lysine 9 acetylation  (H3K9ac)      | Transcriptional activation [@nishida2006histone]                                                    |
"
cat(tabl) 


tabl <- "  
|                File                |                               Description                              |
|:----------------------------------:|:----------------------------------------------------------------------:|
| fc.signal.bigwig                   | Bigwig File containing  fold enrichment signal tracks                  |
| pval.signal.bigwig                 | Bigwig File containing -log10(p-value) signal tracks                   |
| hotspot.fdr0.01.broad.bed.gz       | Broad domains on enrichment for  DNase-seq for consolidated epigenomes |
| hotspot.broad.bed.gz               | Broad domains on enrichment for DNase-seq for consolidated epigenomes  |
| broadPeak.gz                       | Broad ChIP-seq peaks for consolidated  epigenomes                      |
| gappedPeak.gz                      | Gapped ChIP-seq peaks for consolidated   epigenomes                    |
| narrowPeak.gz                      | Narrow ChIP-seq peaks for consolidated epigenomes                      |
| hotspot.fdr0.01.peaks.bed.gz       | Narrow DNasePeaks for   consolidated epigenomes                        |
| hotspot.all.peaks.bed.gz           | Narrow DNasePeaks for  consolidated epigenomes                         |
| .macs2.narrowPeak.gz               | Narrow DNasePeaks for consolidated epigenomes                          |
| coreMarks_mnemonics.bed.gz        | 15 state chromatin segmentations                                       |
| mCRF_FractionalMethylation.bigwig | MeDIP/MRE(mCRF) fractional methylation calls                           |
| RRBS_FractionalMethylation.bigwig | RRBS fractional methylation calls                                      |
| WGBS_FractionalMethylation.bigwig | Whole genome bisulphite fractional methylation calls                   |
"
cat(tabl) 

dev.off()

sessionInfo()