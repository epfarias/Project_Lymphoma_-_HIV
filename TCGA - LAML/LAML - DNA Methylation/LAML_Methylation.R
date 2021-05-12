############# LAML - DNAMethylation ################

# TCGABiolinks: An R/Bioconductor package for integrative analysis with GDC data
# URL: http://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html
# version 2.18.0

## Installing packages
packages_bioconductor <- c("TCGAbiolinks","maftools","BSgenome.Hsapiens.UCSC.hg38","SummarizedExperiment")
packages_cran <- c("DT", "tidyverse", "stringr", "data.table", "pheatmap","NMF","dplyr")

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

setwd("~/TCGA - LAML/LAML - DNAMethylation")

#-----------------------------------
# STEP 1: Search, download, prepare |
#-----------------------------------
# 1.1 - DNA methylation
# ----------------------------------
query.met <- GDCquery(project = "TCGA-LAML", 
                      legacy = TRUE,
                      data.category = "DNA methylation",
                      platform = "Illumina Human Methylation 450")
GDCdownload(query.met)

laml.met <- GDCprepare(query = query.met,
                      save = TRUE, 
                      save.filename = "lamlDNAmet.rda",
                      summarizedExperiment = TRUE)

#-----------------------------------
# 1.2 - RNA expression
# ----------------------------------
query.exp <- GDCquery(project = "TCGA-LAML", 
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type = "results")
GDCdownload(query.exp)
laml.exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = "lamlExp.rda")

# na.omit
laml.met <- subset(laml.met,subset = (rowSums(is.na(assay(laml.met))) == 0))

# Volcano plot
laml.met <- TCGAanalyze_DMC(laml.met, groupCol = "subtype_MethyLevel",
                           group1 = "CIMP-high",
                           group2="CIMP-low",
                           p.cut = 10^-5,
                           diffmean.cut = 0.25,
                           legend = "State",
                           plot.filename = "CIMP-highvsCIMP-low_metvolcano.png")

#-------------------------------------------------
# 1.3 - DEA - Expression analysis - volcano plot
# ------------------------------------------------
laml.exp.aux <- subset(laml.exp, 
                      select = colData(laml.exp)$subtype_MethyLevel %in% c("CIMP-high","CIMP-low"))

idx <- colData(laml.exp.aux)$subtype_MethyLevel %in% c("CIMP-high")
idx2 <- colData(laml.exp.aux)$subtype_MethyLevel %in% c("CIMP-low")

dataPrep <- TCGAanalyze_Preprocessing(object = laml.exp.aux, cor.cut = 0.6)

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfo,
                                      method = "gcContent")

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  qnt.cut = 0.25,
                                  method='quantile')

dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,idx],
                            mat2 = dataFilt[,idx2],
                            Cond1type = "CIMP-high",
                            Cond2type = "CIMP-low",
                            method = "glmLRT")

TCGAVisualize_volcano(x = dataDEGs$logFC,
                      y = dataDEGs$FDR,
                      filename = "LAML_volcanoexp.png",
                      x.cut = 3,
                      y.cut = 10^-5,
                      names = rownames(dataDEGs),
                      color = c("black","red","darkgreen"),
                      names.size = 2,
                      xlab = " Gene expression fold change (Log2)",
                      legend = "State",
                      title = "Volcano plot (CIMP-high vs CIMP-low)",
                      width = 10)

#------------------------------------------
# 1.4 - Starburst plot
# -----------------------------------------
# If true the argument names of the geenes in circle 
# (biologically significant genes, has a change in gene
# expression and DNA methylation and respects all the thresholds)
# will be shown
# these genes are returned by the function see starburst object after the function is executed
starburst <- TCGAvisualize_starburst(met = laml.met, 
                                     exp = dataDEGs,
                                     genome = "hg19",
                                     group1 = "CIMP-high",
                                     group2 = "CIMP-low",
                                     filename = "starburst.png",
                                     met.platform = "450K",
                                     met.p.cut = 10^-5,
                                     exp.p.cut = 10^-5,
                                     diffmean.cut = 0.25,
                                     logFC.cut = 3,
                                     names = FALSE, 
                                     height=10,
                                     width=15,
                                     dpi=300)
