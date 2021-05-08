url <- "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/dlbc_tcga_pan_can_atlas_2018/data_RNA_Seq_v2_expression_median.txt"
destfile <- "data/dlbc_exp.txt"
download.file(url, destfile)


url <- "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/dlbc_tcga_pan_can_atlas_2018/data_clinical_patient.txt"
destfile <- "data/data_clinical.txt"
download.file(url, destfile)


data_clinical <- read.delim("data/data_clinical.txt", comment.char="#", na.strings=c("NA", ""))

dlbc_exp <- read.delim("data/dlbc_exp.txt", comment.char="#", na.strings=c("NA", ""))

dlbc_exp <- dlbc_exp[!is.na(dlbc_exp$Hugo_Symbol),]
dlbc_exp <- dlbc_exp[!duplicated(dlbc_exp$Hugo_Symbol),]
rownames(dlbc_exp) <- as.character(dlbc_exp$Hugo_Symbol)
dlbc_exp$Entrez_Gene_Id <- NULL
dlbc_exp$Hugo_Symbol <- NULL


#data_clinical
