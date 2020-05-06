###################################################################
## TCGA-GTEX DATA MANIPULATION
## 2020, Oguzhan Begik written for Begik et al, 2020 Genome Biology
###################################################################
#Original file (With Sample ID and expression in ENSEMBL IDs)
#TcgaTargetGtex_rsem_gene_tpm
#Phenotype file
#TcgaTargetGTEX_phenotype.txt (Sample ID and phenotype)

#Data from the study is from the UCSC RNA-seq Compendium, where TCGA and GTEx 
#samples are re-analyzed (re-aligned to hg38 genome and expressions are called using RSEM and Kallisto methods)
#by the same RNA-seq pipeline. Because all samples are processed using an uniform bioinformatic pipeline, 
#batch effect due to different computational processing is eliminated.
#Data (file names: *.rsem_genes.results) are downloaded, tpm values are extracted, log2(x+0.001) transformed, and combined.
#extract modomics from the expression file
#rwp_tRNA_all_ensembl contains ensembl ids of all RMWs
#TcgaTargetGtex_rsem_gene_tpm contains log2(tpm+0.001) transformed expression values for 

#How to run the script
#Rscript cancer_script1_datamanipulation.R <TCGA.GTEX.file> <ENSEMBL.file>

#Example:  Rscript cancer_script1_datamanipulation.R RMLP.TcgaTargetGtex_rsem_gene_tpm_withheader.tsv human_id_symbol_class.tsv


###Requirements
#libraries: 
#dplyr
#plyr
#reshape2
#data.table





# 1. Arguments introduced for the execution
########################################

args <- commandArgs(trailingOnly = TRUE) #Argument for first input
input1 <- args[1]#1st variable 
input2 <- args[2] #2nd variable



# 1.Importing the data and manipulating
########################################



#Load the library
library(dplyr)
library(plyr)
library(reshape2)
library(data.table)






data<-read.delim(input1)#RMLP.TcgaTargetGtex_rsem_gene_tpm_withheader.tsv
data$sample <- gsub("\\..*","",data$sample) #remove everything after "." in ENSEMBL IDs
ensembl <- read.table(input2, sep="\t",header=TRUE)#import the ensembl file that contains ENSEMBL ID and matching GeneNames
colnames(ensembl)<- gsub("gene_id","sample",colnames(ensembl)) 
joined<- join(data, ensembl, by="sample") #use join function to add a column of gene names corresponding to the ENSEMBL ID to the last column
replaced <- joined %>% select(Symbol, everything()) #place the last column to first column
replaced <- replaced %>% select(Class, everything()) #place the last column to first column
replaced$Class<-NULL
replaced$sample<-NULL
transposed<-t(replaced)
colnames(transposed)<-transposed[1,]
transposed<-transposed[-1,]
row.names(transposed) <- gsub("\\.","-",row.names(transposed))
transposed2<- as.data.frame(transposed)
transposed2$sample<-row.names(transposed2)
phenotype<- read.table("TcgaTargetGTEX_phenotype.txt",sep="\t",header=TRUE) #read the phenotype table
phenotype$primary.disease.or.tissue <- NULL #remove unnecessary columns
phenotype$X_primary_site<- NULL #remove unnecessary columns
phenotype$X_gender<- NULL#remove unnecessary columns
phenotype$X_study<- NULL#remove unnecessary columns
colnames(phenotype) <- gsub("X_","",colnames(phenotype)) #change column names to simpler ones
colnames(phenotype) <- gsub("_","",colnames(phenotype)) #change column names to simpler ones
withphenotype<- join(transposed2, phenotype, by="sample") #join two tables by their matching sample names (VLOOKUP)
withphenotype <- withphenotype %>% select(sample,detailedcategory,sampletype, everything())
withphenotype2 <- data.frame(sapply(withphenotype[,-c(1:3)], function(x) as.numeric(as.character(x)))) #MAKE NUMERIC
withphenotype3<- cbind(withphenotype[,c(1:3)], withphenotype2) 
write.table(withphenotype3,file="RMLP.TcgaTargetGtex_tpm_with_phenotype.tsv",row.names=FALSE, quote=FALSE, sep="\t",col.names=TRUE) #export


gtex<- withphenotype3[withphenotype3$sample %like% "GTEX", ] # GTEX sample
tcga<- withphenotype3[withphenotype3$sample %like% "TCGA", ] # TCGA sample
k562<- withphenotype3[withphenotype3$sample %like% "K-562", ] # K562 cell lines (TARGET) sample

# Extract table for each cancer type seperately
ACC_TCGA<-tcga[tcga$detailedcategory %like% "Adrenocortical Cancer", ]
BLCA_TCGA<-tcga[tcga$detailedcategory %like% "Bladder Urothelial Carcinoma", ]
BRCA_TCGA<-tcga[tcga$detailedcategory %like% "Breast Invasive Carcinoma", ]
CESC_TCGA<-tcga[tcga$detailedcategory %like% "Cervical & Endocervical Cancer", ]
CHOL_TCGA<-tcga[tcga$detailedcategory %like% "Cholangiocarcinoma", ]
COAD_TCGA<-tcga[tcga$detailedcategory %like% "Colon Adenocarcinoma", ]
DLBC_TCGA<-tcga[tcga$detailedcategory %like% "Diffuse Large B-Cell Lymphoma", ]
ESCA_TCGA<-tcga[tcga$detailedcategory %like% "Esophageal Carcinoma", ]
GBM_TCGA<-tcga[tcga$detailedcategory %like% "Glioblastoma Multiforme", ]
HNSC_TCGA<-tcga[tcga$detailedcategory %like% "Head & Neck Squamous Cell Carcinoma", ]
KICH_TCGA<-tcga[tcga$detailedcategory %like% "Kidney Chromophobe", ]
KIRC_TCGA<-tcga[tcga$detailedcategory %like% "Kidney Clear Cell Carcinoma", ]
KIRP_TCGA<-tcga[tcga$detailedcategory %like% "Kidney Papillary Cell Carcinoma", ]
LAML_TCGA<-tcga[tcga$detailedcategory %like% "Acute Myeloid Leukemia", ]
LGG_TCGA<-tcga[tcga$detailedcategory %like% "Brain Lower Grade Glioma", ]
LIHC_TCGA<-tcga[tcga$detailedcategory %like% "Liver Hepatocellular Carcinoma", ]
LUAD_TCGA<-tcga[tcga$detailedcategory %like% "Lung Adenocarcinoma", ]
LUSC_TCGA<-tcga[tcga$detailedcategory %like% "Lung Squamous Cell Carcinoma", ]
OV_TCGA<-tcga[tcga$detailedcategory %like% "Ovarian Serous Cystadenocarcinoma", ]
PAAD_TCGA<-tcga[tcga$detailedcategory %like% "Pancreatic Adenocarcinoma", ]
PCPG_TCGA<-tcga[tcga$detailedcategory %like% "Pheochromocytoma & Paraganglioma", ]
PRAD_TCGA<-tcga[tcga$detailedcategory %like% "Prostate Adenocarcinoma", ]
READ_TCGA<-tcga[tcga$detailedcategory %like% "Rectum Adenocarcinoma", ]
SARC_TCGA<-tcga[tcga$detailedcategory %like% "Sarcoma", ]
SKCM_TCGA<-tcga[tcga$detailedcategory %like% "Skin Cutaneous Melanoma", ]
STAD_TCGA<-tcga[tcga$detailedcategory %like% "Stomach Adenocarcinoma", ]
TGCT_TCGA<-tcga[tcga$detailedcategory %like% "Testicular Germ Cell Tumor", ]
THCA_TCGA<-tcga[tcga$detailedcategory %like% "Thyroid Carcinoma", ]
THYM_TCGA<-tcga[tcga$detailedcategory %like% "Thymoma", ]
UCEC_TCGA<-tcga[tcga$detailedcategory %like% "Uterine Corpus Endometrioid Carcinoma", ]
UCS_TCGA<-tcga[tcga$detailedcategory %like% "Uterine Carcinosarcoma", ]


# Add a repeating column
ACC_TCGA$Type <- rep("ACC",nrow(ACC_TCGA))
BLCA_TCGA$Type <- rep("BLCA",nrow(BLCA_TCGA))
BRCA_TCGA$Type <- rep("BRCA",nrow(BRCA_TCGA))
CESC_TCGA$Type <- rep("CESC",nrow(CESC_TCGA))
CHOL_TCGA$Type <- rep("CHOL",nrow(CHOL_TCGA))
COAD_TCGA$Type <- rep("COAD",nrow(COAD_TCGA))
DLBC_TCGA$Type <- rep("DLBC",nrow(DLBC_TCGA))
ESCA_TCGA$Type <- rep("ESCA",nrow(ESCA_TCGA))
GBM_TCGA$Type <- rep("GBM",nrow(GBM_TCGA))
HNSC_TCGA$Type <- rep("HNSC",nrow(HNSC_TCGA))
KICH_TCGA$Type <- rep("KICH",nrow(KICH_TCGA))
KIRC_TCGA$Type <- rep("KIRC",nrow(KIRC_TCGA))
KIRP_TCGA$Type <- rep("KIRP",nrow(KIRP_TCGA))
LAML_TCGA$Type <- rep("LAML",nrow(LAML_TCGA))
LGG_TCGA$Type <- rep("LGG",nrow(LGG_TCGA))
LIHC_TCGA$Type <- rep("LIHC",nrow(LIHC_TCGA))
LUAD_TCGA$Type <- rep("LUAD",nrow(LUAD_TCGA))
LUSC_TCGA$Type <- rep("LUSC",nrow(LUSC_TCGA))
OV_TCGA$Type <- rep("OV",nrow(OV_TCGA))
PAAD_TCGA$Type <- rep("PAAD",nrow(PAAD_TCGA))
PCPG_TCGA$Type <- rep("PCPG",nrow(PCPG_TCGA))
PRAD_TCGA$Type <- rep("PRAD",nrow(PRAD_TCGA))
READ_TCGA$Type <- rep("READ",nrow(READ_TCGA))
SARC_TCGA$Type <- rep("SARC",nrow(SARC_TCGA))
SKCM_TCGA$Type <- rep("SKCM",nrow(SKCM_TCGA))
STAD_TCGA$Type <- rep("STAD",nrow(STAD_TCGA))
TGCT_TCGA$Type <- rep("TGCT",nrow(TGCT_TCGA))
THCA_TCGA$Type <- rep("THCA",nrow(THCA_TCGA))
THYM_TCGA$Type <- rep("THYM",nrow(THYM_TCGA))
UCEC_TCGA$Type <- rep("UCEC",nrow(UCEC_TCGA))
UCS_TCGA$Type <- rep("UCS",nrow(UCS_TCGA))

# Extract table for each tissue type seperately
ACC_GTEX<-gtex[gtex$detailedcategory %like% "Adrenal", ]
BLCA_GTEX<-gtex[gtex$detailedcategory %like% "Bladder", ]
BRCA_GTEX<-gtex[gtex$detailedcategory %like% "Breast", ]
CESC_GTEX<-gtex[gtex$detailedcategory %like% "Cervix", ]
COAD_GTEX<-gtex[gtex$detailedcategory %like% "Colon", ]
DLBC_GTEX<-gtex[gtex$detailedcategory %like% "Blood", ]
ESCA_GTEX<-gtex[gtex$detailedcategory %like% "Esophagus - Mucosa", ]
GBM_GTEX<-gtex[gtex$detailedcategory %like% "Brain - Cortex|Frontal Cortex", ]
KICH_GTEX<-gtex[gtex$detailedcategory %like% "Kidney", ]
KIRC_GTEX<-gtex[gtex$detailedcategory %like% "Kidney", ]
KIRP_GTEX<-gtex[gtex$detailedcategory %like% "Kidney", ]
LAML_GTEX<-k562[k562$detailedcategory %like% "Cells", ]
LGG_GTEX<-gtex[gtex$detailedcategory %like% "Brain - Cortex|Frontal Cortex", ]
LIHC_GTEX<-gtex[gtex$detailedcategory %like% "Liver", ]
LUAD_GTEX<-gtex[gtex$detailedcategory %like% "Lung", ]
LUSC_GTEX<-gtex[gtex$detailedcategory %like% "Lung", ]
OV_GTEX<-gtex[gtex$detailedcategory %like% "Ovary", ]
PAAD_GTEX<-gtex[gtex$detailedcategory %like% "Pancreas", ]
PRAD_GTEX<-gtex[gtex$detailedcategory %like% "Prostate", ]
READ_GTEX<-gtex[gtex$detailedcategory %like% "Colon", ]
SKCM_GTEX<-gtex[gtex$detailedcategory %like% "Skin", ]
STAD_GTEX<-gtex[gtex$detailedcategory %like% "Stomach", ]
TGCT_GTEX<-gtex[gtex$detailedcategory %like% "Testis", ]
THCA_GTEX<-gtex[gtex$detailedcategory %like% "Thyroid", ]
THYM_GTEX<-gtex[gtex$detailedcategory %like% "Blood", ]
UCEC_GTEX<-gtex[gtex$detailedcategory %like% "Uterus", ]
UCS_GTEX<-gtex[gtex$detailedcategory %like% "Uterus", ]


# Add a repeating column
ACC_GTEX$Type <- rep("ACC",nrow(ACC_GTEX))
BLCA_GTEX$Type <- rep("BLCA",nrow(BLCA_GTEX))
BRCA_GTEX$Type <- rep("BRCA",nrow(BRCA_GTEX))
CESC_GTEX$Type <- rep("CESC",nrow(CESC_GTEX))
COAD_GTEX$Type <- rep("COAD",nrow(COAD_GTEX))
DLBC_GTEX$Type <- rep("DLBC",nrow(DLBC_GTEX))
ESCA_GTEX$Type <- rep("ESCA",nrow(ESCA_GTEX))
GBM_GTEX$Type <- rep("GBM",nrow(GBM_GTEX))
KICH_GTEX$Type <- rep("KICH",nrow(KICH_GTEX))
KIRC_GTEX$Type <- rep("KIRC",nrow(KIRC_GTEX))
KIRP_GTEX$Type <- rep("KIRP",nrow(KIRP_GTEX))
LAML_GTEX$Type <- rep("LAML",nrow(LAML_GTEX))
LGG_GTEX$Type <- rep("LGG",nrow(LGG_GTEX))
LIHC_GTEX$Type <- rep("LIHC",nrow(LIHC_GTEX))
LUAD_GTEX$Type <- rep("LUAD",nrow(LUAD_GTEX))
LUSC_GTEX$Type <- rep("LUSC",nrow(LUSC_GTEX))
OV_GTEX$Type <- rep("OV",nrow(OV_GTEX))
PAAD_GTEX$Type <- rep("PAAD",nrow(PAAD_GTEX))
PRAD_GTEX$Type <- rep("PRAD",nrow(PRAD_GTEX))
READ_GTEX$Type <- rep("READ",nrow(READ_GTEX))
SKCM_GTEX$Type <- rep("SKCM",nrow(SKCM_GTEX))
STAD_GTEX$Type <- rep("STAD",nrow(STAD_GTEX))
TGCT_GTEX$Type <- rep("TGCT",nrow(TGCT_GTEX))
THCA_GTEX$Type <- rep("THCA",nrow(THCA_GTEX))
THYM_GTEX$Type <- rep("THYM",nrow(THYM_GTEX))
UCEC_GTEX$Type <- rep("UCEC",nrow(UCEC_GTEX))
UCS_GTEX$Type <- rep("UCS",nrow(UCS_GTEX))


#Bind the TCGA and GTEX subtables
tcga2<-rbind(ACC_TCGA,BLCA_TCGA,BRCA_TCGA,CESC_TCGA,CHOL_TCGA,COAD_TCGA,DLBC_TCGA,ESCA_TCGA,GBM_TCGA,HNSC_TCGA,KICH_TCGA,KIRC_TCGA,KIRP_TCGA,LAML_TCGA,LGG_TCGA,LIHC_TCGA,LUAD_TCGA,LUSC_TCGA,OV_TCGA,PAAD_TCGA,PCPG_TCGA,PRAD_TCGA,READ_TCGA,SARC_TCGA,SKCM_TCGA,STAD_TCGA,TGCT_TCGA,THCA_TCGA,THYM_TCGA,UCEC_TCGA,UCS_TCGA)
gtex2<- rbind(ACC_GTEX,BLCA_GTEX,BRCA_GTEX,CESC_GTEX,COAD_GTEX,DLBC_GTEX,ESCA_GTEX,GBM_GTEX,KICH_GTEX,KIRC_GTEX,KIRP_GTEX,LAML_GTEX,LGG_GTEX,LIHC_GTEX,LUAD_GTEX,LUSC_GTEX,OV_GTEX,PAAD_GTEX,PRAD_GTEX,READ_GTEX,SKCM_GTEX,STAD_GTEX,TGCT_GTEX,THCA_GTEX,THYM_GTEX,UCEC_GTEX,UCS_GTEX)
#Bind both TCGA and GTEX
tcga_gtex<-rbind(tcga2,gtex2)
tcga_gtex <- tcga_gtex %>% select(Type, everything()) #place the last column to first column
#write.table(tcga_gtex,file="TCGA_GTEX_FINAL.forR.tsv",quote=FALSE,sep="\t",row.names=FALSE)


#tcga_gtex<- read.delim("TCGA_GTEX_FINAL.forR.tsv")
tcga_gtex2<-tcga_gtex[order(tcga_gtex$MRM3),] #Sort by any gene to see if there is any data point that has off pattern
tcga_gtex2<- subset(tcga_gtex2, tcga_gtex2$MRM3 > -9.9657) #Remove the rows that has -9.9658 values
tcga_gtex2<-tcga_gtex2[order(tcga_gtex2$sampletype),] #Sort by sampletype
tcga_gtex2<-tcga_gtex2[order(tcga_gtex2$Type),] #Sort by Type


tcga_gtex2$sampletype <- gsub("Additional - New Primary","Tumor",tcga_gtex2$sampletype) #change column names to simpler ones
tcga_gtex2$sampletype <- gsub("Additional Metastatic","Tumor",tcga_gtex2$sampletype) #change column names to simpler ones
tcga_gtex2$sampletype <- gsub("Additional Tumor","Tumor",tcga_gtex2$sampletype) #change column names to simpler ones
tcga_gtex2$sampletype <- gsub("Metastatic","Tumor",tcga_gtex2$sampletype) #change column names to simpler ones
tcga_gtex2$sampletype <- gsub("Primary Blood Derived Cancer - Peripheral Blood","Tumor",tcga_gtex2$sampletype) #change column names to simpler ones
tcga_gtex2$sampletype <- gsub("Primary Tumor","Tumor",tcga_gtex2$sampletype) #change column names to simpler ones
tcga_gtex2$sampletype <- gsub("Recurrent Tumor","Tumor",tcga_gtex2$sampletype) #change column names to simpler ones
tcga_gtex2$sampletype <- gsub("Solid Tissue Normal","Normal",tcga_gtex2$sampletype) #change column names to simpler ones
tcga_gtex2$sampletype <- gsub("Normal Tissue","Normal",tcga_gtex2$sampletype) #change column names to simpler ones
tcga_gtex2$sampletype <- gsub("Cell Line","Normal",tcga_gtex2$sampletype) #change column names to simpler ones



#R, Initially we had log2(TPM+0.001) but we want log2(TPM+1). So we will do the transformation
dat<-tcga_gtex2
#dat2 is 2 to the power X, which equals to TPM + 0.001
dat2 <- 2^dat[,-c(1:4)]
#dat2 is TPM (with an extra 0.001 for each datapoint)
#Export this as well
tpmfile <- cbind (dat[,1:4],dat2) #merge the initial columns before exporting
write.table(tpmfile,file="TCGA_GTEX_FINAL.TPM.tsv",quote=FALSE, sep="\t",row.names=FALSE)
#dat3 is log2(TPM+1)
dat3<- log2(dat2+1)
#Merge the processed log2(TPM+1) with the first initial two columns 
logfile <- cbind (dat[,1:4],dat3)
write.table(logfile,file="TCGA_GTEX_FINAL.log2.tsv",quote=FALSE, sep="\t",row.names=FALSE)


#IF you want to remove some cancer types from your file
remove <- c("DLBC", "CHOL", "THYM")
tpmfile2<- tpmfile[!tpmfile$Type %in% remove, ]
write.table(tpmfile2,file="TCGA_GTEX_FINAL.TPM.without3cancer.tsv",quote=FALSE, sep="\t",row.names=FALSE)
logfile2<- logfile[!logfile$Type %in% remove, ]
write.table(logfile2,file="TCGA_GTEX_FINAL.log2.without3cancer.tsv",quote=FALSE, sep="\t",row.names=FALSE)




