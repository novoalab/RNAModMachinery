### Cancer and Normal RNA Expression analysis using the stage information

#Rscript cancer_script6_stageexpression.R RMLP.TcgaTargetGtex_rsem_gene_tpm_withheader.tsv human_id_symbol_class.tsv clinical.tsv TcgaTargetGTEX_phenotype.txt
args <- commandArgs(trailingOnly = TRUE) #Argument for first input
input1 <- args[1]#1st variable 
data<-read.delim(input1)#RMLP.TcgaTargetGtex_rsem_gene_tpm_withheader.tsv
library(dplyr) #load the package for data manipulation
library(plyr)
data$sample <- gsub("\\..*","",data$sample) #remove everything after "." in ENSEMBL IDs
args <- commandArgs(trailingOnly = TRUE) #Argument for first input
input2 <- args[2]#2nd variable #human_id_symbol_class.tsv 
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
library(data.table)
tcga<- transposed2[transposed2$sample %like% "TCGA", ] #7792 GTEX sample
tcga$sample2<- substr(tcga$sample,1,nchar(tcga$sample)-3)
args <- commandArgs(trailingOnly = TRUE) #Argument for first input
input3 <- args[3]#2nd variable 
clinical<- read.delim(input3)#clinical.tsv
clinical<- clinical[,c("submitter_id", "ajcc_pathologic_stage")]
colnames(clinical)<- c("sample2", "stage")
withclinical<- join(tcga, clinical, by="sample2") #join two tables by their matching sample names (VLOOKUP)
withclinical <- withclinical %>% select(sample,sample2,stage, everything())
withclinical2 <- data.frame(sapply(withclinical[,-c(1:3)], function(x) as.numeric(as.character(x)))) #MAKE NUMERIC
withclinical3<- cbind(withclinical[,c(1:3)], withclinical2) 
withclinical3$stage<- gsub("Stage IA","Stage I",withclinical3$stage)
withclinical3$stage<- gsub("Stage IS","Stage I",withclinical3$stage)
withclinical3$stage<- gsub("Stage IB","Stage I",withclinical3$stage)
withclinical3$stage<- gsub("Stage IIA","Stage II",withclinical3$stage)
withclinical3$stage<- gsub("Stage IIB","Stage II",withclinical3$stage)
withclinical3$stage<- gsub("Stage IIC","Stage II",withclinical3$stage)
withclinical3$stage<- gsub("Stage IIIA","Stage III",withclinical3$stage)
withclinical3$stage<- gsub("Stage IIIB","Stage III",withclinical3$stage)
withclinical3$stage<- gsub("Stage IIIC","Stage III",withclinical3$stage)
withclinical3$stage<- gsub("Stage IVA","Stage IV",withclinical3$stage)
withclinical3$stage<- gsub("Stage IVB","Stage IV",withclinical3$stage)
withclinical3$stage<- gsub("Stage IVC","Stage IV",withclinical3$stage)

write.table(withclinical3,file="RMLP.TcgaTargetGtex_tpm_with_clinical.tsv",row.names=FALSE, quote=FALSE, sep="\t",col.names=TRUE) #export


args <- commandArgs(trailingOnly = TRUE) #Argument for first input
input4 <- args[4]#2nd variable 
phenotype<- read.table(input4,sep="\t",header=TRUE) #TcgaTargetGTEX_phenotype.txt
phenotype$primary.disease.or.tissue <- NULL #remove unnecessary columns
phenotype$X_primary_site<- NULL #remove unnecessary columns
phenotype$X_gender<- NULL#remove unnecessary columns
phenotype$X_study<- NULL#remove unnecessary columns
colnames(phenotype) <- gsub("X_","",colnames(phenotype)) #change column names to simpler ones
colnames(phenotype) <- gsub("_","",colnames(phenotype)) #change column names to simpler ones

withphenotypeclinical<- join(withclinical3, phenotype, by="sample") #join two tables by their matching sample names (VLOOKUP)
withphenotypeclinical <- withphenotypeclinical %>% select(sample,detailedcategory,sampletype, everything())
tcga<-withphenotypeclinical

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


tcga2<-rbind(ACC_TCGA,BLCA_TCGA,BRCA_TCGA,CESC_TCGA,CHOL_TCGA,COAD_TCGA,DLBC_TCGA,ESCA_TCGA,GBM_TCGA,HNSC_TCGA,KICH_TCGA,KIRC_TCGA,KIRP_TCGA,LAML_TCGA,LGG_TCGA,LIHC_TCGA,LUAD_TCGA,LUSC_TCGA,OV_TCGA,PAAD_TCGA,PCPG_TCGA,PRAD_TCGA,READ_TCGA,SARC_TCGA,SKCM_TCGA,STAD_TCGA,TGCT_TCGA,THCA_TCGA,THYM_TCGA,UCEC_TCGA,UCS_TCGA)

tcga2$sampletype <- gsub("Additional - New Primary","Tumor",tcga2$sampletype) #change column names to simpler ones
tcga2$sampletype <- gsub("Additional Metastatic","Tumor",tcga2$sampletype) #change column names to simpler ones
tcga2$sampletype <- gsub("Additional Tumor","Tumor",tcga2$sampletype) #change column names to simpler ones
tcga2$sampletype <- gsub("Metastatic","Tumor",tcga2$sampletype) #change column names to simpler ones
tcga2$sampletype <- gsub("Primary Blood Derived Cancer - Peripheral Blood","Tumor",tcga2$sampletype) #change column names to simpler ones
tcga2$sampletype <- gsub("Primary Tumor","Tumor",tcga2$sampletype) #change column names to simpler ones
tcga2$sampletype <- gsub("Recurrent Tumor","Tumor",tcga2$sampletype) #change column names to simpler ones
tcga2$sampletype <- gsub("Solid Tissue Normal","Normal",tcga2$sampletype) #change column names to simpler ones
tcga2$sampletype <- gsub("Normal Tissue","Normal",tcga2$sampletype) #change column names to simpler ones
tcga2$sampletype <- gsub("Cell Line","Normal",tcga2$sampletype) #change column names to simpler ones

tcga2 <- tcga2 %>% select(Type, everything())

tcga2$sample2 <- NULL

tcga2normal<- subset(tcga2, tcga2$sampletype=="Normal")
tcga2normal$stage<- rep("Normal",nrow(tcga2normal))
tcga2tumor<-  subset(tcga2, tcga2$sampletype=="Tumor")
tcga3<- rbind (tcga2normal, tcga2tumor)





###GTEX
#args <- commandArgs(trailingOnly = TRUE) #Argument for first input
#input2 <- args[2]#2nd variable #human_id_symbol_class.tsv 
withphenotype<- join(transposed2, phenotype, by="sample") #join two tables by their matching sample names (VLOOKUP)
withphenotype <- withphenotype %>% select(sample,detailedcategory,sampletype, everything())
withphenotype2 <- data.frame(sapply(withphenotype[,-c(1:3)], function(x) as.numeric(as.character(x)))) #MAKE NUMERIC
withphenotype3<- cbind(withphenotype[,c(1:3)], withphenotype2) 

library(data.table)
gtex<- withphenotype3[withphenotype3$sample %like% "GTEX", ] #7792 GTEX sample
k562<- withphenotype3[withphenotype3$sample %like% "K-562", ] #7792 GTEX sample

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


gtex2<- rbind(ACC_GTEX,BLCA_GTEX,BRCA_GTEX,CESC_GTEX,COAD_GTEX,DLBC_GTEX,ESCA_GTEX,GBM_GTEX,KICH_GTEX,KIRC_GTEX,KIRP_GTEX,LAML_GTEX,LGG_GTEX,LIHC_GTEX,LUAD_GTEX,LUSC_GTEX,OV_GTEX,PAAD_GTEX,PRAD_GTEX,READ_GTEX,SKCM_GTEX,STAD_GTEX,TGCT_GTEX,THCA_GTEX,THYM_GTEX,UCEC_GTEX,UCS_GTEX)
gtex2$stage<- rep("Normal",nrow(gtex2))
gtex2 <- gtex2 %>% select(Type,sample,detailedcategory,sampletype, stage, everything()) #place the last column to first column


tcga_gtex<-rbind(tcga3,gtex2)
tcga_gtex2<-tcga_gtex[order(tcga_gtex$MRM3),] #Sort by any gene to see if there is any data point that has off pattern
tcga_gtex2<- subset(tcga_gtex2, tcga_gtex2$MRM3 > -9.9657) #Remove the rows that has -9.9658 values
tcga_gtex2<-tcga_gtex2[order(tcga_gtex2$sampletype),] #Sort by sampletype
tcga_gtex2<-tcga_gtex2[order(tcga_gtex2$Type),] #Sort by Type





#R, Initially we had log2(TPM+0.001) but we want log2(TPM+1). So we will do the transformation
dat<-tcga_gtex2
#dat2 is 2 to the power X, which equals to TPM + 0.001
dat2 <- 2^dat[,-c(1:5)]
#dat2 is TPM (with an extra 0.001 for each datapoint)
#Export this as well
tpmfile <- cbind (dat[,1:5],dat2) #merge the initial columns before exporting
#write.table(tpmfile,file="TCGA_GTEX_FINAL.TPM.tsv",quote=FALSE, sep="\t",row.names=FALSE)
#dat3 is log2(TPM+1)
dat3<- log2(dat2+1)
#Merge the processed log2(TPM+1) with the first initial two columns 
logfile <- cbind (dat[,1:5],dat3)
#write.table(logfile,file="TCGA_GTEX_FINAL.log2.tsv",quote=FALSE, sep="\t",row.names=FALSE)



#IF you want to remove some cancer types from your file
remove <- c("DLBC", "CHOL", "THYM")
tpmfile2<- tpmfile[!tpmfile$Type %in% remove, ]
#write.table(tpmfile2,file="TCGA_GTEX_FINAL.TPM.without3cancer.tsv",quote=FALSE, sep="\t",row.names=FALSE)
logfile2<- logfile[!logfile$Type %in% remove, ]
#write.table(logfile2,file="TCGA_GTEX_FINAL.log2.without3cancer.tsv",quote=FALSE, sep="\t",row.names=FALSE)



stages<-c("Stage IV","Stage III","Stage II","Stage I","Normal")
final<- logfile2[logfile2$stage %in% stages, ]
final$sampletype <- gsub("Normal Tissue","Normal",final$sampletype) #change column names to simpler ones
final$sampletype <- gsub("Cell Line","Normal",final$sampletype) #change column names to simpler ones




#Plot violin plots for LAGE3
library(ggplot2)
pdf(file="HENMT1.violin.stage.pdf",height=50,width=7)
print(ggplot(data=final, aes(x=stage, y=HENMT1, fill=stage)) +
	geom_violin(aes(fill = stage),scale="width",draw_quantiles = c(0.5))+
	geom_jitter(height = 0, width = 0.1,size=0.2)+
	facet_grid(factor(Type) ~ .)+
	scale_fill_manual(values=c("#fc5185","#ff9e74","#e3c4a8","#3fc1c9","#364f6b","#1b262c"))+
     theme(panel.background = element_blank(),
         panel.border=element_rect(fill=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background=element_blank(),
         axis.text.x=element_text(colour="black",angle=90),
         axis.text.y=element_text(colour="black"),
         axis.ticks=element_line(colour="black"),
         plot.margin=unit(c(1,1,1,1),"line")))
dev.off()






##Eva suggested the following :

#What about taking the median for each cancer type, making that normal is value 1, to just then see if the trend is to go up or go down?
#So a line plot of the median values for each cancer type
#Normalized to the normal such that all cancers start in normal=1


final2<- subset(final, final$Type!="ACC")
final2<- subset(final2, final2$Type!="CESC")
final2<- subset(final2, final2$Type!="GBM")
final2<- subset(final2, final2$Type!="LAML")
final2<- subset(final2, final2$Type!="LGG")
final2<- subset(final2, final2$Type!="LAML")
final2<- subset(final2, final2$Type!="OV")
final2<- subset(final2, final2$Type!="UCS")
final2<- subset(final2, final2$Type!="UCEC")
final2<- subset(final2, final2$Type!="PCPG")
final2<- subset(final2, final2$Type!="SARC")
final2<- subset(final2, final2$Type!="PRAD")
dat<-final2



library(ggpubr) 
library(EnvStats)
library(scales)
library(ggplot2)


for (i in c(5:dim(dat)[2])) {
pdf(file=paste(colnames(dat)[i],"violinstageplot.logtpm.pdf",sep="."),height=35,width=5)
print(ggplot(data=dat, aes_string(x=colnames(dat)[5], y=noquote(paste(colnames(dat)[i])), fill=colnames(dat)[5])) +
	geom_violin(aes(fill = stage),scale="width",draw_quantiles = c(0.5))+
	geom_jitter(height = 0, width = 0.1,size=0.2)+
	facet_grid(factor(Type) ~., scale="free")+
	scale_fill_manual(values=c("#fc5185","#ff9e74","#e3c4a8","#3fc1c9","#364f6b","#1b262c"))+
     theme(panel.background = element_blank(),
         panel.border=element_rect(fill=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background=element_blank(),
         axis.text.x=element_text(colour="black",angle=90),
         axis.text.y=element_text(colour="black"),
         axis.ticks=element_line(colour="black"),
         plot.margin=unit(c(1,1,1,1),"line")))
dev.off()
}




#Line plots for each gene
#library(data.table)
#for (i in c(6:dim(dat)[2])) {
#genea<- final2[,c("Type", "sampletype", "stage", names(dat[i]))]
#genea2 <- data.table(genea)
#genea2<- genea2[, median:=median(dat[,i]), by=list(Type, stage)]
#genea3<- subset(genea2, !duplicated(subset(genea2, select=c(Type, stage))))
#genea3[,c(names(dat[i]))]<- NULL
#geneanorm<- subset(genea3, genea3$sampletype=="Normal")
#names(geneanorm)[4]<- c("normedian")
#genea4<-join(genea3,geneanorm, by="Type")
#genea5<- genea4[,c(1,2,3,4,7)]
#genea5$normalized<- genea5$median/genea5$normedian
#pdf(file=paste(names(dat[i]), "line.stage.grouped.pdf", sep="_"),height=5,width=5)
#print(ggplot(data=genea5, aes(x=stage, y=normalized, group=Type)) +
#  geom_line(aes(color=Type))+
#  geom_point(size=0.5)+
#      theme(panel.background = element_blank(),
#         panel.border=element_rect(fill=NA),
#        panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.background=element_blank(),
#         axis.text.x=element_text(colour="black",angle=90),
#         axis.text.y=element_text(colour="black"),
#        axis.ticks=element_line(colour="black"),
#         plot.margin=unit(c(1,1,1,1),"line")))
#dev.off()
#}






#Do it for only HENMT1
henmt<- final2[,c("Type", "sampletype", "stage","HENMT1")]
library(data.table)
henmt2 <- data.table(henmt)
henmt2<- henmt2[, median:=median(HENMT1), by=list(Type, stage)]
henmt3<- subset(henmt2, !duplicated(subset(henmt2, select=c(Type, stage))))
henmt3$HENMT1<- NULL
henmtnorm<- subset(henmt3, henmt3$sampletype=="Normal")
names(henmtnorm)[4]<- c("normedian")
henmt4<-join(henmt3,henmtnorm, by="Type")
henmt5<- henmt4[,c(1,2,3,4,7)]
henmt5$normalized<- henmt5$median/henmt5$normedian


pdf(file="HENMT1.line.stage.grouped.pdf",height=5,width=5)
print(ggplot(data=henmt5, aes(x=stage, y=normalized, group=Type)) +
  geom_line(aes(color=Type))+
  geom_point(size=0.5)+
      theme(panel.background = element_blank(),
         panel.border=element_rect(fill=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background=element_blank(),
         axis.text.x=element_text(colour="black",angle=90),
         axis.text.y=element_text(colour="black"),
         axis.ticks=element_line(colour="black"),
         plot.margin=unit(c(1,1,1,1),"line")))
dev.off()





#Do it for only lage3
lage<- final2[,c("Type", "sampletype", "stage","LAGE3")]
library(data.table)
lage2 <- data.table(lage)
lage2<- lage2[, median:=median(LAGE3), by=list(Type, stage)]
lage3<- subset(lage2, !duplicated(subset(lage2, select=c(Type, stage))))
lage3$LAGE3<- NULL


lagenorm<- subset(lage3, lage3$sampletype=="Normal")
names(lagenorm)[4]<- c("normedian")
lage4<-join(lage3,lagenorm, by="Type")
lage5<- lage4[,c(1,2,3,4,7)]
lage5$normalized<- lage5$median/lage5$normedian


pdf(file="LAGE3.line.stage.grouped.pdf",height=5,width=5)
print(ggplot(data=lage5, aes(x=stage, y=normalized, group=Type)) +
  geom_line(aes(color=Type))+
  geom_point(size=0.5)+
      theme(panel.background = element_blank(),
         panel.border=element_rect(fill=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background=element_blank(),
         axis.text.x=element_text(colour="black",angle=90),
         axis.text.y=element_text(colour="black"),
         axis.ticks=element_line(colour="black"),
         plot.margin=unit(c(1,1,1,1),"line")))
dev.off()

