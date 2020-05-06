###################################################################
## EXTRACTING SPECIFIC GENES FROM GTEx TPM File
## 2020, Oguzhan Begik written for Begik et al, 2020 Genome Biology
###################################################################

###Requirements
#libraries: 
#dplyr
####How to use the script
#Rscript gtex_manipulation.R GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct human_id_symbol_class.tsv



# 1. Arguments introduced for the execution
########################################

args <- commandArgs(trailingOnly = TRUE) #Argument for first input
input1 <- args[1]#1st variable 
input2 <- args[2] #2nd variable



# 2. Importing and manipulating the data
########################################

#Load the library
library(dplyr)


#File name I will use GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct
data<-read.delim(input1) #Import data
colnames(data)<- gsub("\\.","",colnames(data))#Remove ... and . from colnames introduced by R
data$gene_id<-gsub("\\..*","",data$gene_id) #Remove .* from ENSEMBL IDs

#We need to make following renamings for ease in data analysis
#R script for replacements
colnames(data)<- gsub("AdiposeSubcutaneous","Adipose_SUB",colnames(data))
colnames(data)<- gsub("AdiposeVisceralOmentum","Adipose_VOM",colnames(data))
colnames(data)<- gsub("AdrenalGland","Adrenal",colnames(data))
colnames(data)<- gsub("ArteryAorta","Artery_AOR",colnames(data))
colnames(data)<- gsub("ArteryCoronary","Artery_COR",colnames(data))
colnames(data)<- gsub("ArteryTibial","Artery_TIB",colnames(data))
colnames(data)<- gsub("BrainAmygdala","Brain_AMY",colnames(data))
colnames(data)<- gsub("BrainAnteriorcingulatecortexBA24","Brain_ACC",colnames(data))
colnames(data)<- gsub("BrainCaudatebasalganglia","Brain_CBG",colnames(data))
colnames(data)<- gsub("BrainCerebellarHemisphere","Cerebellum_HEM",colnames(data))
colnames(data)<- gsub("BrainCerebellum","Cerebellum_CER",colnames(data))
colnames(data)<- gsub("BrainCortex","Brain_COR",colnames(data))
colnames(data)<- gsub("BrainFrontalCortexBA9","Brain_FRO",colnames(data))
colnames(data)<- gsub("BrainHippocampus","Brain_HIP",colnames(data))
colnames(data)<- gsub("BrainHypothalamus","Brain_HYP",colnames(data))
colnames(data)<- gsub("BrainNucleusaccumbensbasalganglia","Brain_NBG",colnames(data))
colnames(data)<- gsub("BrainPutamenbasalganglia","Brain_PBG",colnames(data))
colnames(data)<- gsub("BrainSpinalcordcervicalc1","Brain_SPI",colnames(data))
colnames(data)<- gsub("BrainSubstantianigra","Brain_SUB",colnames(data))
colnames(data)<- gsub("BreastMammaryTissue","Mammary",colnames(data))
colnames(data)<- gsub("CellsEBVtransformedlymphocytes","Cells_EBV",colnames(data))
colnames(data)<- gsub("CellsTransformedfibroblasts","Cells_FIB",colnames(data))
colnames(data)<- gsub("CervixEctocervix","Cervix_ECT",colnames(data))
colnames(data)<- gsub("CervixEndocervix","Cervix_END",colnames(data))
colnames(data)<- gsub("ColonSigmoid","Colon_SIG",colnames(data))
colnames(data)<- gsub("ColonTransverse","Colon_TRA",colnames(data))
colnames(data)<- gsub("EsophagusGastroesophagealJunction","Esophagus_GAS",colnames(data))
colnames(data)<- gsub("EsophagusMucosa","Esophagus_MUC",colnames(data))
colnames(data)<- gsub("EsophagusMuscularis","Esophagus_MUS",colnames(data))
colnames(data)<- gsub("HeartAtrialAppendage","Heart_ATR",colnames(data))
colnames(data)<- gsub("HeartLeftVentricle","Heart_LVE",colnames(data))
colnames(data)<- gsub("SkinNotSunExposedSuprapubic","Skin_NON",colnames(data))
colnames(data)<- gsub("SkinSunExposedLowerleg","Skin_SUN",colnames(data))
colnames(data)<- gsub("SmallIntestineTerminalIleum","Small_Intestine",colnames(data))
colnames(data)<- gsub("KidneyCortex","Kidney",colnames(data))
colnames(data)<- gsub("MuscleSkeletal","Skeletal_Muscle",colnames(data))


data2<-data[,3:dim(data)[2]]#Cut the matrix before ordering it by columns (to avoid gene_id and Description columns)
data3<-data2[ , order(names(data2))] #Order the columns by column names
data4<-cbind(data[1],data[2],data3) #Add first two columns back to the matrix



# 3. Importing the gene list and joining them with the initial input
#####################################################################

ensembl <- read.table(input2, sep="\t",header=TRUE)#import the ensembl file that contains ENSEMBL ID and matching GeneNames
joined<- plyr::join(data4, ensembl, by="gene_id") #use join function to add a column of gene names/Class corresponding to the ENSEMBL ID to the last column
joined2<- na.omit(joined, cols="Class") #Remove the rows that contain non-matching genes
#joined2 is tha expression file for RNA Modification Machinery Proteins (RMMs)
replaced <- joined2 %>% dplyr::select(Class, everything()) #place the last column to first column
replaced2 <- replaced %>% dplyr::select(Symbol, everything())#place the last column to first column
replaced2$Description <- NULL#Remove the redundant column (Gene names)
replaced3<-replaced2[order(replaced2$Class),] #Sort by Class


#Removing the Cell data because they introduce bias to tissue-wide expression analysis
replaced3[,"Cells_EBV"]<-NULL #Remove Cells 
replaced3[,"Cells_FIB"]<-NULL #Remove Cells
replaced3[,"WholeBlood"]<-NULL #Remove Whole Blood
write.table(replaced3,file="RMLP.GTEX.TPM.tsv",quote=FALSE, sep="\t",row.names=FALSE) #Export TPM for specific genes

#Take mean of tissues with multiple parts (Brain, Cerebellum etc.)
replaced3$Brain <- rowMeans(replaced3[,grep("^Brain_",colnames(replaced3))],na.rm = TRUE) #Mean of the columns that have matching name
replaced3[,grep("^Brain_",colnames(replaced3))]<- NULL #remove those columns after adding mean column
replaced3$Cerebellum <- rowMeans(replaced3[,grep("^Cerebellum_",colnames(replaced3))],na.rm = TRUE)
replaced3[,grep("^Cerebellum_",colnames(replaced3))]<- NULL
replaced3$Adipose <- rowMeans(replaced3[,grep("^Adipose_",colnames(replaced3))],na.rm = TRUE)
replaced3[,grep("^Adipose_",colnames(replaced3))]<- NULL
replaced3$Artery <- rowMeans(replaced3[,grep("^Artery_",colnames(replaced3))],na.rm = TRUE)
replaced3[,grep("^Artery_",colnames(replaced3))]<- NULL
replaced3$Cervix <- rowMeans(replaced3[,grep("^Cervix_",colnames(replaced3))],na.rm = TRUE)
replaced3[,grep("^Cervix_",colnames(replaced3))]<- NULL
replaced3$Colon <- rowMeans(replaced3[,grep("^Colon_",colnames(replaced3))],na.rm = TRUE)
replaced3[,grep("^Colon_",colnames(replaced3))]<- NULL
replaced3$Esophagus <- rowMeans(replaced3[,grep("^Esophagus_",colnames(replaced3))],na.rm = TRUE)
replaced3[,grep("^Esophagus_",colnames(replaced3))]<- NULL
replaced3$Heart <- rowMeans(replaced3[,grep("^Heart_",colnames(replaced3))],na.rm = TRUE)
replaced3[,grep("^Heart_",colnames(replaced3))]<- NULL
replaced3$Skin <- rowMeans(replaced3[,grep("^Skin_",colnames(replaced3))],na.rm = TRUE)
replaced3[,grep("^Skin_",colnames(replaced3))]<- NULL

#Order the columns by alphabetical order again
replaced4<-replaced3[,4:dim(replaced3)[2]]#Cut the matrix before ordering it by columns (to avoid gene_id and Description columns)
replaced5<-replaced4[ , order(names(replaced4))] #Order the columns by column names
final<-cbind(replaced3[1],replaced3[2],replaced3[3], replaced5) #Add first two columns back to the matrix
write.table(final,file="RMLP.GTEX.TissueAveraged.TPM.tsv",quote=FALSE, sep="\t",row.names=FALSE) #Export TPM for specific genes in averaged tissues
