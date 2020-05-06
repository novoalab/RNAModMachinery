###################################################################
## EXTRACTING SPECIFIC GENES FROM ENCODE (MOUSE) TPM File
## 2020, Oguzhan Begik written for Begik et al, 2020 Genome Biology
###################################################################

###Requirements
#libraries: 
#dplyr
#stringr
####How to use the script
#Rscript encode_manipulation.R mm65.long.gene.with.expr.cshl.tsv mouse_id_symbol_class.tsv




# 1. Arguments introduced for the execution
########################################

args <- commandArgs(trailingOnly = TRUE) #Argument for first input
input1 <- args[1]#1st variable 
input2 <- args[2] #2nd variable




# 2. Importing and manipulating the data
########################################

#Load the library
library(dplyr)
library(stringr)



data<- read.delim(input1)
colnames(data)<- gsub("\\.","",colnames(data))#Remove ... and . from colnames introduced by R
#Removing other expression values from each cell (reformatting)
data2 <- as.data.frame(data[,1])
for (i in 2:ncol(data)) {
	data2[paste(colnames(data[i]))] <- str_split_fixed(data[,i],":",3)[,1]	
}
data2$gene_id<-rownames(data) #add gene_ids from the data
#Leaving only TPM values

data3 <- data2 %>% select(gene_id, everything()) #Place it in the first column
names(data3)[2] <- "symbol" #rename second column
write.table(data3,file="encode_TPM.tsv",quote=FALSE, sep="\t") #Exporting TPM values 
#Replace column names with proper tissue names
colnames(data3)<- gsub("LID20728LID20729Adrenaladult8wkscelllongPolyA","Adrenal",colnames(data3))
colnames(data3)<- gsub("LID20730LID20731Duodenumadult8wkscelllongPolyA","Duodenum",colnames(data3))
colnames(data3)<- gsub("LID20732LID20733Stomachadult8wkscelllongPolyA","Stomach",colnames(data3))
colnames(data3)<- gsub("LID20819LID20820SmIntestineadult8wkscelllongPolyA","Small_Intestine",colnames(data3))
colnames(data3)<- gsub("LID20821LID20822Ovaryadult8wkscelllongPolyA","Ovary",colnames(data3))
colnames(data3)<- gsub("LID20868LID20869Testisadult8wkscelllongPolyA","Testis",colnames(data3))
colnames(data3)<- gsub("LID20870LID20871Heartadult8wkscelllongPolyA","Heart",colnames(data3))
colnames(data3)<- gsub("LID20872LID20873Kidneyadult8wkscelllongPolyA","Kidney",colnames(data3))
colnames(data3)<- gsub("LID20920LID20921Lungadult8wkscelllongPolyA","Lung",colnames(data3))
colnames(data3)<- gsub("LID20922LID20923Thymusadult8wkscelllongPolyA","Thymus",colnames(data3))
colnames(data3)<- gsub("LID20924LID20925MammaryGlandadult8wkscelllongPolyA","Mammary",colnames(data3))
colnames(data3)<- gsub("LID21038LID21039Spleenadult8wkscelllongPolyA","Spleen",colnames(data3))
colnames(data3)<- gsub("LID21040LID21041Colonadult8wkscelllongPolyA","Colon",colnames(data3))
colnames(data3)<- gsub("LID21042LID21043Liveradult8wkscelllongPolyA","Liver",colnames(data3))
colnames(data3)<- gsub("LID21179LID21180GenitalFatPadadult8wkscelllongPolyA","Genital_Fat",colnames(data3))
colnames(data3)<- gsub("LID21181LID21182SubcFatPadadult8wkscelllongPolyA","Subc_Fat",colnames(data3))
colnames(data3)<- gsub("LID21183LID21184LgIntestineadult8wkscelllongPolyA","Large_Intestine",colnames(data3))
colnames(data3)<- gsub("LID46946LID46947CNSE115celllongPolyA","CNS_E11.5",colnames(data3))
colnames(data3)<- gsub("LID46948LID46949CNSE14celllongPolyA","CNS_E14",colnames(data3))
colnames(data3)<- gsub("LID46950LID46951CNSE18celllongPolyA","CNS_E18",colnames(data3))
colnames(data3)<- gsub("LID46983LID46984Placentaadult8wkscelllongPolyA","Placenta",colnames(data3))
colnames(data3)<- gsub("LID46985LID46986LimbE145celllongPolyA","Limb_E14.5",colnames(data3))
colnames(data3)<- gsub("LID46987LID46988WholeBrainE145celllongPolyA","Whole_Brain_E14.5",colnames(data3))
colnames(data3)<- gsub("LID47030LID47031Bladderadult8wkscelllongPolyA","Bladder",colnames(data3))
colnames(data3)<- gsub("LID47032LID47033Cortexadult8wkscelllongPolyA","Brain_COR",colnames(data3))
colnames(data3)<- gsub("LID47036LID47037Cerebellumadult8wkscelllongPolyA","Cerebellum",colnames(data3))
colnames(data3)<- gsub("LID47081LID47082FrontalLobeadult8wkscelllongPolyA","Brain_FRO",colnames(data3))
colnames(data3)<- gsub("LID47144LID47145LiverE14celllongPolyA","Liver_E14",colnames(data3))
colnames(data3)<- gsub("LID47146LID47147LiverE145celllongPolyA","Liver_E14.5",colnames(data3))
colnames(data3)<- gsub("LID47148LID47149LiverE18celllongPolyA","Liver_E18",colnames(data3))
colnames(data3)<- gsub("SID38132SID38133CH12adult8wkscelllongPolyA","CH12",colnames(data3))
colnames(data3)<- gsub("SID38134SID38135Liveradult8wkscelltotal","Liver_total",colnames(data3))

write.table(data3,file="encode_TPM_renamed.tsv",quote=FALSE, sep="\t") #Exporting Encode table with Renamed tissue names



# 3. Importing the gene list and joining them with the initial input
#####################################################################
ensembl <- read.table(input2, sep="\t",header=TRUE)#import the ensembl file that contains ENSEMBL ID and matching GeneNames
joined<- join(data3, ensembl, by="gene_id") #use join function to add a column of gene names/Class corresponding to the ENSEMBL ID to the last column
joined2<- na.omit(joined, cols="Class") #Remove the rows that contain non-matching genes
#joined2 is tha expression file for RNA Modification Related Proteins (RMRP)
replaced <- joined2 %>% select(Class, everything()) #place the last column to first column
replaced2 <- replaced %>% select(Symbol, everything())#place the last column to first column
replaced3<-replaced2[order(replaced2$Class),] #Sort by Class

#Removing the embryo samples, ch12, and total RNA seq liver
replaced3[,"Liver_total"]<-NULL  
replaced3[,"CH12"]<-NULL #
replaced3[,"Liver_E18"]<-NULL 
replaced3[,"Liver_E14.5"]<-NULL 
replaced3[,"Liver_E14"]<-NULL 
replaced3[,"Whole_Brain_E14.5"]<-NULL 
replaced3[,"Limb_E14.5"]<-NULL 
replaced3[,"CNS_E18"]<-NULL 
replaced3[,"CNS_E14"]<-NULL 
replaced3[,"CNS_E11.5"]<-NULL 

write.table(replaced3,file="RMLP.encode.TPM.tsv",quote=FALSE, sep="\t",row.names=FALSE)

#Take mean of tireplaced4<- ssues with multiple sections (Brain, Cerebellum etc.)
replaced4 <- data.frame(sapply(replaced3[,5:ncol(replaced3)], function(x) as.numeric(as.character(x))))
replaced5<- cbind(replaced3[,1:4],replaced4)
replaced5$Brain <- rowMeans(replaced5[,grep("^Brain_",colnames(replaced5))],na.rm = TRUE) #Mean of the columns that have matching name
replaced5[,grep("^Brain_",colnames(replaced5))]<- NULL #remove those columns after adding mean column
#Order the columns by alphabetical order again
replaced6<-replaced5[,5:dim(replaced5)[2]]#Cut the matrix before ordering it by columns (to avoid gene_id and Description columns)
replaced7<-replaced6[ , order(names(replaced6))] #Order the columns by column names
final<-cbind(replaced3[1:4], replaced7) #Add first four columns back to the matrix
final$symbol<-NULL
write.table(final,file="RMLP.encode.TPM.brainav.tsv",quote=FALSE, sep="\t",row.names=FALSE)
