###################################################################
## Expression profile similarity between datasets
## 2020, Oguzhan Begik written for Begik et al, 2020 Genome Biology
###################################################################

###Requirements
#libraries: 
#dplyr
#MASS
#reshape
#scales
#ggplot2


#how to use the script
#Rscript gtex_tissuewide.R RMLP.GTEX.TissueAveraged.TPM.tsv

# 1. Arguments introduced for the execution
########################################

args <- commandArgs(trailingOnly = TRUE) #Argument for first input
input1 <- args[1]#1st variable 
input2 <- args[2] #2nd variable



# 2. Importing and manipulating the data
########################################

#Load the library
library(dplyr)
library(reshape)
library(MASS)#for RLM function
library(scales)
library(ggplot2)




gtex_input<- read.delim(input1) #RMLP.GTEX.TissueAveraged.TPM.tsv
rownames(gtex_input)<- gtex_input$Symbol
gtex_input<- gtex_input[,-c(1:3)]
encode_input<- read.delim(input2)#RMLP.encode.TPM.brainav.tsv
rownames(encode_input)<- encode_input$Symbol
encode_input<- encode_input[,-c(1:4)]


gtex_cn = colnames(gtex_input) 
encode_cn = colnames(encode_input)
common_t = merge(gtex_cn,encode_cn,by.x=1,by.y=1) 
common_t = as.character(common_t[,1]) 
gtex = gtex[,common_t] 
encode = encode[,common_t] 
gtex_d = rownames(gtex) 
encode_d = rownames(encode) 
common_d = merge(gtex_d,encode_d,by.x=1,by.y=1) 
common_d = as.character(common_d[,1]) 
gtex = gtex[common_d,] 
encode = encode[common_d,] 


# 3. CALCULATE RESIDUALS FOR COMPARISON OF DATASETS
########################################
gtexl<-log(gtex+1)
gtexres<-gtexl
res_vec <- vector() #res_vec file is a vector file that will contain all the residuals (tissue.vs.all) for all of the genes in all tissues
	for (i in 1:dim(gtexl)[2]) {
		res <- rlm(gtexl[,i] ~0+ as.vector(rowMeans(gtexl))) #rlm function is a robust linear model which is less affected by outliers. Starts from 0,0 point
		gtexres[,i]<-res$residuals
		#res_vec = c(res$residuals,res_vec) #this contains residuals for every gene in every tissuevsall combination
}

encodel<-log(encode+1)
encoderes<-encodel
res_vec <- vector() #res_vec file is a vector file that will contain all the residuals (tissue.vs.all) for all of the genes in all tissues
	for (i in 1:dim(encodel)[2]) {
		res <- rlm(encodel[,i] ~0+ as.vector(rowMeans(encodel))) #rlm function is a robust linear model which is less affected by outliers. Starts from 0,0 point
		encoderes[,i]<-res$residuals
		#res_vec = c(res$residuals,res_vec) #this contains residuals for every gene in every tissuevsall combination
	}


#cormat<- cor(gtexres,encoderes)
cormat <- round(cor(gtexres,encoderes),2)
melted_cormat <- melt(cormat)



# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

upper_tri <- get_upper_tri(cormat)
lower_tri<-get_lower_tri(cormat)
melted_cormat <- melt(lower_tri, na.rm = TRUE)

# 3. Plotting the correlation matrix
########################################
pdf("gtex_encode_correlation_matrix.pdf",height=5,width=5)
ggplot(data = melted_cormat, aes(X1, X2, fill = value))+
 xlab("GTEX")+
 ylab("ENCODE")+
 geom_tile(color = "white")+
 geom_text(aes(X1, X2, label = value), color = "black", size = 1.5) +
 scale_fill_gradient2(low = "#2c7bb6", high = "#d7191c", mid = "floralwhite", 
   midpoint = 0, limit = c(-0.5,0.5), space = "Lab", 
   name="Pearson\nCorrelation",na.value="white",oob=squish) +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 10, hjust = 1))+
 coord_fixed()
 dev.off()

