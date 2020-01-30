###DATA SIMILARITY 
#Library required
#library(MASS)
#library(reshape)
#library(scales)
#library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
input1 <- args[1] #1st variable
gtex<- read.delim(input1) #RMLP.GTEX.TissueAveraged.TPM.tsv
rownames(gtex)<- gtex$Symbol
gtex<- gtex[,-c(1:3)]
input2 <- args[2] #1st variable
encode<- read.delim(input2)#RMLP.encode.TPM.brainav.tsv
rownames(encode)<- encode$Symbol
encode<- encode[,-c(1:4)]


rd1<-gtex
ua <-encode 

##
rd1_cn = colnames(rd1) 
ua_cn = colnames(ua) 
common_t = merge(rd1_cn,ua_cn,by.x=1,by.y=1) 
common_t = as.character(common_t[,1]) 
rd1 = rd1[,common_t] 
ua = ua[,common_t] 
rd1_d = rownames(rd1) 
ua_d = rownames(ua) 
common_d = merge(rd1_d,ua_d,by.x=1,by.y=1) 
common_d = as.character(common_d[,1]) 
rd = rd1[common_d,] 
ua = ua[common_d,] 


gtexM<-rd
encodeM<-ua


####CALCULATE RESIDUALS FOR COMPARISON OF DATASETS
library(MASS)#for RLM function
gtexl<-log(gtexM+1)
gtexres<-gtexl
res_vec <- vector() #res_vec file is a vector file that will contain all the residuals (tissue.vs.all) for all of the genes in all tissues
	for (i in 1:dim(gtexl)[2]) {
		res <- rlm(gtexl[,i] ~0+ as.vector(rowMeans(gtexl))) #rlm function is a robust linear model which is less affected by outliers. Starts from 0,0 point
		gtexres[,i]<-res$residuals
		#res_vec = c(res$residuals,res_vec) #this contains residuals for every gene in every tissuevsall combination
}

encodel<-log(encodeM+1)
encoderes<-encodel
res_vec <- vector() #res_vec file is a vector file that will contain all the residuals (tissue.vs.all) for all of the genes in all tissues
	for (i in 1:dim(encodel)[2]) {
		res <- rlm(encodel[,i] ~0+ as.vector(rowMeans(encodel))) #rlm function is a robust linear model which is less affected by outliers. Starts from 0,0 point
		encoderes[,i]<-res$residuals
		#res_vec = c(res$residuals,res_vec) #this contains residuals for every gene in every tissuevsall combination
	}


#cormat<- cor(gtexres,encoderes)
library(reshape)
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


library(scales)
library(ggplot2)

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

