####################################################################
## Tumor and Normal Median Values and Heatmap
## 2020, Oguzhan Begik written for Begik et al, 2020 Genome Biology
###################################################################

##Tumor vs Normal Median Heatmap
	#Required libraries
		#library(ggplot2)
		#library(ggpubr) 
		#library(EnvStats)
		#library(ComplexHeatmap)
		#library(circlize)
		#library(plyr)


#LIBRARIES NEEDED
library(ggplot2)
library(ggpubr) 
library(EnvStats)
library(ComplexHeatmap)
library(circlize)
library(plyr)

args <- commandArgs(trailingOnly = TRUE) #Argument for first input
input1 <- args[1]#1st variable 





#INPUT
logfile<-read.table(input1, header=T, sep="\t") #TCGA_GTEX_FINAL.log2.without3cancer.tsv



#Define a function to extract medians for each tissue type
#Function to calculate statistics of the data: 
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(median = median(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- plyr::rename(data_sum, c("median" = varname))
 return(data_sum)
}

#Use log(TPM+1) file
dat<-logfile
###For loop for the analysis
test_file <- vector()
dat2 <- data_summary(dat, varname="ADAD1", #Just to get the first two columns #Replace the varname with any column name
                    groupnames=c("sampletype", "Type"))
test_file <- noquote(cbind(as.character(dat2[,1]),as.character(dat2[,2])))
for (i in c(5:dim(dat)[2])){
vname<- noquote(paste(colnames(dat)[i]))
dat2 <- data_summary(dat, varname=vname, 
                    groupnames=c("sampletype", "Type"))
test_file <- cbind(test_file, dat2[,3])
}
header <- colnames(dat[,-c(2,4)])
final <- rbind(header,test_file)
write.table(final, file= "medianlog_tumor&normal.tpm.tsv", quote=FALSE,row.names=F, col.names=F,sep="\t")
medianfile<-read.table("medianlog_tumor&normal.tpm.tsv", header=T, sep="\t")




###Heatmap with median log(TPM+1) values Tumor & Normal
medianfile2<-medianfile[,-c(1,2)]#Remove first two columns
####HEATMAP WITH COMPLEX HEATMAP PACKAGE
#https://bioconductor.org/packages/release/bioc/vignettes/ComplexHeatmap/inst/doc/s2.single_heatmap.html
Type<- as.data.frame(medianfile$Type)
colnames(Type)<- c("Type")
ha = rowAnnotation(df = Type, col = list(Type = c("Normal" = "#66CC99", "Tumor" = "#FF6961")),
    width = unit(0.2, "cm"))
medianfile3<-scale(medianfile2) #Normalize by column (gene)
pdf("heatmap.tumorvsnormal_median.logtpm.pdf",height=20,width=25)
Heatmap(medianfile3, name = "z-scale log(TPM+1)", 
        col = colorRamp2(c(-2.5,-1,0,1,2.5), c("#2c7bb6","#abd9e9","floralwhite","#fdae61", "#d7191c"),space = "RGB"),
        #cluster_rows = TRUE, 
        cluster_columns = TRUE,
        column_title = "mRNA Expression", 
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        column_names_gp = gpar(fontsize = 13, fontface = "bold"),
        row_title = "RNA Modification Enzymes", row_title_rot = 90,
        row_title_gp = gpar(fontsize = 8, fontface = "bold"),
        cluster_rows = FALSE,
        show_row_names = FALSE,
        row_names_gp = gpar(fontsize = 10), #row names size
        #column_order = 1:dim(data4)[2],#Keep the column order, make clustering FALSE for this
        row_dend_side = "right", #Dendogram on the right side
        #row_order = 1:dim(data4)[1], #Keep the row order, make clustering FALSE for this
        show_column_dend = TRUE, #
        column_dend_side = "top",
        column_dend_height = unit(2.5, "cm"),
        column_names_side = "bottom",
        split = medianfile$detailedcategory, #Splitting by Class
        gap = unit(1, "mm"), #Gap
)+ha
dev.off()