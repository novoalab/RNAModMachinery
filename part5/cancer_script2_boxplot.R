####################################################################
## BOXPLOT OF EACH RMP GENE (NORMAL VS TUMOR)
## 2020, Oguzhan Begik written for Begik et al, 2020 Genome Biology
###################################################################


#Required libraries
#library(ggplot2)
#library(ggpubr) 
#library(EnvStats)



#LIBRARIES NEEDED
library(ggplot2)
library(ggpubr) 
library(EnvStats)

#INPUT
args <- commandArgs(trailingOnly = TRUE) #Argument for first input
input1 <- args[1]#1st variable 
logfile<-read.table(input1, header=T, sep="\t") #TCGA_GTEX_FINAL.log2.without3cancer.tsv
#tpmfile<-read.table("TCGA_GTEX_FINAL.TPM.without3cancer.tsv", header=T, sep="\t")

dat<-logfile #Use logTPM file as dat file
#Print boxplot with outliers with multiple groups ##TPM
for (i in c(5:dim(dat)[2])) {
pdf(file=paste(colnames(dat)[i],"boxplot.logtpm.pdf",sep="."),height=5,width=11)
print(ggplot(data=dat, aes_string(x=colnames(dat)[1], y=noquote(paste(colnames(dat)[i])), fill=colnames(dat)[4])) +
  geom_boxplot(outlier.shape=16,aes(colour = sampletype))+#colored outliers
  scale_colour_manual(values=c("#66CC99","#FF6961"))+
  scale_fill_manual(values=c("#66CC99","#FF6961"))+
  stat_summary(fun.y = median, geom="point",colour="darkred", size=0.75, position = position_dodge(width = .75)) +
  theme(
  # Remove panel border
  panel.border = element_blank(),  
  # Remove panel grid lines
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Remove panel background
  panel.background = element_blank(),
  # Add axis line
  axis.line = element_line(colour = "grey"),
  axis.text.x = element_text(color="black", 
                           size=12, angle=45,hjust = 1),
  axis.text.y = element_text(color="black", 
                           size=12, angle=45)
  ))
dev.off()
}

