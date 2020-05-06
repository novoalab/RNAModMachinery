####################################################################
## Tumor and Normal Log2FC Calculation
## 2020, Oguzhan Begik written for Begik et al, 2020 Genome Biology
###################################################################
#Log2 FC calculation
#The expression data are first log2(TPM+1) transformed for differential analysis and 
#the log2FC is defined as median(Tumor) - median(Normal).
#Genes with higher |log2FC| values  than pre-set thresholds (1.5) 
#are considered differentially expressed genes.


  #Required libraries
    #library(ggplot2)
    #library(ggpubr) 
    #library(EnvStats)
    #library(ComplexHeatmap)
    #library(circlize)
    #library(plyr)
    #library(reshape2)
    #library(ggrepel) 


#Rscript cancer_script4_log2FC.R medianlog_tumor_normal.tpm.tsv RMLP.TcgaTargetGtex_rsem_gene_tpm_withheader.tsv

# 1. Arguments introduced for the execution
########################################

args <- commandArgs(trailingOnly = TRUE) #Argument for first input
input1 <- args[1]#1st variable 
input2 <- args[2]#2nd variable 




# 2.Importing the data and manipulating
########################################


#Load the library
    library(ggplot2)
    library(ggpubr) 
    library(EnvStats)
    library(ComplexHeatmap)
    library(circlize)
    library(plyr)
    library(reshape2)
    library(ggrepel) 



medianfile<-read.table(input1, header=T, sep="\t") #medianlog_tumor_normal.tpm.tsv
medianfileN <- medianfile[with(medianfile, medianfile$Type %in% "Normal"),]
medianfileT <- medianfile[with(medianfile, medianfile$Type %in% "Tumor"),]
log2FC<- medianfileT[,-c(1,2)]-medianfileN[,-c(1,2)]
log2FC_2<-cbind(medianfileN[,2],log2FC)
log2FC_3<-log2FC_2
rownames(log2FC_3)<-log2FC_2[,1]
log2FC_3<-log2FC_3[,-1]
log2FC_4<-t(log2FC_3)
log2FC_5<- melt(log2FC_4) #To export
write.table(log2FC_5,file="log2FC_TvsN.tsv",quote=FALSE,sep="\t",row.names=FALSE)  #Export melted file containing log(ratio)



## Labeling significant genes
log2fc_significant <- log2FC_5
colnames(log2fc_significant)<- c("Gene", "Tissue", "log2FC")
log2fc_significant2<-subset(log2fc_significant, log2FC> 1 | log2FC < -1)
write.table(log2fc_significant2,file="log2FC_TvsN_significants.tsv",quote=FALSE,sep="\t",row.names=FALSE)  #Export melted file containing log(ratio)

up_log2FC<-log2fc_significant2[(log2fc_significant2$log2FC >1),]
down_log2FC<-log2fc_significant2[(log2fc_significant2$log2FC < -1),]
up_log2FC2<-data.frame(table(up_log2FC[,1]))
down_log2FC2<-data.frame(table(down_log2FC[,1]))
up_log2FC3<-up_log2FC2[(up_log2FC2$Freq >1),]
down_log2FC3<-down_log2FC2[(down_log2FC2$Freq >1),]
up_log2FC4 <- up_log2FC3[order(up_log2FC3$Freq),] 
down_log2FC4 <- down_log2FC3[order(down_log2FC3$Freq),]
up_log2FC4$Var1 <- factor(up_log2FC4$Var1, levels = unique(up_log2FC4$Var1))
down_log2FC4$Var1 <- factor(down_log2FC4$Var1, levels = unique(down_log2FC4$Var1))










# 3.Barplot of UP-regulated and DOWN-regulated RMPs
####################################################


pdf("log2FC_UP_Barplot.pdf",height=10,width=5)
p <- ggplot(up_log2FC4, aes(x=Var1,y=Freq)) + 
       geom_bar(stat="identity",width = 0.75)+
       scale_fill_manual(values = c("#216583", "#f76262"))+
       coord_flip()+
       theme(panel.background = element_blank(),
         panel.border=element_rect(fill=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background=element_blank(),
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         axis.ticks=element_line(colour="black"),
         plot.margin=unit(c(1,1,1,1),"line"))+
         ylim(0,15)
print(p) 
dev.off()


pdf("log2FC_DOWN_Barplot.pdf",height=10,width=5)
p <- ggplot(down_log2FC4, aes(x=Var1,y=Freq)) + 
       geom_bar(stat="identity",width = 0.75)+
       scale_fill_manual(values = c("#216583", "#f76262"))+
       coord_flip()+
       theme(panel.background = element_blank(),
         panel.border=element_rect(fill=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background=element_blank(),
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         axis.ticks=element_line(colour="black"),
         plot.margin=unit(c(1,1,1,1),"line"))+
       ylim(0,15)
print(p) 
dev.off()








# 4.Barplot of Tissues with their significant gene counts
####################################################

log2fc_significant <- log2FC_5
colnames(log2fc_significant)<- c("Gene", "Tissue", "log2FC")
log2fc_significant2<-subset(log2fc_significant, log2FC> 1 | log2FC < -1)
write.table(log2fc_significant2,file="log2FC_TvsN_significants.tsv",quote=FALSE,sep="\t",row.names=FALSE)  #Export melted file containing log(ratio)

up_log2FC<-log2fc_significant2[(log2fc_significant2$log2FC >1),]
down_log2FC<-log2fc_significant2[(log2fc_significant2$log2FC < -1),]
up_log2FC2<-data.frame(table(up_log2FC[,2]))
down_log2FC2<-data.frame(table(down_log2FC[,2]))
up_log2FC3<-up_log2FC2[(up_log2FC2$Freq >1),]
down_log2FC3<-down_log2FC2[(down_log2FC2$Freq >1),]
up_log2FC4 <- up_log2FC3[order(up_log2FC3$Freq),] 
down_log2FC4 <- down_log2FC3[order(down_log2FC3$Freq),]
up_log2FC4$Var1 <- factor(up_log2FC4$Var1, levels = unique(up_log2FC4$Var1))
down_log2FC4$Var1 <- factor(down_log2FC4$Var1, levels = unique(down_log2FC4$Var1))

pdf("log2FC_UP_TissueBarplot.pdf",height=10,width=5)
p <- ggplot(up_log2FC4, aes(x=Var1,y=Freq)) + 
       geom_bar(stat="identity",width = 0.75)+
       scale_fill_manual(values = c("#216583", "#f76262"))+
       coord_flip()+
       theme(panel.background = element_blank(),
         panel.border=element_rect(fill=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background=element_blank(),
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         axis.ticks=element_line(colour="black"),
         plot.margin=unit(c(1,1,1,1),"line"))+
         ylim(0,65)
print(p) 
dev.off()

pdf("log2FC_DOWN_TissueBarplot.pdf",height=10,width=5)
p <- ggplot(down_log2FC4, aes(x=Var1,y=Freq)) + 
       geom_bar(stat="identity",width = 0.75)+
       scale_fill_manual(values = c("#216583", "#f76262"))+
       coord_flip()+
       theme(panel.background = element_blank(),
         panel.border=element_rect(fill=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background=element_blank(),
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         axis.ticks=element_line(colour="black"),
         plot.margin=unit(c(1,1,1,1),"line"))+
       ylim(0,65)
print(p) 
dev.off()










# 5.Heatmap
####################################################

####This is for the Class information in the heatmap
data<-read.delim(input2) #RMLP.TcgaTargetGtex_rsem_gene_tpm_withheader.tsv
data$sample <- gsub("\\..*","",data$sample) #remove everything after "." in ENSEMBL IDs
ensembl <- read.table("human_id_symbol_class.tsv", sep="\t",header=TRUE)#import the ensembl file that contains ENSEMBL ID and matching GeneNames
colnames(ensembl)<- gsub("gene_id","sample",colnames(ensembl)) 
joined<- join(data, ensembl, by="sample")


pdf("heatmap_log2FC.pdf",height=15,width=15)
heat<-Heatmap(log2FC_4, name = "log(ratio)", 
        col = colorRamp2(c(-1.5,-1,0,1,1.5), c("#2c7bb6","#abd9e9","floralwhite","#fdae61", "#d7191c"),space = "RGB"), 
        #cluster_rows = TRUE, 
        cluster_columns = TRUE,
        column_title = "log2FC Tumor vs Normal", 
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        column_names_gp = gpar(fontsize = 7, fontface = "bold"),
        row_title = "RNA Modification Enzymes", row_title_rot = 90,
        row_title_gp = gpar(fontsize = 8, fontface = "bold"),
        cluster_rows = TRUE,
        show_row_names = TRUE,
        #top_annotation = boxplotf,
        row_names_gp = gpar(fontsize = 8), #row names size
        #column_order = 1:dim(data4)[2],#Keep the column order, make clustering FALSE for this
        row_dend_side = "right", #Dendogram on the right side
        #row_order = 1:dim(data4)[1], #Keep the row order, make clustering FALSE for this
        show_column_dend = TRUE, #
        column_dend_side = "top",
        column_names_side = "bottom",
        #split = joined$Class, #Splitting by Class
        gap = unit(1, "mm"), #Gap
)
heat
#+ha_mix_right
dev.off()


###Boxplot alphabetical order
boxdata<-log2FC_5
boxdata<-boxdata[order(boxdata$value),] #Sort by any gene to see if there is any data point that has off pattern
pdf("boxplot_log2FC_alphabetical.pdf",height=10,width=20)
print(ggplot(boxdata, aes(x = Var2, y = value, label=Var1)) +
        geom_boxplot(fill='#dedede', color="#211717")+
		geom_text_repel(data = subset(boxdata,value > 1.5| value < -1.5),segment.size  = 0.4,segment.color = "grey50")+
		xlab("Tissues")+
  		ylab("Log2FC")+
		theme(panel.background = element_blank(),
         panel.border=element_rect(fill=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background=element_blank(),
         axis.text.x=element_text(colour="black",size = 30, angle = 45,hjust = 1),
         axis.text.y=element_text(colour="black",size=30),
         axis.ticks=element_line(colour="black"),
         axis.title.x=element_text(size=30,face="bold"),
         axis.title.y=element_text(size=30,face="bold"),
         plot.margin=unit(c(1,1,1,1),"line")))
dev.off()





# 5.Boxplot with the same order as Heatmap
####################################################


cols<- as.data.frame(colnames(log2FC_4)) 
cols$order<- rownames(cols)
order<- as.data.frame(column_order(heat))
colnames(order)<- c("order")
colsorder<- plyr::join(order,cols , by="order") 
colsorder$order<-rownames(colsorder)
colnames(colsorder)<- c("order","Tissue")
log2FC_5$Var2 <- factor(log2FC_5$Var2, levels = colsorder$Tissue)
##Boxplot of dysregulation scores
pdf("boxplot_log2FC_heatmaporder.pdf",height=10,width=20)
print(ggplot(log2FC_5, aes(x = Var2, y = value, label=Var1)) +
        geom_boxplot(fill='#dedede', color="#211717")+
        #geom_point(data=subset(dysregulation_scores, diff>0 & res>0),col="#d7191c",size=3)+ #Except the dysregulated genes
        #geom_point(data=subset(dysregulation_scores, diff>0 & res<0),col="#2c7bb6",size=3)+ #Except the dysregulated genes
		geom_text_repel(data = subset(boxdata,value > 1.5| value < -1.5),segment.size  = 0.4,segment.color = "grey50")+
		xlab("Tissues")+
  		ylab("Log2FC")+
		theme(panel.background = element_blank(),
         panel.border=element_rect(fill=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background=element_blank(),
         axis.text.x=element_text(colour="black",size = 30, angle = 45,hjust = 1),
         axis.text.y=element_text(colour="black",size=30),
         axis.ticks=element_line(colour="black"),
         axis.title.x=element_text(size=30,face="bold"),
         axis.title.y=element_text(size=30,face="bold"),
         plot.margin=unit(c(1,1,1,1),"line")))
dev.off()



