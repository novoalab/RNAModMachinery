####################################################################
## Dysregulation Analysis
## 2020, Oguzhan Begik written for Begik et al, 2020 Genome Biology
###################################################################


  #Required libraries
  #Required libraries
    #library(ggplot2)
    #library(ggpubr) 
    #library(ComplexHeatmap)
    #library(circlize)
    #library(plyr)
    #library(reshape2)
    #library(ggrepel) 
    #library(MASS)#for RLM function




# 1. Arguments introduced for the execution
########################################

args <- commandArgs(trailingOnly = TRUE) #Argument for first input
input1 <- args[1]#1st variable 
input2 <- args[2] #2nd variable






# 2.Importing the data and manipulating
########################################


library(MASS)
library(ggplot2)
library(ggrepel) 
library(MASS)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(plyr)



median_all<- read.delim(input1) #all_genes_logmedian_scatter_format.tsv
median_all$ID <- gsub("\\.*","",median_all$ID )
median_all$Origin<- rep("ALL",nrow(median_all))


#Input for RMP
median_RMP<-read.delim(input2) #medianlog_tumor_normal.tpm.tsv
median_RMP_N <- median_RMP[with(median_RMP, median_RMP$Type %in% "Normal"),] #Normal Mean Data
median_RMP_T <- median_RMP[with(median_RMP, median_RMP$Type %in% "Tumor"),] #Tumor Mean Data
median_RMP_N2<- median_RMP_N[,-c(1:2)]
median_RMP_T2<- median_RMP_T[,-c(1:2)]
rownames(median_RMP_N2)<- median_RMP_N$detailedcategory
rownames(median_RMP_T2)<- median_RMP_T$detailedcategory
median_RMP_N3<-t(median_RMP_N2)#Transpose
median_RMP_T3<-t(median_RMP_T2)#Transpose
Nexp<- melt(median_RMP_N3) #Reformat the matrix
Nexp$comb <- paste(Nexp[,1],Nexp[,2]) #Combine Gene-Tissue pair
Texp<-melt(median_RMP_T3)#Reformat the matrix
Texp$comb <- paste(Texp[,1],Texp[,2])#Combine Gene-Tissue pair
scatter<- plyr::join(Nexp, Texp, by="comb") #Join Normal and Tumor pairs
scatter<- scatter[,c(2,1,3,7)]
colnames(scatter)<- c("Type", "ID", "Normal", "Tumor") #Rename the columns
scatter$Origin<- rep("RMP",nrow(scatter))
final_merged_forthreshold<- rbind(median_all,scatter)


####PLOTTING ONLY THE RMPS####
#Calculate threshold

specific_genes<-vector()  #Empty vector for dysregulated genes
dysregulation_scores<- vector() #empty vector for dysregulation scores
for (tissue in unique(final_merged_forthreshold$Type)){ #for every single tissue
  res_vec <- vector() #res_vec file is a vector file that will contain all the residuals (tumor.vs.normal) for all of the genes in all tissues
  	subset <- final_merged_forthreshold[with(final_merged_forthreshold, final_merged_forthreshold$Type %in% tissue),] #extract the data for a specific tissue
	res<- rlm(subset$Tumor ~0 + subset$Normal) #linear model for that tissue
	res_vec= c(res$residuals,res_vec)#this contains residuals for every gene in every tissuevsall combination
	threshold <- 2.5*sd(res_vec) #The threshold is 2.5 times the standard deviation of all the residuals
	##Seperate plots and calculations for each tissue
	#subset <- scatter[with(scatter, scatter$Tissue %in% tissue),]
	#res<- rlm(subset$Tumor ~0 + subset$Normal)
	subset$res<- res$residuals
	subset$diff<- abs(subset$res)-threshold
	subset2<-subset[with(subset, subset$Origin %in% "RMP"),]
	spec<-subset(subset2, diff>0) #extract dysregulated genes in each tissue
	specific_genes<- rbind(spec,specific_genes) #add these genes to the initial data
	dysregulation_scores<- rbind(dysregulation_scores,subset2) #add these genes to the initial data
pdf(file=paste(tissue,"tumorvsnormalplot.pdf",sep="."),height=5,width=5)
print(ggplot(subset2, aes(x=Normal, y=Tumor,label=ID)) + 
  ggtitle(tissue)+
  geom_point(data=subset2, col="black",size=0.5)+ #All data points will be black
  geom_point(data=subset(subset2, diff>0 & res>0),col="#d7191c",size=2)+ #Except the dysregulated genes
  geom_point(data=subset(subset2, diff>0 & res<0),col="#2c7bb6",size=2)+ #Except the dysregulated genes
  geom_text_repel(data=subset(subset2, diff>0 ),segment.size  = 0.4,segment.color = "grey50",size=5)+ #Add text to the dysregulated genes
  geom_smooth(method=rlm, formula = y ~0 + x, size=0.5,color="black",se=F)+ #abline will be from rlm function that passes through 0,0
  xlab("log(TPM) Normal")+
  ylab("log(TPM) Tumor")+
  stat_cor(method = "pearson")+
  theme(panel.background = element_blank(),
         panel.border=element_rect(fill=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background=element_blank(),
         axis.text.x=element_text(colour="black",size=17),
         axis.text.y=element_text(colour="black",size=17),
         axis.ticks=element_line(colour="black"),
         plot.margin=unit(c(1,1,1,1),"line")))
dev.off()
}
write.table(specific_genes, file="dysregulated_genes.tsv",quote=FALSE, row.names=FALSE,sep="\t")
write.table(dysregulation_scores, file="dysregulation_scores.tsv",quote=FALSE, row.names=FALSE,sep="\t")






######PLOTTING ALL THE GENES AND LABELLING RMPS

specific_genes<-vector()  #Empty vector for dysregulated genes
dysregulation_scores<- vector() #empty vector for dysregulation scores
for (tissue in unique(final_merged_forthreshold$Type)){ #for every single tissue
  res_vec <- vector() #res_vec file is a vector file that will contain all the residuals (tumor.vs.normal) for all of the genes in all tissues
  	subset <- final_merged_forthreshold[with(final_merged_forthreshold, final_merged_forthreshold$Type %in% tissue),] #extract the data for a specific tissue
	res<- rlm(subset$Tumor ~0 + subset$Normal) #linear model for that tissue
	res_vec= c(res$residuals,res_vec)#this contains residuals for every gene in every tissuevsall combination
	threshold <- 2.5*sd(res_vec) #The threshold is 2.5 times the standard deviation of all the residuals
	##Seperate plots and calculations for each tissue
	library(ggplot2)
	library(ggrepel) #to repel the labels overlapping
	#subset <- scatter[with(scatter, scatter$Tissue %in% tissue),]
	#res<- rlm(subset$Tumor ~0 + subset$Normal)
	subset$res<- res$residuals
	subset$diff<- abs(subset$res)-threshold
	subset2<-subset[with(subset, subset$Origin %in% "RMP"),]
	spec<-subset(subset2, diff>0) #extract dysregulated genes in each tissue
	specific_genes<- rbind(spec,specific_genes) #add these genes to the initial data
	dysregulation_scores<- rbind(dysregulation_scores,subset2) #add these genes to the initial data
pdf(file=paste(tissue,"tumorvsnormalplot_ALLGENES.pdf",sep="."),height=5,width=5)
print(ggplot(subset, aes(x=Normal, y=Tumor,label=ID)) + 
  ggtitle(tissue)+
  geom_point(data=subset, col="grey",alpha = 0.6,size=0.02)+ #All data points will be transparent
  geom_point(data=subset2, col="black",size=0.6)+ #All data points will be black
  geom_point(data=subset(subset2, diff>0 & res>0),col="#F72419",size=2)+ #Except the dysregulated genes
  geom_point(data=subset(subset2, diff>0 & res<0),col="#5584FF",size=2)+ #Except the dysregulated genes
  geom_text_repel(data=subset(subset2, diff>0 & res>0 ),segment.size  = 0.4,col="#D61107", box.padding = unit(1, "lines"),
  point.padding = unit(0.3, "lines"),segment.color = "grey50",size=6)+ #Add text to the dysregulated genes
  geom_text_repel(data=subset(subset2, diff>0 & res<0),segment.size  = 0.4,col="#0047FF", box.padding = unit(1, "lines"),	
  point.padding = unit(0.3, "lines"),segment.color = "grey50",size=6)+ #Add text to the dysregulated genes
  #geom_smooth(method=rlm, formula = y ~0 + x, size=0.5,color="black",se=F)+ #abline will be from rlm function that passes through 0,0
  xlab("log(TPM) Normal")+
  ylab("log(TPM) Tumor")+
  stat_cor(method = "pearson",size=5)+
  theme(panel.background = element_blank(),
         panel.border=element_rect(fill=NA),
         plot.title=element_text(hjust=0.5, size=20),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background=element_blank(),
         axis.text.x=element_text(colour="black",size=17),
         axis.title.x=element_text(size=20),
         axis.text.y=element_text(colour="black",size=17),
         axis.title.y=element_text(size=20),
         axis.ticks=element_line(colour="black"),
         plot.margin=unit(c(1,1,1,1),"line")))
dev.off()
}
write.table(specific_genes, file="dysregulated_genes.tsv",quote=FALSE, row.names=FALSE,sep="\t")
write.table(dysregulation_scores, file="dysregulation_scores.tsv",quote=FALSE, row.names=FALSE,sep="\t")







######Making barplot for the significant gene count########################
dys_sign <- specific_genes[,c(1,2,6,7)]
up_dys<- dys_sign[(dys_sign$res>0),]
down_dys<- dys_sign[(dys_sign$res<0),]

up_dys2<-data.frame(table(up_dys[,2]))
down_dys2<-data.frame(table(down_dys[,2]))

up_dys3<-up_dys2[(up_dys2$Freq >1),]
down_dys3<-down_dys2[(down_dys2$Freq >1),]

up_dys4 <- up_dys3[order(up_dys3$Freq),] 
down_dys4 <- down_dys3[order(down_dys3$Freq),]

up_dys4$Var1 <- factor(up_dys4$Var1, levels = unique(up_dys4$Var1))
down_dys4$Var1 <- factor(down_dys4$Var1, levels = unique(down_dys4$Var1))


pdf("Dys_UP_Barplot2.pdf",height=5,width=5)
p <- ggplot(up_dys4, aes(x=Var1,y=Freq)) + 
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

pdf("Dys_DOWN_Barplot2.pdf",height=5,width=5)
p <- ggplot(down_dys4, aes(x=Var1,y=Freq)) + 
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




###Heatmap and boxplot###########################

##For the heatmap class distinctions
#data<-read.delim("RMLP.TcgaTargetGtex_rsem_gene_tpm_withheader.tsv")#RMLP.TcgaTargetGtex_rsem_gene_tpm_withheader.tsv
#library(dplyr) #load the package for data manipulation
#library(plyr)
#data$sample <- gsub("\\..*","",data$sample) #remove everything after "." in ENSEMBL IDs
#ensembl <- read.table("human_id_symbol_class.tsv", sep="\t",header=TRUE)#import the ensembl file that contains ENSEMBL ID and matching GeneNames
#colnames(ensembl)<- gsub("gene_id","sample",colnames(ensembl)) 
#joined<- join(data, ensembl, by="sample")

#HEATMAP OF DYSREGULATION SCORES

dysregulation_scores2<- dysregulation_scores[,-c(3,4,5,7)]
heatdysr<- dcast(dysregulation_scores2, ID~Type)
rownames(heatdysr)<- heatdysr$Gene
heatdysr<- heatdysr[-1]
pdf("heatmap.dysregulationscores.pdf",height=15,width=15)
heat<- Heatmap(heatdysr, name = "Dysregulation Scores", 
        col = colorRamp2(c(-1.2,-1,0,1,1.2), c("#2c7bb6","#abd9e9","floralwhite","#fdae61", "#d7191c"),space = "RGB"), 
        #cluster_rows = TRUE, 
        cluster_columns = TRUE,
        column_title = "Dysregulation Scores", 
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        column_names_gp = gpar(fontsize = 10, fontface = "bold"),
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
        column_names_side = "bottom"
)
heat
dev.off()

#To make boxplot the same order as heatmap
cols<- as.data.frame(colnames(heatdysr)) 
cols$order<- rownames(cols)
order<- as.data.frame(column_order(heat))
colnames(order)<- c("order")
colsorder<- plyr::join(order,cols , by="order") 
colsorder$order<-rownames(colsorder)
colnames(colsorder)<- c("order","Tissue")
dysregulation_scores$Type <- factor(dysregulation_scores$Type, levels = colsorder$Tissue)
##Boxplot of dysregulation scores
pdf("boxplot.dysregulationscores.pdf",height=10,width=20)
print(ggplot(dysregulation_scores, aes(x = Type, y = res)) +
        geom_boxplot(fill='#dedede', color="#211717")+
        geom_point(data=subset(dysregulation_scores, diff>0 & res>0),col="#d7191c",size=3)+ #Except the dysregulated genes
        geom_point(data=subset(dysregulation_scores, diff>0 & res<0),col="#2c7bb6",size=3)+ #Except the dysregulated genes
    geom_text_repel(data=subset(dysregulation_scores, diff>0),
      mapping = aes(label = ID),vjust = 0, nudge_y = 0.05,
      segment.size  = 0.4,segment.color = "grey50",size=5)+ #Add text to the dysregulated genes
    xlab("Tissues")+
      ylab("Dysregulation Scores")+
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

##Boxplot of dysregulation scores NO LABELS
pdf("boxplot.dysregulationscores.nolabel.pdf",height=5,width=20)
print(ggplot(dysregulation_scores, aes(x = Type, y = res)) +
        geom_boxplot(fill='#dedede', color="#211717")+
        geom_point(data=subset(dysregulation_scores, diff>0 & res>0),col="#d7191c",size=2)+ #Except the dysregulated genes
        geom_point(data=subset(dysregulation_scores, diff>0 & res<0),col="#2c7bb6",size=2)+ #Except the dysregulated genes
        xlab("Tissues")+
        ylab("Dysregulation Scores")+
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




#####$$$$$$$$$$$$$HEATMAP WITH SCALED SCORES#############$$$$$$$$$$$$$$$$$$$
dysregulation_scores2<- dysregulation_scores[,-c(3,4,5,7)]
heatdysr<- dcast(dysregulation_scores2, ID~Type)
rownames(heatdysr)<- heatdysr$ID
heatdysr<- heatdysr[-1]
heatdysr<-scale(heatdysr,center = TRUE, scale = TRUE) 
pdf("heatmap.dysregulationscores_scaled_centeredbycolumn.pdf",height=15,width=15)
heat<- Heatmap(heatdysr, name = "z-scaled Dysregulation Scores", 
        col = colorRamp2(c(-2,-1,0,1,2), c("#2c7bb6","#abd9e9","floralwhite","#fdae61", "#d7191c"),space = "RGB"), 
        #cluster_rows = TRUE, 
        cluster_columns = TRUE,
        column_title = "z-scaled Dysregulation Scores", 
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        column_names_gp = gpar(fontsize = 10, fontface = "bold"),
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
)
heat
dev.off()

#To make boxplot the same order as heatmap
cols<- as.data.frame(colnames(heatdysr)) 
cols$order<- rownames(cols)
order<- as.data.frame(column_order(heat))
colnames(order)<- c("order")
colsorder<- plyr::join(order,cols , by="order") 
colsorder$order<-rownames(colsorder)
colnames(colsorder)<- c("order","Type")
dysregulation_scores$Type <- factor(dysregulation_scores$Type, levels = colsorder$Type)
##Boxplot of dysregulation scores
pdf("boxplot.dysregulationscores_scaledcenteredbycolumnorder.pdf",height=10,width=20)
print(ggplot(dysregulation_scores, aes(x = Type, y = res)) +
        geom_boxplot(fill='#dedede', color="#211717")+
        geom_point(data=subset(dysregulation_scores, diff>0 & res>0),col="#d7191c",size=3)+ #Except the dysregulated genes
        geom_point(data=subset(dysregulation_scores, diff>0 & res<0),col="#2c7bb6",size=3)+ #Except the dysregulated genes
    geom_text_repel(data=subset(dysregulation_scores, diff>0),
      mapping = aes(label = ID),vjust = 0, nudge_y = 0.05,
      segment.size  = 0.4,segment.color = "grey50",size=5)+ #Add text to the dysregulated genes
    xlab("Tissues")+
      ylab("Dysregulation Scores")+
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


##Boxplot of dysregulation scores NO LABELS
pdf("boxplot.dysregulationscores.scaledcenteredbycolumnorder.pdf",height=5,width=20)
print(ggplot(dysregulation_scores, aes(x = Type, y = res)) +
        geom_boxplot(fill='#dedede', color="#211717")+
        geom_point(data=subset(dysregulation_scores, diff>0 & res>0),col="#d7191c",size=2)+ #Except the dysregulated genes
        geom_point(data=subset(dysregulation_scores, diff>0 & res<0),col="#2c7bb6",size=2)+ #Except the dysregulated genes
        xlab("Tissues")+
        ylab("Dysregulation Scores")+
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




