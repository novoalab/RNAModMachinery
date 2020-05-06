###################################################################
## Tissue-wide expression plots using GTEx values
## 2020, Oguzhan Begik written for Begik et al, 2020 Genome Biology
###################################################################

##Required libraries:
#MASS
#ggplot2
#ggrepel
#ComplexHeatmap
#circlize
#ggplot2
#grid
#gridExtra



#how to use the script
#Rscript encode_tissuewide.R RMLP.encode.TPM.brainav.tsv

# 1. Arguments introduced for the execution
########################################

args <- commandArgs(trailingOnly = TRUE) #Argument for first input
input1 <- args[1]#1st variable 


# 2. Importing and manipulating the data
########################################

library(MASS)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(grid)
library(gridExtra)



##Tissue Wide Expression Plots for encode
data<- read.delim(input1)#"RMLP.encode.TPM.brainav.tsv"
scatter<- melt(data, c("Symbol","Class","gene_id")) #data for scatterplot
colnames(scatter)<- c("Symbol","Class","gene_id","Tissue", "value") #change column names
scatter$value<- log(scatter$value+1)#Take log of the file (with a pseudocount)


# 3. SCATTER PLOTS
########################################


#Calculate THRESHOLD for the specificity
genemean<- aggregate(scatter[, 5], list(scatter$Symbol), mean) #Row mean grouped by Gene
colnames(genemean)<- c("Symbol", "genemean") 
scatter2<- plyr::join(scatter, genemean, by="Symbol") #Add Rowmeans to the original data
scatter2$abs<- scatter2$value- scatter2$genemean #absolute distance of gene's expresion in tissue A from mean expression of this gene ins all tissues
res_vec <- vector() #res_vec file is a vector file that will contain all the residuals (tissue.vs.all) for all of the genes in all tissues
for (tissue in unique(scatter2$Tissue)){ #for every single tissue
	subset <- scatter2[with(scatter2, scatter2$Tissue %in% tissue),] #extract the data for a specific tissue
	res<- rlm(subset$value ~0 + subset$genemean) #linear model for that tissue 
	res_vec= c(res$residuals,res_vec)#this contains residuals for every gene in every tissuevsall combination
}
threshold <- 2.5*sd(res_vec) #The threshold is 2.5 times the standard deviation of all the residuals


##Seperate plots and calculations for each tissue
specific_genes<-vector() 
for (tissue in unique(scatter2$Tissue)){ #for each tissue
subset <- scatter2[with(scatter2, scatter2$Tissue %in% tissue),] #extract the data for a specific tissue
res<- rlm(subset$value ~0 + subset$genemean)#linear model for that tissue 
subset$res<- res$residuals #add residual values to the matrix
subset$diff<- abs(subset$res)-threshold #difference between gene's residual and threshold
spec<-subset(subset, diff>0 & abs>0) #extract specific genes in each tissue
specific_genes<- rbind(spec,specific_genes) #add these genes to the initial data
pdf(file=paste(tissue,"specificityplot.pdf",sep="."),height=5,width=5)
print(ggplot(subset, aes(x=genemean, y=value,label=Symbol)) + 
	  scale_x_continuous(limits= c(0,5))+
	  scale_y_continuous(limits= c(0,5))+
	  geom_point(data=subset, col="black",size=0.5)+ #All data points will be black
	  geom_point(data=subset(subset, diff>0 & abs>0),col="red",size=2)+ #Except the specific genes
	  geom_text_repel(data=subset(subset, diff>0 & abs>0),segment.size  = 0.4,segment.color = "grey50",)+ #Add text to the specific genes
	  geom_smooth(method=rlm, formula = y ~0 + x, size=0.5,fullrange=TRUE)+ #abline will be from rlm function that passes through 0,0
	  xlab("mRNA mean abundance All Human Tissues")+
	  ylab(paste("mRNA mean abundance",tissue,sep=" "))+
	  theme(panel.background = element_blank(),
         panel.border=element_rect(fill=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background=element_blank(),
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         axis.ticks=element_line(colour="black"),
        plot.margin=unit(c(1,1,1,1),"line")))
dev.off()
}
write.table(specific_genes, file="specific_genes_encode.tsv",quote=FALSE, row.names=FALSE,sep="\t")




# 4. HEATMAP
########################################
rownames(data)<- data[,1] #assign gene names as rownames  
data2<- data[,-c(1:3)] #remove the first three columns for the heatmap
data3<- log(data2+1)#Take log of the file (with a pseudocount)
data4 <- t(scale(t(data3)))#Normalize by row (by gene)

pdf("heatmap.encode.logzscaled.pdf",height=12,width=8)
Heatmap(data4, name = "z-scale log(TPM)", 
	#col = colorRamp2(c(-3,0,4), c("cadetblue3","floralwhite", "maroon4"),space = "RGB"), 
    #cluster_rows = TRUE, 
    col = colorRamp2(c(-3,-1.5,0,1.5,3), c("#2c7bb6","#abd9e9","floralwhite","#fdae61", "#d7191c"),space = "RGB"),
    cluster_columns = TRUE,
    column_title = "mRNA Expression of Human Tissues (encode)", 
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_names_gp = gpar(fontsize = 7, fontface = "bold"),
    row_title = "RNA Modification Enzymes", row_title_rot = 90,
    row_title_gp = gpar(fontsize = 8, fontface = "bold"),
    cluster_rows = TRUE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 5), #row names size
    #column_order = 1:dim(data4)[2],#Keep the column order, make clustering FALSE for this
    row_dend_side = "right", #Dendogram on the right side
    #row_order = 1:dim(data4)[1], #Keep the row order, make clustering FALSE for this
    show_column_dend = TRUE, #
    column_dend_side = "top",
    column_names_side = "bottom",
    split = data$Class, #Splitting by Class
    gap = unit(1, "mm"), #Gap
    )
dev.off()

# 4. PCA
########################################

##PCA Plotting 
data_pca<-prcomp(data4,center=TRUE) #PCA 
data_out <- as.data.frame(data_pca$x) #X table of PCA
data_out$Class<- data$Class
data_out$Symbol<- rownames(data)
data_out$Class <- factor(data_out$Class , levels = unique(data_out$Class))#keep the order of species
##Calculation for percentage of variance explained by each component
eigs <- data_pca$sdev^2#Calculate percentage for PC values
percentage<- round(eigs/sum(eigs)*100,2)#Calculate percentage for PC values
percentage <- paste( colnames(data_out), "(", paste( as.character(percentage), "%", ")", sep="") ) #Calculate percentage for PC values

## PLOT FOR X ALL LABELS
pdf("pca_encode.all.labels.pdf",height=5,width=7)
print(ggplot(data_out,aes(x=PC1,y=PC2,color=Class,label=Symbol))+
	scale_color_manual(values = c("#117A65","#D2B4DE","#F1948A","#B03A2E","#85C1E9","#17202A","#7B7D7D"))+
	geom_hline(yintercept = 0, lty = 2) +
	geom_vline(xintercept = 0, lty = 2)+
	geom_point(alpha = 0.8, size = 1.2)+
	geom_text(aes(label=Symbol),hjust=0, vjust=-1)+
	theme(panel.background = element_blank(),
		panel.border=element_rect(fill=NA),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		strip.background=element_blank(),
		axis.text.x=element_text(colour="black"),
		axis.text.y=element_text(colour="black"),
		axis.ticks=element_line(colour="black"),
		plot.margin=unit(c(1,1,1,1),"line"))+
	xlab(percentage[1]) + ylab(percentage[2]) #Labels containing percentages
	)
dev.off()


## PLOT FOR X SOME LABELS
pdf("pca_encode.pdf",height=5,width=7)
print(ggplot(data_out,aes(x=PC1,y=PC2,color=Class,label=Symbol))+
	scale_color_manual(values = c("#117A65","#D2B4DE","#F1948A","#B03A2E","#85C1E9","#17202A","#7B7D7D"))+
	geom_hline(yintercept = 0, lty = 2) +
	geom_vline(xintercept = 0, lty = 2)+
	geom_point(alpha = 0.8, size = 1.2)+
	geom_text(data=subset(data_out, PC1 < -3| PC2 > 2 | PC2 < -2),
		aes(label=Symbol),hjust=0, vjust=-1)+
	theme(panel.background = element_blank(),
		panel.border=element_rect(fill=NA),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		strip.background=element_blank(),
		axis.text.x=element_text(colour="black"),
		axis.text.y=element_text(colour="black"),
		axis.ticks=element_line(colour="black"),
		plot.margin=unit(c(1,1,1,1),"line"))+
	xlab(percentage[1]) + ylab(percentage[2]) #Labels containing percentages
	)
dev.off()



##PLOT FOR LOADINGS ALL LABELS
data_out_r <- as.data.frame(data_pca$rotation) #rotation data (loadings)
data_out_r$Symbol <- row.names(data_out_r) 
pdf(file="pca_encode.loadings.all.labels.pdf",height=5,width=5)
print(ggplot(data_out_r,aes(x=PC1,y=PC2,label=Symbol))+
	geom_point()+
	geom_hline(yintercept = 0, lty = 2) +
	geom_vline(xintercept = 0, lty = 2)+
	geom_point(alpha = 0.8, size = 2)+
	geom_text(aes(label=Symbol),hjust=0, vjust=-1)+
	theme(panel.background = element_blank(),
		panel.border=element_rect(fill=NA),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		strip.background=element_blank(),
		axis.text.x=element_text(colour="black"),
		axis.text.y=element_text(colour="black"),
		axis.ticks=element_line(colour="black"),
		plot.margin=unit(c(1,1,1,1),"line"))+
	xlab(percentage[1]) + ylab(percentage[2])
	)
dev.off()



##PLOT FOR LOADINGS SOME LABELS
data_out_r <- as.data.frame(data_pca$rotation) #rotation data (loadings)
data_out_r$Symbol <- row.names(data_out_r) 
pdf(file="pca_encode.loadings.some.labels.pdf",height=5,width=5)
print(ggplot(data_out_r,aes(x=PC1,y=PC2,label=Symbol))+
	geom_point()+
	geom_hline(yintercept = 0, lty = 2) +
	geom_vline(xintercept = 0, lty = 2)+
	geom_point(alpha = 0.8, size = 2)+
	geom_text(data=subset(data_out_r, PC1 > 0.2),
	aes(label=Symbol),hjust=0, vjust=-1)+
	theme(panel.background = element_blank(),
		panel.border=element_rect(fill=NA),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		strip.background=element_blank(),
		axis.text.x=element_text(colour="black"),
		axis.text.y=element_text(colour="black"),
		axis.ticks=element_line(colour="black"),
		plot.margin=unit(c(1,1,1,1),"line"))+
	xlab(percentage[1]) + ylab(percentage[2])
	)
dev.off()

