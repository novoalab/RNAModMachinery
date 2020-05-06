###################################################################
## ANALYSING THE SPERMATOGENESIS DATA (Green et al, 2019)
## 2020, Oguzhan Begik written for Begik et al, 2020 Genome Biology
###################################################################

###Requirements
#libraries: 
#dplyr
#plyr
#cluster
#ggplot2
#grid
#gridExtra
#ComplexHeatmap
#circlize
#reshape2
####How to use the script
#Rscript spermatogenesis.R <spermatogenesis.expression.data> <ensembl_file>
#Rscript spermatogenesis.R spermatogenesis_scRNA_averageexpression.tsv gene_hgnc_ensmus.tsv





# 1. Arguments introduced for the execution
########################################

args <- commandArgs(trailingOnly = TRUE) #Argument for first input
input1 <- args[1]#1st variable 
input2 <- args[2] #2nd variable



# 2. Importing and manipulating the data
########################################

#Load the library
library(dplyr)
library(plyr)
library(cluster) 
library(ggplot2)
library(grid)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(reshape2)




data<-read.delim(input1) #spermatogenesis_scRNA_averageexpression.tsv
ensembl <- read.table(input2, sep="\t",header=TRUE)#import the ensembl file that contains ENSEMBL ID and matching GeneNames
library(dplyr)
library(plyr)
joined<- join(data, ensembl, by="Gene") #use join function to add a column of gene names/Class corresponding to the ENSEMBL ID to the last column
joined2<-joined[order(joined$ENSEMBL),]
joined3<- na.omit(joined2, cols="ENSEMBL")
replaced <- joined3 %>% select(ENSEMBL, everything()) #place the last column to first column
replaced2 <- replaced %>% select(HGNC, everything())#place the last column to first column
HGNC<-as.character(replaced2[,1])
Spermatogonia <- replaced2[,c("GC1")]
Prelep_Spermatocyte<- rowMeans(replaced2[,c("GC2","GC3")])
Spermatocytes<- rowMeans(replaced2[,c("GC4","GC5","GC6","GC7","GC8")])
Spermatids<- rowMeans(replaced2[,c("GC9","GC10","GC11")])
Elongating_Spermatids<- replaced2[,c("GC12")]
clustered<-cbind(HGNC,Spermatogonia,Prelep_Spermatocyte,Spermatocytes,Spermatids,Elongating_Spermatids)
head(clustered)
clustered<-as.data.frame(clustered)
clustered[,-1] <- data.frame(sapply(clustered[,-1], function(x) as.numeric(as.character(x))))
row.names(clustered)<- clustered$HGNC
clustered<- clustered[,-1]
head(clustered)
clustered<-cbind(HGNC,Spermatogonia,Prelep_Spermatocyte,Spermatocytes,Spermatids,Elongating_Spermatids)
clustered<-as.data.frame(clustered)
clustered[,-1] <- data.frame(sapply(clustered[,-1], function(x) as.numeric(as.character(x)))) #make the values numeric
row.names(clustered)<- clustered$HGNC
data<- clustered[,-1]
head(data)




# 3. K means clustering
########################################

pdf("wss.pdf",height=5,width=5)
wss <- (nrow(data)-1)*sum(apply(data,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(data, 
  centers=i)$withinss)
  plot(1:15, wss, type="b", xlab="Number of Clusters",
  ylab="Within groups sum of squares")
  dev.off()


cluster<- as.vector(kmeans(data,4)$cluster)
data2<- cbind(cluster,data) #new table with clusters
data_pca<-prcomp(data2[,-1],center=TRUE) #PCA 
data_out <- as.data.frame(data_pca$x) #X table of PCA
data_out$cluster <- data2$cluster
data_out<-data_out[order(data_out$cluster),]
data_out$symbol<- rownames(data_out)
data_out$symbol <- factor(data_out$symbol , levels = unique(data_out$symbol))#keep the order of species
eigs <- data_pca$sdev^2#Calculate percentage for PC values
percentage<- round(eigs/sum(eigs)*100,2)#Calculate percentage for PC values
percentage <- paste( colnames(data_out), "(", paste( as.character(percentage), "%", ")", sep="") ) #Calculate percentage for PC values





# 4.  PCA Plotting
########################################

pdf("pca_spermatogenesis_clustered_less_label.pdf",height=5,width=5.2)
print(ggplot(data_out,aes(x=PC1,y=PC2,color=cluster,label=symbol))+
     geom_point()+
     geom_text(data=subset(data_out, PC1 < -0.5),
        aes(label=symbol),hjust=0, vjust=-1)+
     geom_hline(yintercept = 0, lty = 2) +
     geom_vline(xintercept = 0, lty = 2)+
     geom_point(alpha = 1, size = 4)+
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


 #PCA Plotting X
pdf("pca_spermatogenesis_clustered.pdf",height=5,width=5.2)
print(ggplot(data_out,aes(x=PC1,y=PC2,color=cluster,label=symbol))+
     geom_point()+
     geom_text(aes(label=symbol),hjust=0, vjust=-1)+
     geom_hline(yintercept = 0, lty = 2) +
     geom_vline(xintercept = 0, lty = 2)+
     geom_point(alpha = 1, size = 4)+
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






data_out_r <- as.data.frame(data_pca$rotation) #rotation data (loadings)
data_out_r$Stage <- row.names(data_out_r) 
data_out_r$Stage <- factor(data_out_r$Stage , levels = unique(data_out_r$Stage))#keep the order of species
 pdf(file="pca_spermatogenesis.loadings.pdf",height=5,width=7)
 print(ggplot(data_out_r,aes(x=PC1,y=PC2,label=Stage,color=Stage))+
     geom_point()+
     scale_color_manual(values=c("#fc5185","#ff9e74","#e3c4a8","#3fc1c9","#364f6b"))+
     geom_hline(yintercept = 0, lty = 2) +
     geom_vline(xintercept = 0, lty = 2)+
     geom_point(alpha = 1, size = 4)+
     theme(panel.background = element_blank(),
         panel.border=element_rect(fill=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background=element_blank(),
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         axis.ticks=element_line(colour="black"),
         plot.margin=unit(c(1,1,1,1),"line"))+
     xlab(percentage[1]) + ylab(percentage[2])+
     geom_text(hjust=0, vjust=-1)
     )
dev.off()




# 5. Heatmap
########################################


datah<- data
datah$symbol<- rownames(data)
datah2<- join(datah, data_out, by="symbol") #use join function to add a column of gene names/Class corresponding to the ENSEMBL ID to the last column
datah3<- datah2[,-c(7:11)]
datah4 <- datah3 %>% select(symbol,cluster, everything())
datah5<-datah4[order(datah4$cluster),]
datah6<-datah5[,-c(1:2)] 
 pdf("spermatogenesis.heatmap.pdf",height=12,width=8)
 Heatmap(datah6, name = "log(expression)", 
     col = colorRamp2(c(0,0.25,0.5,1,2), c("#2c7bb6","#abd9e9","floralwhite","#fdae61", "#d7191c"),space = "RGB"), 
     #cluster_rows = TRUE, 
     cluster_columns = FALSE,
     column_title = "mRNA Expression in Mouse Spermatogenesis", 
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
     split = datah5$cluster, #Splitting by Class
     gap = unit(1, "mm"), #Gap
     )
dev.off()

# 6. Violin Plot
########################################

datah7<- melt(datah5,c("symbol","cluster"))
pdf("spermatogenesis.violinplot.pdf",height=5,width=10)
print(ggplot(datah7, aes(variable, value)) + 
     geom_violin(aes(fill = variable),scale="width",draw_quantiles = c(0.5))+
    #stat_summary(fun.y=mean, geom="point", shape=19, size=1)+
     geom_jitter(height = 0, width = 0.1,size=0.2)+
     facet_grid(~factor(cluster))+
     scale_fill_manual(values=c("#fc5185","#ff9e74","#e3c4a8","#3fc1c9","#364f6b"))+
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


# 7. Barplots
########################################

data2$symbol<- rownames(data2)
data3<- melt(data2,c("symbol","cluster"))
pdf(file="spermatogenesis_barplot.pdf",height=30,width=30)
print(ggplot(data=data3, aes(x=variable, y=value, fill=variable))+
  geom_bar(stat= "identity",position=position_dodge())+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values=c("#fc5185","#ff9e74","#e3c4a8","#3fc1c9","#364f6b"))+
  theme(
  # Remove panel border
  panel.border = element_blank(),  
  # Remove panel grid lines
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Remove panel background
  panel.background = element_blank(),
  # Add axis line
  axis.line = element_line(colour = "grey"))+
        facet_wrap(.~symbol,scales="free")
        )
dev.off()

