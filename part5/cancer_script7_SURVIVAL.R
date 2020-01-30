#####SURVIVAL ANALYSIS
###############################
#Rscript cancer_script7_SURVIVAL.R TCGA_GTEX_FINAL.log2.without3cancer.tsv TCGA_survival_data
###Data for survival plots
library(tidyr)
args <- commandArgs(trailingOnly = TRUE) #Argument for first input
input1 <- args[1]#1st variable 
logfile<-read.table(input1, header=T, sep="\t") #TCGA_GTEX_FINAL.log2.without3cancer.tsv
surv<- logfile[with(logfile, logfile$sampletype %in% "Tumor"),]


###Reformat the logfile into high and low values (Expression-wise)
#High = Higher than mean
#Low = Lower than mean
array<- surv
array2<- seq(ncol(array)-4)
for (t in unique(surv$Type)){
  arraysub<-surv[surv$Type==t,-c(1:4)]
  colmean<- colMeans(arraysub)
  for (r in c(1:nrow(arraysub))){
    cond <- as.numeric(arraysub[r,]) > colmean
    cond<- gsub("TRUE", "HIGH", cond)
    cond<- gsub("FALSE", "LOW", cond)
    array2<- rbind(array2, cond)
}
}
array2 <- array2[-1,]
colnames(array2)<-colnames(surv[,-c(1:4)])
row.names(array2)<-NULL
array3<- cbind(surv[,c(1:4)],array2)



args <- commandArgs(trailingOnly = TRUE) #Argument for second input
input2 <- args[2]#2nd variable 
survivaldata<-read.delim(input2) #TCGA_survival_data


library(plyr)
joined<- join(array3, survivaldata, by="sample") #use join function to add a column of gene names/Class corresponding to the ENSEMBL ID to the last column
joined2<- na.omit(joined, cols="PFI.time") #Remove the rows that contain non-matching genes
library(dplyr) #load the package for data manipulation
#Testing
survivalinput <- joined %>% dplyr::select(OS,OS.time,DSS,DSS.time,DFI,DFI.time,PFI,PFI.time, everything()) #place the last column to first column


###plot survival curves
library("survminer")
require("survival")
for (t in unique(survivalinput$Type)){
  arraysub<-survivalinput[survivalinput$Type==t,]
  for (c in colnames(arraysub[,-c(1:12)])){
  formula<-paste("Surv(OS.time, OS) ~",c)
  form<- as.formula(formula)
  fit <- surv_fit(form, data = arraysub)
  pdf(file=paste(c,t,"survival.pdf",sep="."),height=5,width=5,onefile=FALSE)
  try(print(ggsurvplot(
  fit, 
  data = arraysub, 
  size = 1,                 # change line size
  palette = 
    c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  ggtheme = theme_bw(),     # Change ggplot2 theme
  xlab = paste(t,"Time in days",sep="-")
)))
  dev.off()
}
}




#Extract survival curve p value
library("survminer")
require("survival")
my_function <- function(x){  #We create a function 
  form <- paste("Surv(OS.time, OS) ~",x) #We state the formula for the p value calculation
  form <- as.formula(form) #We translate it into a formula
  fit <- surv_fit(form, data = survivalinput, group.by="Type") #Survival fits grouped by Type
  return(surv_pvalue(fit)) #P value
}
output <- list()
output <- sapply(colnames(survivalinput[,-c(1:12)]), my_function)
for (i in c(1:ncol(output))){
  for (j in c(1:nrow(output))){
  output[j,i]<- output[j,i][[1]][,2]
}
}


output2<- as.data.frame(output)
output2<- as.data.frame(unlist(output2))
output2$Gene<- gsub("\\..*","",rownames(output2))
rownames(output2)<- gsub(":.*","",rownames(output2))
output2$Type<-rownames(output2)
output2$Type<- gsub(".*\\.","",output2$Type)
rownames(output2)<-NULL
output2 <- output2 %>% dplyr::select(Gene,Type, everything()) #place the last column to first column
colnames(output2)<- c("Gene", "Type", "value")
write.table(output, file="survival_pvalues.tsv", quote=FALSE, sep="\t")
output3<-output2
output3$value<- 1/output3$value
output4<-output3
output4$value<- log(output4$value+1)
write.table(output4, file="survival_transformed_values.tsv",row.names=FALSE, quote=FALSE, sep="\t")
output5<- dcast(output4, Gene~Type)
write.table(output5, file="survival_transformed_values_reformatted.tsv",row.names=FALSE, quote=FALSE, sep="\t")


output4<- read.delim("survival_transformed_values.tsv")
library(tidyverse)
output_sign<-output4 %>% filter(value > 4.618)
output_sign$value <- factor(output_sign$value, levels = unique(output_sign$value))
output_sign2<-data.frame(table(output_sign[,1]))
output_sign3<-output_sign2 %>% filter(Freq > 2)
output_sign4 <- output_sign3[order(output_sign3$Freq),] 
output_sign4$Freq <- factor(output_sign4$Freq, levels = unique(output_sign4$Freq))
output_sign4$Var1 <- factor(output_sign4$Var1, levels = output_sign4$Var1[order(output_sign4$Freq)])




library(ggplot2)
pdf("survival_count_barplot_0.01.pdf",height=5,width=5)
p <- ggplot(output_sign4, aes(x=Var1,y=Freq)) + 
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
         plot.margin=unit(c(1,1,1,1),"line"))
print(p) 
dev.off()






##Box plots of transformed survival pvalues
library(ggplot2)
pdf("boxplot.survpvalueGene.pdf",height=10,width=20)
print(ggplot(output4, aes(x = Gene, y = value)) +
    geom_boxplot(fill="#8c96c6", color="darkred")+
    theme(axis.text.x = element_text(angle = 90))+
    geom_text(data = subset(output4,value > 6.91),
      mapping = aes(label = Type),vjust = 0, nudge_y = 0.05)+
    scale_y_continuous(limits=c(0,60), expand = c(0, 0))
    #annotate("segment", x =79, xend = 79, y = 23, yend = 25, arrow = arrow()) +
    #annotate("text", x = 79, y = 22, label = "LGG(NSUN7) at >50"))
    )
dev.off()





#heatmap of transformed survival pvalues
library(ComplexHeatmap)
library(circlize)
#Type<- as.data.frame(meanfile$Type)
#colnames(Type)<- c("Type")
rownames(output5) <- output5[,1]
output6<- output5
output6[,1]<-NULL
pdf("heatmap_survivaltransformedpvalues_2.pdf",height=20,width=15)
Heatmap(output6, name = "Transformed log(P-value)", 
        col = colorRamp2(c(0,2,4.6,7.5,10), c("#fcf9ec","#fcf9ec","#ffe2e2","#ff8260", "#900048"),space = "RGB"),
        cluster_columns = TRUE,
        column_title = "Cancer Tissues", 
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        column_names_gp = gpar(fontsize = 13, fontface = "bold"),
        row_title = "RNA Modification Enzymes", row_title_rot = 90,
        row_title_gp = gpar(fontsize = 8, fontface = "bold"),
        cluster_rows = TRUE,
        row_dend_width = unit(4, "cm"),
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 10), #row names size
        #column_order = 1:dim(data4)[2],#Keep the column order, make clustering FALSE for this
        row_dend_side = "right", #Dendogram on the right side
        #row_order = 1:dim(data4)[1], #Keep the row order, make clustering FALSE for this
        show_column_dend = TRUE, #
        column_dend_side = "top",
        column_dend_height = unit(2.5, "cm"),
        column_names_side = "bottom")
dev.off()

