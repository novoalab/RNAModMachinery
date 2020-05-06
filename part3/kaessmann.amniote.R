###################################################################
## Analysis of Kaessmann Amniote Data
## 2020, Oguzhan Begik written for Begik et al, 2020 Genome Biology
###################################################################

###Requirements
#library(plyr)
#library(dplyr) 
#library(ggplot2)
#library(grid)
#ibrary(gridExtra)
#library(MASS)



#How to run the script
#Rscript kaessman.amniote.R NormalizedRPKM_ConstitutiveExons_Amniote1to1Orthologues.txt human_id_symbol_class.tsv



# 1. Arguments introduced for the execution
########################################

args <- commandArgs(trailingOnly = TRUE) #Argument for first input
input1 <- args[1]#1st variable 
input2 <- args[2] #2nd variable



# 2. Importing and manipulating the data
########################################

#Load the library
library(plyr)
library(dplyr) 
library(ggplot2)
library(grid)
library(gridExtra)
library(MASS)


data<- read.delim(input1) #Amniote RPKM data
data[2:9]<-NULL  #Columns that contain other species ENSEMBL ID
ensembl <- read.table(input2, sep="\t",header=TRUE)#import the ensembl file that contains ENSEMBL ID and matching GeneNames
colnames(ensembl)[colnames(ensembl)=="gene_id"] <- "hsa" #Rename the gene_id column to hsa
joined<- join(data, ensembl, by="hsa") #use join function to add a column of gene names/Class corresponding to the ENSEMBL ID to the last column
joined2<- na.omit(joined, cols="Class") #Remove the rows that contain non-matching genes
#joined2 is tha expression file for RNA Modification Related Proteins (RMLP)
replaced <- joined2 %>% select(Class,Symbol, everything()) #place the last column to first column
replaced2<-replaced[order(replaced$Class),] #Sort by ClasS
replaced2$hsa<-NULL
replaced2$Class<-NULL
replaced3<-replaced2
rownames(replaced3)<- replaced3$Symbol
replaced3$Symbol<-NULL
replaced4<-t(replaced3) #Transpose the data
#Renames the species and tissues
rownames(replaced4)<- gsub("hsa.","Human_", rownames(replaced4))
rownames(replaced4)<- gsub("ptr.","Chimpanzee_", rownames(replaced4))
rownames(replaced4)<- gsub("ppa.","Bonobo_", rownames(replaced4))
rownames(replaced4)<- gsub("ggo.","Gorilla_", rownames(replaced4))
rownames(replaced4)<- gsub("ppy.","Orangutan_", rownames(replaced4))
rownames(replaced4)<- gsub("mml.","Macaque_", rownames(replaced4))
rownames(replaced4)<- gsub("mmu.","Mouse_", rownames(replaced4))
rownames(replaced4)<- gsub("mdo.","Opossum_", rownames(replaced4))
rownames(replaced4)<- gsub("oan.","Platypus_", rownames(replaced4))
rownames(replaced4)<- gsub("gga.","Chicken_", rownames(replaced4))
rownames(replaced4)<- gsub("br","Brain", rownames(replaced4))
rownames(replaced4)<- gsub("cb","Cerebellum", rownames(replaced4))
rownames(replaced4)<- gsub("ht","Heart", rownames(replaced4))
rownames(replaced4)<- gsub("kd","Kidney", rownames(replaced4))
rownames(replaced4)<- gsub("lv","Liver", rownames(replaced4))
rownames(replaced4)<- gsub("ts","Testis",rownames(replaced4))




####DATA MANIPULATION FOR PCA ANALYSIS
replaced5 <- log(replaced4+1) #log transformation
replaced6<- t(scale(t(replaced5))) #Scaling the data by the genes


data_pca<-prcomp(replaced6,center=TRUE) #PCA 
data_out <- as.data.frame(data_pca$x) #X table of PCA
data_out$species <- gsub("_.*","", row.names(replaced6)) #Create a column for species
data_out$tissue <- gsub(".*_","", row.names(replaced6)) #create a column for tissues
data_out$tissue <- sub("\\..*","", data_out$tissue)#create a column for tissues
data_out$species <- factor(data_out$species , levels = unique(data_out$species))#keep the order of species


##Calculation for percentage of variance explained by each component
eigs <- data_pca$sdev^2#Calculate percentage for PC values
percentage<- round(eigs/sum(eigs)*100,2)#Calculate percentage for PC values
percentage <- paste( colnames(data_out), "(", paste( as.character(percentage), "%", ")", sep="") ) #Calculate percentage for PC values

## PLOT FOR X
pdf("pca_amniote.pdf",height=5,width=6.5)
print(ggplot(data_out,aes(x=PC1,y=PC2,shape=species,color=tissue ))+
	geom_point()+
	scale_shape_manual(values=c(22,23,24,14,15,16,17,18,20,21))+ #Shapes
	scale_color_manual(values=c("#4575b4","#91bfdb","#8c510a","#af8dc3","#7fbf7b","#d73027"))+
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

##PLOT FOR LOADINGS
data_out_r <- as.data.frame(data_pca$rotation) #rotation data (loadings)
data_out_r$Symbol <- row.names(data_out_r) 
data_out_r<- join(data_out_r, ensembl, by="Symbol") #Include the Class information



pdf(file="pca_amniote.loadings.pdf",height=5,width=7)
print(ggplot(data_out_r,aes(x=PC1,y=PC2,label=Symbol,color=Class))+
	scale_color_manual(values = c("#117A65","#D2B4DE","#F1948A","#B03A2E","#85C1E9","#17202A","#7B7D7D"))+
	geom_point()+
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

##PLOT FOR LOADINGS SOME LABELS
data_out_r <- as.data.frame(data_pca$rotation) #rotation data (loadings)
data_out_r$Symbol <- row.names(data_out_r) 
data_out_r<- join(data_out_r, ensembl, by="Symbol") #Include the Class information


pdf(file="pca_amniote.loadings.somelabels.pdf",height=5,width=7)
print(ggplot(data_out_r,aes(x=PC1,y=PC2,label=Symbol,color=Class))+
	scale_color_manual(values = c("#117A65","#D2B4DE","#F1948A","#B03A2E","#85C1E9","#17202A","#7B7D7D"))+
	geom_point()+
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
	geom_text(data=subset(data_out_r, PC2 > 0.2 | PC2 < -0.1),hjust=0, vjust=-1)
	)
dev.off()



##BARPLOTS FOR EVERY GENE
new<-as.data.frame(replaced4)
new<-log(new+1)  #Log transformation
new$species <- gsub("_.*","", row.names(replaced6)) #Create a column for species
new$tissue <- gsub(".*_","", row.names(replaced6))#Create a column for tissue
new$tissue <- sub("\\..*","", new$tissue)#Create a column for tissue
new <- new %>% dplyr::select(species,tissue, everything())
new$species <- factor(new$species , levels = unique(new$species)) #keep the order
dat<-new


###For loop for the barplot 
for (i in c(3:dim(dat)[2])) {
	vname<- noquote(paste(colnames(new)[i]))
	#Plot the ggplot with error bars (both min and max)
	pdf(file=paste(vname, "barplot_amniote_logrpkm.pdf", sep=".") ,height=10,width=10)
	print(ggplot(new, aes_string(x=colnames(new)[2], y=colnames(new)[i], fill=colnames(new)[2])) + 
		geom_bar(stat="summary",fun.y="mean")+
		facet_wrap(~species,scales = "free")+
		scale_fill_manual(values = c("#4575b4","#91bfdb","#8c510a","#af8dc3","#7fbf7b","#d73027"))+
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
}




##TISSUE-VS-ALL ANALYSIS
#A FOR LOOP FOR THE ANALYSIS
for (sp in unique(new$species)) {
	new2<- new[with(new, new$species %in% sp),]
	new3 <- aggregate(. ~ tissue, new2[,-1], function(x) c(mean = mean(x))) #mean of same tissues
	new4<-t(new3)
	colnames(new4)<-new4[1,]
	new4<-new4[-1,]
	dat<-as.data.frame(new4)
	dat <- data.frame(sapply(dat, function(x) as.numeric(as.character(x))))
	row.names(dat)<-row.names(new4)

	#Calculate THRESHOLD for the specificity 
	res_vec <- vector() #res_vec file is a vector file that will contain all the residuals (tissue.vs.all) for all of the genes in all tissues
	for (i in 1:dim(dat)[2]) {
		res <- rlm(dat[,i] ~0+ as.vector(rowMeans(dat))) #rlm function is a robust linear model which is less affected by outliers. Starts from 0,0 point
		res_vec = c(res$residuals,res_vec) #this contains residuals for every gene in every tissuevsall combination
	}
	threshold <- 2.5*sd(res_vec) #The threshold is 2.5 times the standard deviation of all the residuals
	#Plotting highly expressed genes for each tissue
	specific_genes<-vector() #a vector that will store the specific genes information
	specific_tissue<-vector() #a vector that will store the tissue information
	for (i in 1:ncol(dat)) { 
		matrix1 <- matrix(seq(ncol(dat)), nrow = 1) #seq(X) should be same dimension as number of tissues
		matrix <- as.data.frame((matrix1))
		head<-colnames(dat)
		colnames(matrix)<- head
		pdf(file=paste(colnames(dat)[i],sp,"tissuespecific.pdf",sep="."),height=5,width=5)
		res <- rlm(dat[,i] ~0+ as.vector(rowMeans(dat))) #Forced to start from 0,0 point
		rstandard<-residuals(res)
		col = ifelse(rstandard>threshold, "red", "black")
		plot(dat[,i]~as.vector(rowMeans(dat)),xlab="mRNA mean abundance All Tissues", ylab=paste("mRNA abundance",colnames(dat)[i],sep=" "),pch=20, col=col, main= paste(colnames(dat)[i],sp,sep=" "))
		row_index <- c(as.numeric(names(rstandard[rstandard>threshold])))
		at_genes_with_specificity<-dat[row_index,]
			if (length (names(rstandard[rstandard>threshold]))> 0 ) {
				for (r in row_index) {
				vec = c ()
					if (dat[r,i] > rowMeans(dat[r,])) {
						print (length (dat[r,]))
		  				vec <- c(vec,r)
					}
					print (vec)
		 			matrix<-rbind(matrix, dat[unique(vec),])
				}
				matrix <- matrix[2:nrow(matrix),]
				row.names(matrix)
    			text(as.vector(rowMeans(matrix)),matrix[,i], labels=row.names(matrix), pos=2) 
   				specific_genes<-c(row.names(matrix),specific_genes) #feeding the vector with specific genes in each tissue
    			number<-as.numeric(NROW(c(matrix[,i])))
    			specific_tissue_input <-c(replicate(number, colnames(dat)[i]))
    			specific_tissue<-c(specific_tissue_input,specific_tissue)#feeding the vector with tissue information of each specific gene
    			abline(res)
			}
			else {
     		abline(res)
			}
	dev.off()
	}
}



