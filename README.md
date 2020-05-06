# RNA Modification Machinery Analysis
Scripts used to extract curated lists of RNA modification enzymes and assess their tissue-specificity across multiple species and tissues, as well as cancer and normal tissues, used in the paper:

Begik O, Lucas MC, Liu H, Ramirez JM, Mattick JS and Novoa EM. Integrative analyses of the RNA modification machinery reveal tissue- and cancer-specific signatures. Genome Biology, May 2020. doi: https://doi.org/10.1101/830968


## DATA ACCOMPANYING THIS PAPER CAN BE FOUND [HERE](https://public-docs.crg.es/enovoa/public/website/Begik_RMP2020.html)

## PART1: Search and extract sets of RNA modification ezymes in selected species across the tree of life

### FINDING RNA MODIFICATION ENZYMES in a specific SPECIES
Required: HMMER
``` 
find_homologs.sh <pfam_hmm> <fasta_reference_proteome> 

# Example: find_homologs.sh A_deamin.hmm.txt Saccharomyces_cerevisiae.fasta
``` 
- This script takes pfam profile as an input and fasta reference proteome
- Output will be proteins that have similiar functional domains, e.g. A_deamin.hmm.txtSaccharomyces_cerevisiae.fasta
- Manual curation is performed to select candidates (based on literature, etc)

### ALIGNING THE ORTHOLOG PROTEINS AND BUILDING PHYLOGENETIC TREE FOR EACH GROUP OF ENZYMES

#### Extract fasta sequences for the proteins of interest in a list of orthologs proteins of a gene group
Required: perl
```
extract_fasta.sh <uniprot_ID_list>

# Example: extract_fasta.sh A_deamin
```
#### Align protein sequences with MAFFT align 
Required: mafft
```
mafft.sh <FASTA>

# Example: mafft.sh uniprotIDlist.txt.named.fasta
```

#### Build build a phylogenetic tree with iqtree software using 
Required : iqtree
```
iqtree.sh <MAFFT_OUTPUT>

# Example: iqtree.sh uniprotIDlist.txt.named.fasta.mafft
```

## PART2: EXPRESSION ANALYSIS WITH HUMAN AND MOUSE DATASETS

### EXTRACTING RNA MODIFICATION MACHINERY PROTEIN (RMP) EXPRESSION FROM RNA SEQ DATASETS 
Initially obtained list of Main RNA Writer Proteins and we added non-catalytic subunits, readers, erasers and other tRNA writer proteins from V de Crécy-Lagard et al - 2019. Therefore, we have 146 RMMs at the end. 

#### GTEx 
Extract TPM values for the RNA modification enzymes from the GTEx TPM dataset
```

Rscript gtex_manipulation.R <GTEX.Expression.File> <ENSEMBL_GeneSymbol_Class.File>

# Example: Rscript gtex_manipulation.R GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct human_id_symbol_class.tsv
```

#### ENCODE 
Extract TPM values for the RNA modification enzymes from the ENCODE TPM dataset
```
Rscript encode_manipulation.R <ENCODE.Expression.File> <ENSEMBL_GeneSymbol_Class.File>

# Example: Rscript encode_manipulation.R mm65.long.gene.with.expr.cshl.tsv mouse_id_symbol_class.tsv
```

### GTEx Tissue-Wide Expression Analysis
Scripts for tissue-specificity analysis and plots 

``` 
Rscript gtex_tissuewide.R  <TPM table>

# Example: Rscript gtex_tissuewide.R RMLP.GTEX.TissueAveraged.TPM.tsv
``` 

### ENCODE Tissue-Wide Expression Analysis
Scripts for tissue-specificity analysis and plots 

``` 
Rscript encode_tissuewide.R  <TPM table>

# Example: Rscript encode_tissuewide.R RMLP.encode.TPM.brainav.tsv
``` 


### Data similarity analysis between GTEx (Human) and ENCODE (Mouse)
Scripts for Pearson correlation analysis between two datasets
``` 
Rscript gtex_vs_encode_similarity.R  <GTEX data> <ENCODE data> 

# Example: Rscript gtex_vs_encode_similarity.R RMLP.GTEX.TissueAveraged.TPM.tsv RMLP.encode.TPM.brainav.tsv
``` 

## PART3: EXPRESSION ANALYSIS WITH AMNIOTE AND PRIMATE SPECIES
### Analysis for the Amniote Orthologs

``` 
Rscript kaessmann.amniote.R <input.expression.data> <ENSEMBL_GeneSymbol_Class.File>

# Example: Rscript kaessmann.amniote.R NormalizedRPKM_ConstitutiveExons_Amniote1to1Orthologues.txt human_id_symbol_class.tsv
``` 

### Analysis for the Primate Orthologs
``` 
Rscript kaessmann.primate.R <input.expression.data> <ENSEMBL_GeneSymbol_Class.File>

# Example: Rscript kaessmann.primate.R NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt human_id_symbol_class.tsv
``` 


## PART4: EXPRESSION ANALYSIS IN SPERMATOGENESIS
### Analysis of single-cell RNA sequencing data (Green et al., 2019)

``` 
Rscript spermatogenesis.R <spermatogenesis.expression.data> <ensembl_file>

# Example: Rscript spermatogenesis.R spermatogenesis_scRNA_averageexpression.tsv gene_hgnc_ensmus.tsv
``` 

## PART5: TUMOR AND NORMAL TISSUE ANALYSIS
### Extracting log(TPM+1) with reformatting

``` 
Rscript cancer_script1_datamanipulation.R <TCGA.GTEX.file> <ENSEMBL.file>

Example:  Rscript cancer_script1_datamanipulation.R RMLP.TcgaTargetGtex_rsem_gene_tpm_withheader.tsv human_id_symbol_class.tsv
``` 

### Boxplots of individual patients (log)
``` 
Rscript cancer_script2_boxplot.R <TCGA.GTEX.file>

Example:  Rscript cancer_script2_boxplot.R TCGA_GTEX_FINAL.log2.without3cancer.tsv
``` 

### Heatmap plot of average log(TPM) values

``` 
Rscript cancer_script3_mean_heatmap.R <TCGA.GTEX.file>

Example:  Rscript cancer_script3_mean_heatmap.R TCGA_GTEX_FINAL.log2.without3cancer.tsv
``` 

### Calculation of log2FC values for RMPs in tumor/normal pairs and plots

``` 
Rscript cancer_script4_log2FC.R <MedianLog Tumor and Normal TPM file>  <Original TPM File>
 
Example: Rscript cancer_script4_log2FC.R medianlog_tumor_normal.tpm.tsv RMLP.TcgaTargetGtex_rsem_gene_tpm_withheader.tsv
``` 

### Calculation of Dysregulation scores in tumor/normal pairs and plots
``` 
Rscript cancer_script5_Dysregulation.R <MedianLog Tumor and Normal TPM file for all genes> <MedianLog Tumor and Normal TPM file for RMPs>

Example:  Rscript cancer_script4_log2FC.R all_genes_logmedian_scatter_format.tsv medianlog_tumor_normal.tpm.tsv
``` 

### Plotting the expression values with stage information
``` 
Rscript cancer_script6_Stage.R <Gene Expression File> <ENSEMBL_GeneSymbol_Class.File> <clinical information> <Phenotype data>

Example:  Rscript cancer_script6_stageexpression.R RMLP.TcgaTargetGtex_rsem_gene_tpm_withheader.tsv human_id_symbol_class.tsv clinical.tsv TcgaTargetGTEX_phenotype.txt
``` 

### Plotting the survival curves
``` 
Rscript cancer_script7_SURVIVAL.R <Gene Expression File> <Survival data>

Example: Rscript cancer_script7_SURVIVAL.R TCGA_GTEX_FINAL.log2.without3cancer.tsv TCGA_survival_data
``` 

