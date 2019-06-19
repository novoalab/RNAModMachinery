# RNA Modification Machinery
Scripts used to extract curated lists of RNA modification enzymes and assess their tissue-specificity across multiple species and tissues, as well as cancer and normal tissues

## PART1: Search and extract sets of RNA modification ezymes in selected species across the tree of life

### FINDING RNA MODIFICATION ENZYMES in a specific SPECIES

``` 
find_homologs.sh <pfam_hmm> <fasta_reference_proteome> 

# Example: find_homologs.sh A_deamin.hmm.txt Saccharomyces_cerevisiae.fasta
``` 
- This script takes pfam profile as an input amd fasta reference proteome
- Output will be proteins that have similiar functional domains, e.g. A_deamin.hmm.txtSaccharomyces_cerevisiae.fasta
- Manual curation is performed to select candidates (based on literature, etc)

### ALIGNING THE ORTHOLOG PROTEINS AND BUILDING PHYLOGENETIC TREE FOR EACH GROUP OF ENZYMES

#### Extract fasta sequences for the proteins of interest in a list of orthologs proteins of a gene group
```
extract_fasta.sh <uniprot_ID_list>

# Example: extract_fasta.sh uniprotIDlist.txt
```
#### Align protein sequences with MAFFT align 
```
mafft.sh <FASTA>

# Example: mafft.sh uniprotIDlist.txt.named.fasta
```

#### Build build a phylogenetic tree with iqtree software using 
```
iqtree.sh <MAFFT_OUTPUT>

# Example: iqtree.sh uniprotIDlist.txt.named.fasta.mafft
```

## PART2: EXPRESSION ANALYSIS WITH HUMAN AND MOUSE DATASETS

### EXTRACTING RNA MODIFICATION MACHINERY (RMM) PROTEIN EXPRESSION FROM RNA SEQ DATASETS 
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
Rscript gtex_manipulation.R <ENCODE.Expression.File> <ENSEMBL_GeneSymbol_Class.File>

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
### Extracting RNA modification enzymes
Use Normalized RPKM Constitutive Exons tables for this
``` 
extract_modomics.RPKM.sh  <List of ENSEMBL IDs>  <RPKM table>

# Example: extract_modomics.RPKM.sh NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt.forR

``` 
### PCA for Amniotes
Use Normalized RPKM Constitutive Exons Table with row names are replaced with Gene Names instead of ENSEMBL ID
``` 
Rscript PCA_Amniote.R  <RPKM table>

# Example: Rscript PCA_Amniote.R  Amniote_RPKM_GeneNames


``` 
### PCA for Primates
Use Normalized RPKM Constitutive Exons Table with row names are replaced with Gene Names instead of ENSEMBL ID
``` 
Rscript PCA_Primate.R  <RPKM table>

# Example: Rscript PCA_Primate.R  Primate_RPKM_GeneNames

``` 


## PART4: EXPRESSION ANALYSIS IN SPERMATOGENESIS
### Within groups sum of squares (WSS) determining number of clusters 
Use log( Average of Normalized Expression in the cluster + 1) table for mouse RNA modification enzymes with Average expression values of clusters according to the paper (Green et al. 2018) {GC1 (Spermatogonia), GC2-3(Prelep-Spermatocyte), GC4-8(Spermatocytes), GC9-11(Spermatids), GC12 (Elongating Spermatids)}
cluster library is required
``` 
wss.R  <logExpressionTable>

# Example: Rscript wss.R scRNA_spermatogenesis.modomics.grouped
``` 

### K-means clustering the genes by their expression profile
Use the same input as wss.R and use the best cluster number
``` 
kmeans.R  <logExpressionTable> <Clusters>

# Example: Rscript kmeans.R scRNA_spermatogenesis.modomics.grouped 4
``` 

### Heatmap and PCA plots
Use the heatmap input which contains the Genes in a certain order based on clusters
``` 
heatmap_PCA_spermatogenesis.R <logExpressionTableHeatmapInput>

# Example: Rscript heatmap_PCA_spermatogenesis.R heatmapinput_logexpression
``` 
