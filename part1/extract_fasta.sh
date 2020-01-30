#%%%%%%%%Extracting Fasta from Uniprot lists of families%%%%%%%
#extract_fasta.sh
#Extracting fasta sequences for the proteins of interest in a list of orthologs proteins of a gene group
#Dependencies: perl
#Notes: replace.pl (added) is a perl script to replace Uniprot ID with Gene and species name by using a file called Uniprot_Species_file 
#Uniprot_Species_file contains a table of 3 columns with Uniprot ID, GeneName and Species

#Create new directory where to work
mkdir dir.$1
cd dir.$1
cp ../$1 .

# List of IDs to uniprot fasta
cat $1 | while read line; do
     #wget http://www.uniprot.org/uniprot/$line.fasta
     curl -L -O http://www.uniprot.org/uniprot/$line.fasta > $line.fasta
done

#Merge all fasta sequences
cat *fasta > $1.fasta

#change the name of fasta to GeneName_Species (replace.pl and Uniprot_GeneName_Species files are in the folder)
perl ../replace.pl ../Uniprot_GeneName_Species $1.fasta > $1.named.fasta

#how to run
#extract_fasta.sh uniprotlist.file