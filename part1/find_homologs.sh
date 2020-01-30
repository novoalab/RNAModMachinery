#%%%%%%%%FINDING ORTHOLOGS OF HUMAN RMWs IN GIVEN SPECIES%%%%%%
#--------------------------------------------------------------
#find_homologs.sh
#Dependencies : HMMER
#Find orthologs for RMWs in other species. This script takes pfam profile as an input and lists proteins that have the domain

#Dependencies : HMMER
f1=$(basename $1)
f2=$(basename $2)
hmmsearch $1 $2 > ${f1}.${f2}

#how to run the script
#find_homologs.sh pfam species.fasta

#For loop
#for in1 in pfam_list/*txt; do for in2 in proteomes/*fasta ;do ./find_homologs.sh $in1 $in2;done;done
