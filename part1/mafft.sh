#MAFFT Alignment of multiple protein sequences
#Required : mafft
#mafft_alignment.sh
mafft --thread 10 --threadtb 5 --threadit 0 --reorder --maxiterate 0 --globalpair $1 > $1.mafft

#how to run
#mafft_alignment.sh file.fasta