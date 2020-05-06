#Maximum likelihood tree generation
#Required iqtree
#iqtree.sh
iqtree -s $1 -st AA -m TEST -bb 5000 -alrt 5000 -nt 8 -minsup 0.5
#5000 bootstrapping
#min supporting 50% 
#how to run
#iqtree.sh alignment.file
