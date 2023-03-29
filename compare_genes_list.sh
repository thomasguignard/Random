#!/bin/bash

echo "usage: compare_genes_list.sh genes_list  annotation_file"
echo 'genes not found will be in output file genes_not_found' 

GENE=$1
ANNOT=$2

wc -l $1
wc -l $2 

for i in `cat $GENE`; do if ! grep -qw $i $ANNOT; then echo $i " not found" >> genes_not_found.txt; fi ; done

exit 0
