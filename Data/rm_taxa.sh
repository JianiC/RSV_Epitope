#!/bin/bash

##https://bioinformatics.stackexchange.com/questions/3931/remove-delete-sequences-by-id-from-multifasta


awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' $1 | grep -vwf $2 | tr "\t" "\n" 


	