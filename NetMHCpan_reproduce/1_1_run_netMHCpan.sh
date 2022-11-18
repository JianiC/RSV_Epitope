#!/bin/bash

## example to run: run_netMHCpan.sh -i "input_protein_file" -l "length_of_epitope" -A "list_of_Allel(sep by space ")
print_usage() {
  printf "Usage: ..."
}


#mkdir netMHCpan_out
## Shell script to run NetMHCpan


#####################################################
##### take command line args
#######################################################

while getopts 'i:l:A:' flag

do
	case "${flag}" in
        i) seq_pep=${OPTARG} ;;
        l) epitope_len=${OPTARG} ;; 
        A) A_list=${OPTARG} ;;            
        *) echo "Invalid option: -$flag" 
 #       	print_usage 
  #      	exit 1 ;;
    esac
done


#####################################################
##### default paramters to run netMHCpan
#######################################################
## for multiple allele use allele supetertype
A_supertype_full=()
A_supertype_short=(HLA-A*01:01 HLA-A*02:01 HLA-A*03:01 HLA-A*24:02 HLA-B*07:02 HLA-B*44:03)
#epitope_len="8,9,10,11" 


## no args use default
: ${epitope_len:="8,9,10,11"}

if [ -z "$A_list" ]
then
	A_list=("${A_supertype_short[@]}")
fi


#####################################################
##### run netMHCpan command line 
#######################################################

for allele in "${A_list[@]}"
do
	out=$(sed 's|[-*:]||g' <<< $allele)
	outxls="./netMHCpan_out/$out.xls"

	echo "netMHCpan -f $seq_pep" -l $epitope_len -a $allele -xls -xlsfile $outxls  	##test
#	netMHCpan -f $seq_pep -l $epitope_len -xls -a allele -xlsfile $outxls 
done 
