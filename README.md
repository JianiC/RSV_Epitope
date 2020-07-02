# RSV_Epitope

## Get Data

HRSV_0622.gb/ HRSV_0622.fasta file were downoladed form NCBI Genbank 'nucleotide database' with searching term "HRSV"
25790 records were get there

### metadata
1. Metadata was download from NCBI virus with Virus search "Human orthopneumovirus (HRSV), taxid 11250", saved in "HRSV_0622.csv"
2. genotype refernce were get with the note features from gb file using customized Gbmunge package
3. subtype information were collected with Blastn tool with cutomized HRSV_ref (n=2) database
`blastn -query /Users/jianichen1/Dropbox/RSV-genotyping_20/RSV_epitope/HRSV_0622.fasta -db HRSV_ref -outfmt "6 qseqid stitle sseqid" -max_target_seqs 1 >> HRSV_subtype.txt`
4. all metadata information were combined and save into "HRSV_0622_meta"
5. metadata were further code to remove space and were saved as "/Users/jianichen1/Dropbox/RSV/RSV_epitope/Data_collection/HRSV_0622_meta_code_0627.csv"

### Vaccine strain
1. artificial strain were saved under the dirctory "Artificial"
2. vaccine combinnant mutant strain: recombinant mutant rA2cp
mutant cp-RSV
mutant cpts-248
cp52, cold passaged deletion mutant
mutant cpts-248/404
recombinant D46/D53 strain

### RSV F gene data
1. RSV F gene were extracted from gb file with python script using 'CDS' annotation.
2. duplicate records and low quality sequence (< 80% coverage or > 20% gaps ) were removed.
3. Sequece id were code with corresponding metadata Accession|subtype|country|isolated date|genotype_ref
4. Artifical strains were removed with shell script
`sh rm_taxa.sh HRSV_F_taxa_0624.fasta HRSV_artificial_acession.txt > HRSV_F_0627clean.fasta`
5. The full HRSV_F dataset were seperated into A and B using subtype information with python scripts
6. Gene alignment were perfomed using MAFFT with normal option
7. Sequence with low quality were removed ( causing insertion; no start codon)

### RSV F protein data
1. RSV F nucleotide data were translated to protein in Genious using standard genetic code, and final stop codon has been removed.
2. File were export as fasta and further edit the extension as .pep were epitope analysis.
3. HTML file were converted into csv file with Excel (use group_merge.sh)
