# Diversity and Evolution of Computationally Predicted T cell Epitopes against Human Respiratory Syncytial Virus

## Data_collection

HRSV_0622.gb/ HRSV_0622.fasta file were downoladed form NCBI Genbank 'nucleotide database' with searching term "HRSV" 25790 records were get there

### metadata
1. Metadata was download from NCBI virus with Virus search "Human orthopneumovirus (HRSV), taxid 11250", saved in "HRSV_0622.csv"
2. genotype refernce were get with the note features from gb file using customized Gbmunge package, subtype information were collected with Blastn tool with cutomized HRSV_ref (n=2) database
`blastn -query /Users/jianichen1/Dropbox/RSV-genotyping_20/RSV_epitope/HRSV_0622.fasta -db HRSV_ref -outfmt "6 qseqid stitle sseqid" -max_target_seqs 1 >> HRSV_subtype.txt`

### Vaccine strain
1. artificial strain were saved under the dirctory "Artificial"
2. vaccine combinnant mutant strain: recombinant mutant rA2cp
mutant cp-RSV
mutant cpts-248
cp52, cold passaged deletion mutant
mutant cpts-248/404
recombinant D46/D53 strain

### Eptiope identification

Eitopes information collect from iVAX tool kits were cleaned up
1. Class II: Binding potential score top 5%, Binding alleles >3 for class II

Class I: Bind to any HLA allele top 1% binding prbablity

2. Conservation
3. tolerate epitope: cross-conserved to human protein
4. cross-conserved for both type A and type B

#### Epitope distribution
1. Epitope heatbar, compare with other pathogen
2. distribution map to 3D structure using pymol

### Epitope landscape
1. Epitope distance algorithm with cross-conseved epitope
T cell epitope immune distance (D) between two wild circulating strains (w_(1 )and w_(2 )) can be defined as the sum of Z-scaled binding probabilities of paired epitopes that are unable to induce cross-reactivity (non-cross conserved epitopes) within two protein sequences using equations (1.1 and 1.2).

d is the T cell immunity distance between a pair of epitopes, i and j are the non-cross conserved T cell epitopes from two protein sequences,

a is a class I or class II allele,

p is the predicted binding probability against allele a,

 A is a set of alleles. 


2. MDS + k-means cluster
3. validate with IEDB results

### Vaccine strain evaluation
calculate share epitope propotion

### Phylogenetic analysis
1. ML with Raxml
2. Ancestral sequence recosntruction "Treetime"
