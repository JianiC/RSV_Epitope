# RSV_Epitope

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

#### Epitope distribution
1. Eptipe heatbar
2. distribution map to 3D structure using pymol

### Epitope landscape
1. Epitope distance algorithm
2. MDS
3. IEDB results evaluation

### Vaccine strain evalutation
