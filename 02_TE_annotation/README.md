# Pan_Bulbosum
Code for Article "A haplotype-resolved pangenome of the barley wild relative Hordeum bulbosum"



## EDTA annotation TE

The curated library TREP  (https://trep-db.uzh.ch/) was used as the curated library '--curatedlib', and the Manual selection of TE categories gives classification from EDTA (https://github.com/oushujun/EDTA/blob/master/util/TE_Sequence_Ontology.txt).

The high-confidence, evidence-based de-novo gene annotation of Morex V3 was used to remove genic sequences in the TE annotation '--cds'.

````shell
conda activate EDTA

EDTA.pl --genome $genome.fa -t 30 --anno 1 --cds Hv_Morex.pgsb.Jul2020.HC.cds_longest.fa --curatedlib Hordeum.200630.lib 

````