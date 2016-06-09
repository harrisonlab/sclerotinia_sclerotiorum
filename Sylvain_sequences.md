##Blasting to find 28S,18S,ITS,RPB2,HSP60,G3PDH and B-tubulin in S.subartica
###

```bash
mkdir -p Sylvain_sequences
nano sclerotina_pub_genes/HSP60_blast.fa
```

```bash
ProgDir=/home/ransoe/git_repos/tools/pathogen/blast
Query=analysis/blast_homology/Sylvain_sequences/HSP60_blast.fa 
for Assembly in $(ls assembly/spades/S.*/*/filtered_contigs/*_500bp_renamed.fasta); do
echo $Assembly
qsub $ProgDir/Sylvain_blast.sh $Query dna $Assembly
done
```