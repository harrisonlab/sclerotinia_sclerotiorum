##Blasting to find regions of interest
###28S,18S,ITS,RPB2,HSP60,G3PDH and B-tubulin 
###Files stored in analysis/blast_homology


###Only interested in S.subartica but worth running through with all genomes for future reference

```bash
mkdir -p Sylvain_sequences
nano Sylvain_sequences/HSP60_blast.fa
nano Sylvain_sequences/RPB2_blast.fa
nano Sylvain_sequences/B_tubulin_blast.fa
nano Sylvain_sequences/G3PDH_blast.fa
nano Sylvain_sequences/Ref_Seq_18S.fa
nano Sylvain_sequences/Ref_Seq_28S.fa
nano Sylvain_sequences/Ref_Seq_ITS.fa
```

##Running against whole genome
###Run with Sylvain_blast script (Andrews blast_pipe.sh edited down to basics)

###HSP60
```bash
ProgDir=/home/ransoe/git_repos/tools/pathogen/blast
Query=analysis/blast_homology/Sylvain_sequences/HSP60_blast.fa 
for Assembly in $(ls assembly/spades/S.*/*/filtered_contigs/*_500bp_renamed.fasta); do
echo $Assembly
qsub $ProgDir/Sylvain_blast.sh $Query dna $Assembly
done
```

###RPB2
```bash
ProgDir=/home/ransoe/git_repos/tools/pathogen/blast
Query=analysis/blast_homology/Sylvain_sequences/RPB2_blast.fa 
for Assembly in $(ls assembly/spades/S.*/*/filtered_contigs/*_500bp_renamed.fasta); do
echo $Assembly
qsub $ProgDir/Sylvain_blast.sh $Query dna $Assembly
done
```

###Beta-tubulin
```bash
ProgDir=/home/ransoe/git_repos/tools/pathogen/blast
Query=analysis/blast_homology/Sylvain_sequences/B_tubulin_blast.fa 
for Assembly in $(ls assembly/spades/S.*/*/filtered_contigs/*_500bp_renamed.fasta); do
echo $Assembly
qsub $ProgDir/Sylvain_blast.sh $Query dna $Assembly
done
```

###G3PDH
```bash
ProgDir=/home/ransoe/git_repos/tools/pathogen/blast
Query=analysis/blast_homology/Sylvain_sequences/G3PDH_blast.fa 
for Assembly in $(ls assembly/spades/S.*/*/filtered_contigs/*_500bp_renamed.fasta); do
echo $Assembly
qsub $ProgDir/Sylvain_blast.sh $Query dna $Assembly
done
```

###Ref_Seq_18S
```bash
ProgDir=/home/ransoe/git_repos/tools/pathogen/blast
Query=analysis/blast_homology/Sylvain_sequences/Ref_Seq_18S.fa 
for Assembly in $(ls assembly/spades/S.*/*/filtered_contigs/*_500bp_renamed.fasta); do
echo $Assembly
qsub $ProgDir/Sylvain_blast.sh $Query dna $Assembly
done
```

###Ref_Seq_28S
```bash
ProgDir=/home/ransoe/git_repos/tools/pathogen/blast
Query=analysis/blast_homology/Sylvain_sequences/Ref_Seq_28S.fa 
for Assembly in $(ls assembly/spades/S.*/*/filtered_contigs/*_500bp_renamed.fasta); do
echo $Assembly
qsub $ProgDir/Sylvain_blast.sh $Query dna $Assembly
done
```

###Ref_Seq_ITS
```bash
ProgDir=/home/ransoe/git_repos/tools/pathogen/blast
Query=analysis/blast_homology/Sylvain_sequences/Ref_Seq_ITS.fa 
for Assembly in $(ls assembly/spades/S.*/*/filtered_contigs/*_500bp_renamed.fasta); do
echo $Assembly
qsub $ProgDir/Sylvain_blast.sh $Query dna $Assembly
done
```

##Turn into a gff file only the top 3 hits

```bash
ProgDir=/home/ransoe/git_repos/tools/pathogen/blast
for blast in B_tubulin_blast HSP60_blast RPB2_blast G3PDH_blast Ref_Seq_18S; do
echo $blast
for BlastHits in $(ls analysis/blast_homology/Sylvain_sequences_output/S.*/*/*_$blast.fa_hits.csv); do
Strain=$(echo $BlastHits | rev | cut -d '/' -f2 | rev)
Organism=$(echo $BlastHits | rev | cut -d '/' -f3 | rev)
echo $Strain
echo $Organism
HitsGff=analysis/blast_homology/Sylvain_sequences_output/$Organism/$Strain/"$Strain"_$blast.fa_hits.gff
Column2=Blast_homolog
NumHits=5
$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
done
done
```

```bash
ProgDir=/home/ransoe/git_repos/tools/pathogen/blast
for blast in Ref_Seq_28S Ref_Seq_ITS; do
echo $blast
for BlastHits in $(ls analysis/blast_homology/Sylvain_sequences_output/S.*/*/*_$blast.fa_hits.csv); do
Strain=$(echo $BlastHits | rev | cut -d '/' -f2 | rev)
Organism=$(echo $BlastHits | rev | cut -d '/' -f3 | rev)
echo $Strain
echo $Organism
HitsGff=analysis/blast_homology/Sylvain_sequences_output/$Organism/$Strain/"$Strain"_$blast.fa_hits.gff
Column2=Blast_homolog
NumHits=5
$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
done
done
```

##Running against predicted genes (CDS) using Sylvain_blast_genemodel.sh script

###HSP60
```bash
ProgDir=/home/ransoe/git_repos/tools/pathogen/blast
Query=analysis/blast_homology/Sylvain_sequences/HSP60_blast.fa 
for Assembly in $(ls gene_pred/augustus/*/*/*_aug_out.codingseq); do
echo $Assembly
qsub $ProgDir/Sylvain_blast_genemodel.sh $Query dna $Assembly
done
```

###RPB2
```bash
ProgDir=/home/ransoe/git_repos/tools/pathogen/blast
Query=analysis/blast_homology/Sylvain_sequences/RPB2_blast.fa 
for Assembly in $(ls gene_pred/augustus/*/*/*_aug_out.codingseq); do
echo $Assembly
qsub $ProgDir/Sylvain_blast_genemodel.sh $Query dna $Assembly
done
```


###Beta-tubulin
```bash
ProgDir=/home/ransoe/git_repos/tools/pathogen/blast
Query=analysis/blast_homology/Sylvain_sequences/B_tubulin_blast.fa 
for Assembly in $(ls gene_pred/augustus/*/*/*_aug_out.codingseq); do
echo $Assembly
qsub $ProgDir/Sylvain_blast_genemodel.sh $Query dna $Assembly
done
```

###G3PDH
```bash
ProgDir=/home/ransoe/git_repos/tools/pathogen/blast
Query=analysis/blast_homology/Sylvain_sequences/G3PDH_blast.fa 
for Assembly in $(ls gene_pred/augustus/*/*/*_aug_out.codingseq); do
echo $Assembly
qsub $ProgDir/Sylvain_blast_genemodel.sh $Query dna $Assembly
done
```

###Ref_Seq_18S
```bash
ProgDir=/home/ransoe/git_repos/tools/pathogen/blast
Query=analysis/blast_homology/Sylvain_sequences/Ref_Seq_18S.fa 
for Assembly in $(ls gene_pred/augustus/*/*/*_aug_out.codingseq); do
echo $Assembly
qsub $ProgDir/Sylvain_blast_genemodel.sh $Query dna $Assembly
done
```

###Ref_Seq_28S
```bash
ProgDir=/home/ransoe/git_repos/tools/pathogen/blast
Query=analysis/blast_homology/Sylvain_sequences/Ref_Seq_28S.fa 
for Assembly in $(ls gene_pred/augustus/*/*/*_aug_out.codingseq); do
echo $Assembly
qsub $ProgDir/Sylvain_blast_genemodel.sh $Query dna $Assembly
done
```

###Ref_Seq_ITS
```bash
ProgDir=/home/ransoe/git_repos/tools/pathogen/blast
Query=analysis/blast_homology/Sylvain_sequences/Ref_Seq_ITS.fa 
for Assembly in $(ls gene_pred/augustus/*/*/*_aug_out.codingseq); do
echo $Assembly
qsub $ProgDir/Sylvain_blast_genemodel.sh $Query dna $Assembly
done
```

```bash
ProgDir=/home/ransoe/git_repos/tools/pathogen/blast
for blast in B_tubulin_blast HSP60_blast RPB2_blast G3PDH_blast Ref_Seq_18S Ref_Seq_28S Ref_Seq_ITS; do
echo $blast
for BlastHits in $(ls analysis/blast_homology/Sylvain_sequences_output/S.*/*/*_EMR_singlestrand_aug_out.codingseq_*_$blast.fa_hits.csv); do
Strain=$(echo $BlastHits | rev | cut -d '/' -f2 | rev)
Organism=$(echo $BlastHits | rev | cut -d '/' -f3 | rev)
echo $Strain
echo $Organism
HitsGff=analysis/blast_homology/Sylvain_sequences_output/$Organism/$Strain/CodingSequence_"$Strain"_$blast.fa_hits.gff
Column2=Blast_homolog
NumHits=5
$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
done
done
```