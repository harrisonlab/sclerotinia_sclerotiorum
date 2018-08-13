# Phylogeny

### FTP all other Sclerotinia genomes from ncbi for use within the phylogeny
```bash
Species='Sclerotinia_glacialis'
Strain='LMK_743'
OutDir=assembly/Scl_database/${Species}_${Strain}
mkdir -p $OutDir
wget -P $OutDir ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/162/575/GCA_002162575.1_ASM216257v1/GCA_002162575.1_ASM216257v1_genomic.fna.gz
gunzip $OutDir/*.gz
cp -s $PWD/$OutDir/*.fna $OutDir/genome.ctg.fa
```

```bash
Species='Sclerotinia_borealis'
Strain='F-4128'
OutDir=assembly/Scl_database/${Species}_${Strain}
mkdir -p $OutDir
wget -P $OutDir ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/503/235/GCA_000503235.1_SBOR_1/GCA_000503235.1_SBOR_1_genomic.fna.gz
gunzip $OutDir/*.gz
cp -s $PWD/$OutDir/*.fna $OutDir/genome.ctg.fa
```

```bash
Species='Sclerotinia_sclerotiorum'
Strain='1980'
OutDir=assembly/Scl_database/${Species}_${Strain}
mkdir -p $OutDir
wget -P $OutDir ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/945/GCF_000146945.2_ASM14694v2/GCF_000146945.2_ASM14694v2_cds_from_genomic.fna.gz
gunzip $OutDir/*.gz
cp -s $PWD/$OutDir/*.fna $OutDir/genome.ctg.fa
```

## Run Busco on all phylogenies
```bash
for Assembly in $(ls assembly/Scl_database/*/genome.ctg.fa); do
Organism=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
echo $Organism
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=gene_pred/busco/$Organism/
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

## Copy over all busco outputs from "assembly" of my genomes to gene_pred/busco

## Create a list of all BUSCO IDs
```bash
OutDir=analysis/MinION/popgen/busco_phylogeny
mkdir -p $OutDir
BuscoDb="ascomycota_odb9"
ls -1 /home/groups/harrisonlab/dbBusco/$BuscoDb/hmms/*hmm | rev | cut -f1 -d '/' | rev | sed -e 's/.hmm//' > $OutDir/all_buscos_"$BuscoDb".txt
```

## For each busco gene create a folder and move all single copy busco hits from each assembly to the folder. Then create a fasta file containing all the aligned reads for each busco gene for alignment later.

```bash
printf "" > analysis/MinION/popgen/busco_phylogeny/single_hits.txt
for Busco in $(cat analysis/MinION/popgen/busco_phylogeny/all_buscos_*.txt); do
echo $Busco
OutDir=analysis/MinION/popgen/busco_phylogeny/$Busco
mkdir -p $OutDir
for Fasta in $(ls gene_pred/busco/*/*/single_copy_busco_sequences/$Busco*.fna); do
Organism=$(echo $Fasta | rev | cut -f4 -d '/' | rev)
FileName=$(basename $Fasta)
cat $Fasta | sed "s/:.*.fasta:/:"$Organism":/g" | sed "s/:.*.fa:/:"$Organism":/g" > $OutDir/"$Organism"_"$Busco".fasta
done
cat $OutDir/*_*_"$Busco".fasta > $OutDir/"$Busco"_appended.fasta
SingleBuscoNum=$(cat $OutDir/"$Busco"_appended.fasta | grep '>' | cut -f2 -d ':' | sort | uniq | wc -l)
printf "$Busco\t$SingleBuscoNum\n" >> analysis/MinION/popgen/busco_phylogeny/single_hits.txt
done
```

## If all isolates have a single copy of a busco gene, move the appended fasta to a new folder
```bash
OutDir=analysis/MinION/popgen/busco_phylogeny/alignments
mkdir -p $OutDir
OrganismNum=$(cat analysis/MinION/popgen/busco_phylogeny/single_hits.txt | cut -f2 | sort -nr | head -n1)
for Busco in $(cat analysis/MinION/popgen/busco_phylogeny/all_buscos_*.txt); do
echo $Busco
HitNum=$(cat analysis/MinION/popgen/busco_phylogeny/single_hits.txt | grep "$Busco" | cut -f2)
if [ $HitNum == $OrganismNum ]; then
# cp analysis/popgen/busco_phylogeny/$Busco/"$Busco"_appended.fasta $OutDir/.
    cat analysis/MinION/popgen/busco_phylogeny/$Busco/"$Busco"_appended.fasta \
    | sed "s/$Busco://g" \
    | sed "s/genome.ctg.fa://g" \
    | sed "s/_contigs_unmasked.fa//g" \
    | sed -E "s/:.*//g" \
    | tr '.,:' '_' \
    > $OutDir/"$Busco"_appended.fasta
  fi
done
```

## Submit alignment for single copy busco genes with a hit in each organism
```bash
AlignDir=analysis/MinION/popgen/busco_phylogeny/alignments
CurDir=$PWD
cd $AlignDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/phylogenetics
qsub $ProgDir/sub_mafft_alignment.sh
cd $CurDir
```

## Trimming sequence alignments using Trim-Al
### Note - automated1 mode is optimised for ML tree reconstruction

```bash
OutDir=analysis/MinION/popgen/busco_phylogeny/trimmed_alignments
mkdir -p $OutDir
for Alignment in $(ls analysis/MinION/popgen/busco_phylogeny/alignments/*_appended_aligned.fasta); do
TrimmedName=$(basename $Alignment .fasta)"_trimmed.fasta"
echo $Alignment
trimal -in $Alignment -out $OutDir/$TrimmedName -automated1
done
```

```bash
for Alignment in $(ls analysis/MinION/popgen/busco_phylogeny/trimmed_alignments/*aligned_trimmed.fasta); do
Jobs=$(qstat | grep 'sub_RAxML' | grep 'qw' | wc -l)
while [ $Jobs -gt 2 ]; do
sleep 2s
printf "."
Jobs=$(qstat | grep 'sub_RAxML' | grep 'qw' | wc -l)
done		
printf "\n"
echo $Prefix
Prefix=$(basename $Alignment | cut -f1 -d '_')
OutDir=analysis/MinION/popgen/busco_phylogeny/RAxML/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/phylogenetics
qsub $ProgDir/sub_RAxML.sh $Alignment $Prefix $OutDir
done
```
