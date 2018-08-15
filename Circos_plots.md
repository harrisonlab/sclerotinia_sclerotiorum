#Circos plot

## Fastqc new S.subartica data

```bash
for StrainPath in $(ls -d raw_data/paired/S.subartica/HE1); do
echo $StrainPath
IluminaAdapters=/home/ransoe/git_repos/tools/seq_tools/ncbi_adapters.fa
ProgDir=/home/ransoe/git_repos/tools/seq_tools/rna_qc
Read_F=$(ls $StrainPath/new_F/*.fastq.gz)
Read_R=$(ls $StrainPath/new_R/*.fastq.gz)
echo $Read_F
echo $Read_R
qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
done
```

## Align all Illumina reads to each genome
### Concatenate multiple Illumina runs
```bash
for StrainPath in $(ls -d qc_dna/paired/*/* | grep -v -e 'S.trifoliorum' -e 'old_S.subartica' -e 'DG4'); do
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
echo $F_Read
echo $R_Read
cat $F_Read > $StrainPath/F/"$Strain"_appendedreads_F.fastq.gz
cat $R_Read > $StrainPath/R/"$Strain"_appendedreads_R.fastq.gz
done
```

```bash
for Reference in $(ls repeat_masked/MinION_genomes/*/*/filtered_contigs/*_contigs_unmasked.fa); do
  RefSpecies=$(echo $Reference | rev | cut -f4 -d '/' | rev)
  echo $RefSpecies
for StrainPath in $(ls -d qc_dna/paired/*/* | grep -v -e 'S.trifoliorum' -e 'old_S.subartica' -e 'DG4' -e 'S.subartica'); do
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_appendedreads_F.fastq.gz)
R_Read=$(ls $StrainPath/R/*_appendedreads_R.fastq.gz)
echo $F_Read
echo $R_Read
Prefix="${Organism}_${Strain}"
OutDir=/genome_alignment/bwa/$Organism/$Strain/vs_${RefSpecies}
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
$ProgDir/bwa/sub_bwa.sh $Prefix $Reference $F_Read $R_Read $OutDir
done
done
```

# Circos plots individual genomes

# Synteny between genomes

## Copy over reference sclerotiorum genome

### Fasta
```bash
Species='Sclerotinia_sclerotiorum_v2'
Strain='1980'
OutDir=assembly/Scl_database/${Species}_${Strain}
mkdir -p $OutDir
wget -P $OutDir ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/857/865/GCA_001857865.1_ASM185786v1/GCA_001857865.1_ASM185786v1_genomic.fna.gz
gunzip $OutDir/*.gz
cp -s $PWD/$OutDir/*.fna $OutDir/genome.ctg.fa
```
### gff
```bash
Species='Sclerotinia_sclerotiorum_v2'
Strain='1980'
OutDir=assembly/Scl_database/${Species}_${Strain}
wget -P $OutDir ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/857/865/GCA_001857865.1_ASM185786v1/GCA_001857865.1_ASM185786v1_genomic.gff.gz
gunzip $OutDir/*.gz
cp -s $PWD/$OutDir/*.gff $OutDir/genome.ctg.gff
```




