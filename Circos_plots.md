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

## Align all illumina reads to each genome

###STOPPED HERE

```bash
for Reference in $(ls repeat_masked/MinION_genomes/*/*/filtered_contigs/*_contigs_unmasked.fa); do
for StrainPath in $(ls -d qc_dna/paired/*/* | grep -v -e 'S.trifoliorum' -e 'old_S.subartica' -e 'DG4'); do
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
echo $F_Read
echo $R_Read
Prefix="${Organism}_${Strain}"
OutDir=/genome_alignment/bwa/$Organism/$Strain/vs_${Reference}
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bwa/sub_bwa.sh $Prefix $Reference $F_Read $R_Read $OutDir
done
done
```

## ?
```bash
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
Fus2_genome=assembly/canu_spades_hybrid/F.oxysporum_fsp_cepae/Fus2/filtered_contigs/Fus2_contigs_renamed.fasta
$ProgDir/fasta2circos.py --genome $Fus2_genome --contig_prefix "" > tmp5/Fus2_genome.txt
for GffFile in $(ls analysis/blast_homology/*/Fus2_pacbio_test_merged/*_chr_*_gene_single_copy.aa_hits.gff); do
  echo $GffFile
  Chr=$(echo $GffFile | rev |cut -f1 -d'/' | rev | cut -f6 -d '_')
  $ProgDir/gff2circos_scatterplot.py --gff $GffFile --feature Chr"$Chr"_gene_homolog > tmp5/FoL_chr"$Chr"_genes.txt
done
  circos -conf /home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos/Fus2/Fus2_circos.conf -outputdir ./tmp5
```
