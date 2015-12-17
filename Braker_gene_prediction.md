#Gene prediction using Braker

#RNA-Seq data

Gene prediction is trained using RNA-Seq data 

##Getting data onto the server
### From ncbi
Data used for this test analysis SRA code: SRR1915981 was downloaded using fastq-dump

```bash
mkdir -p raw_rna/external/S_sclerotiorum_test_RNA_data
fastq-dump -O raw_rna/external/S_sclerotiorum_test_RNA_data SRR1915981
```

### From my computer
Example shows fasta files moved from RNA-Seq_T42 Rep_3 on my computer to RNA_Seq/Rep_3 on east malling server.

(1) Change fasta file permissions
```bash
cd /Desktop/RNA-Seq_T42/Rep_3
chmod ugo+r WTCHG_*.fastq.gz
chmod ug+w  WTCHG_*.fastq.gz
```

(2) Copy to east malling server
```bash
scp WTCHG_*.fastq.gz ransoe@149.155.34.72:/home/groups/harrisonlab/project_files/Sclerotinia_spp/RNA_Seq/Rep_3
```

##Unzip files using gunzip
```bash
gunzip WTCHG_*.fastq.gz
```

##Concatenate the technical replicates
```bash
cat *_
```

##QC data
FileF=raw_rna/genbank/P.cactorum/10300/SRR1206032.fastq
  FileR=raw_rna/genbank/P.cactorum/10300/SRR1206033.fastq
  IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
  qsub /home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc/rna_qc_fastq-mcf.sh $FileF $FileR $IlluminaAdapters RNA


##Trim data

##Align data
Alignments of RNAseq reads were made against the 10300 Genome using tophat:

2.1) Alignment

ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
Genome=assembly/abyss/P.cactorum/10300/10300_abyss_53/10300_abyss-scaffolds_500bp_renamed.fa
FileF=qc_rna/raw_rna/genbank/P.cactorum/F/SRR1206032_trim.fq.gz
FileR=qc_rna/raw_rna/genbank/P.cactorum/R/SRR1206033_trim.fq.gz
OutDir=alignment/P.cactorum/10300
qsub $ProgDir/tophat_alignment.sh $Genome $FileF $FileR $OutDir

#Run BRAKER
screen -a set variables

  qlogin
  WorkDir=/tmp/braker
  ProjDir=/home/groups/harrisonlab/project_files/idris
  Assembly=$ProjDir/assembly/abyss/P.cactorum/10300/10300_abyss_53/10300_abyss-scaffolds_500bp_renamed.fa
  OutDir=$ProjDir/gene_pred/braker/P.cactorum/10300
move to working directory

  mkdir -p $WorkDir
  cd $WorkDir
  braker.pl \
    --cores 16 \
    --genome=$Assembly \
    --GENEMARK_PATH=/home/armita/prog/genemark/gm_et_linux_64/gmes_petap \
    --BAMTOOLS_PATH=/home/armita/prog/bamtools/bamtools/bin \
    --species=P.cactorum \
    --bam=$ProjDir/alignment/P.cactorum/10300/accepted_hits.bam
  mkdir -p $OutDir
  cp -r braker/* $OutDir/.

  rm -r $WorkDir