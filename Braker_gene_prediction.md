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
```bash
for RNASeqdata in $(ls RNA_Seq/Rep_*/tech_*/*.fastq); do
echo $RNASeqdata;
ProgDir=~/git_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RNASeqdata;
done

```
###Data doesn't appear to need trimming due to good quality and no adapters. Will however try trimming just to see.

##Data trimming

###This bit of code doesn't work

Trimming was performed on data to trim adapters from sequences and remove poor quality data.
This was done with fastq-mcf.

```bash
for RNASeqdata in $(ls -d RNA_Seq/*/*); do
echo $RNASeqdata;
ILLUMINA_ADAPTERS=/home/ransoe/git_repos/tools/seq_tools/ncbi_adapters.fa
ProgDir=/home/ransoe/git_repos/tools/seq_tools/rna_qc
FileF=$(ls $RNASeqdata/*_*.fastq | grep '1.fastq')
FileR=$(ls $RNASeqdata/*_*.fastq | grep '2.fastq')
echo $FileF
echo $FileR
F_FILE=$(echo $FileF | rev | cut -d "/" -f1 | rev | sed 's/.gz//')
R_FILE=$(echo $FileR| rev | cut -d "/" -f1 | rev | sed 's/.gz//')
echo $F_FILE
echo $R_FILE
F_OUT=$(echo "$F_FILE" | sed 's/.fq/_trim.fq/g' | sed 's/.fastq/_trim.fq/g')
R_OUT=$(echo "$R_FILE" | sed 's/.fq/_trim.fq/g' | sed 's/.fastq/_trim.fq/g')
echo $F_OUT
echo $R_OUT
qsub $ProgDir/rna_qc_fastq-mcf_RNA.sh $Read_F $Read_R $IluminaAdapters RNA
done
```

###Try running on each pair as so. 
```bash
fastq-mcf $ILLUMINA_ADAPTERS WTCHG_143994_02_1.fastq WTCHG_143994_02_2.fastq -o WTCHG_143994_02_1_trim.fq -o WTCHG_143994_02_2_trim.fq -C 1000000 -u -k 20 -t 0.01 -q 30 -p 5
qsub $ProgDir/run_fastqc.sh WTCHG_143994_02_1_trim.fq;
```

###Sort of works but still not needed so going to continue with un-trimmed data


##Align data
Align reads to published genome using tophat. As the reads are for both Sclerotinia and lettuce align to Sclerotinia 1980 genome.
#Test run with Rep_1, tech_1

```bash
Sclerotiniagenome=/home/groups/harrisonlab/project_files/Sclerotinia_spp/Genomes/Sclerotinia/Ssclerotiorum_v2.fasta
ProgDir=/home/ransoe/git_repos/tools/seq_tools/RNAseq
	
for Filepath in $(ls -d RNA_Seq/Rep_1/tech_1); do
	echo $Filepath
	Rep=$(echo $Filepath| rev | cut -d '/' -f2 | rev)
	echo $Rep
	Tech=$(echo $Filepath | rev | cut -d '/' -f1 | rev)
	echo $Tech
	FileF=$(ls $Filepath/*_1.fastq)
	FileR=$(ls $Filepath/*_2.fastq)
	echo $FileF
	echo $FileR
	OutDir=alignment/$Rep/$Tech
	qsub $ProgDir/tophat_alignment_edit.sh $Sclerotiniagenome $FileF $FileR $OutDir
done
```

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