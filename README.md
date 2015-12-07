# sclerotinia_sclerotiorum
Commands used fr the analysis of Sclerotinia spp. genomes

Sclerotinia sclerotiorum
====================


Commands used during analysis of the Sclerotinia sclerotiorum genome. Note - all this work was performed in the directory:
```bash
mkdir -p /home/groups/harrisonlab/project_files/Sclerotinia_spp
cd /home/groups/harrisonlab/project_files/Sclerotinia_spp
```

The following is a summary of the work presented in this Readme:
Data organisation:
  * Preparing data
Draft Genome assembly
  * Data qc
  * Genome assembly
  * Repeatmasking
  * Gene prediction
  * Functional annotation
Genome analysis
  * Homology between predicted genes & published effectors


#Data organisation

Data was copied from the raw_data repository to a local directory for assembly
and annotation.


```bash
  RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/Sclerotinia
  # DG4
  mkdir -p raw_data/paired/S.sclerotiorum/DG4/F
  mkdir -p raw_data/paired/S.sclerotiorum/DG4/R
  cp $RawDatDir/Run\ 1/DG4_S1_L001_R1_001.fastq.gz raw_data/paired/S.sclerotiorum/DG4/F/DG4_run1_F_fq.gz
  cp $RawDatDir/Run\ 1/DG4_S1_L001_R2_001.fastq.gz  raw_data/paired/S.sclerotiorum/DG4/R/DG4_run1_R_fq.gz
  cp $RawDatDir/Run\ 2/DG4_S1_L001_R1_001.fastq.gz raw_data/paired/S.sclerotiorum/DG4/F/DG4_run2_F_fq.gz
  cp $RawDatDir/Run\ 2/DG4_S1_L001_R2_001.fastq.gz  raw_data/paired/S.sclerotiorum/DG4/R/DG4_run2_R_fq.gz
  cp $RawDatDir/Run\ 3/DG4_S1_L001_R1_001.fastq.gz raw_data/paired/S.sclerotiorum/DG4/F/DG4_run3_F_fq.gz
  cp $RawDatDir/Run\ 3/DG4_S1_L001_R2_001.fastq.gz  raw_data/paired/S.sclerotiorum/DG4/R/DG4_run3_R_fq.gz

  # P7
  mkdir -p raw_data/paired/S.sclerotiorum/P7/F
  mkdir -p raw_data/paired/S.sclerotiorum/P7/R
  cp $RawDatDir/Run\ 1/P7_S2_L001_R1_001.fastq.gz raw_data/paired/S.sclerotiorum/P7/F/P7_run1_F_fq.gz
  cp $RawDatDir/Run\ 1/P7_S2_L001_R2_001.fastq.gz raw_data/paired/S.sclerotiorum/P7/R/P7_run1_R_fq.gz
  cp $RawDatDir/Run\ 2/P7_S2_L001_R1_001.fastq.gz raw_data/paired/S.sclerotiorum/P7/F/P7_run2_F_fq.gz
  cp $RawDatDir/Run\ 2/P7_S2_L001_R2_001.fastq.gz raw_data/paired/S.sclerotiorum/P7/R/P7_run2_R_fq.gz
  cp $RawDatDir/Run\ 3/P7_S2_L001_R1_001.fastq.gz raw_data/paired/S.sclerotiorum/P7/F/P7_run3_F_fq.gz
  cp $RawDatDir/Run\ 3/P7_S2_L001_R2_001.fastq.gz raw_data/paired/S.sclerotiorum/P7/R/P7_run3_R_fq.gz

  # HE1
  mkdir -p raw_data/paired/S.subartica/HE1/F
  mkdir -p raw_data/paired/S.subartica/HE1/R
  cp $RawDatDir/Run\ 1/HE1_S4_L001_R1_001.fastq.gz raw_data/paired/S.subartica/HE1/F/HE1_run1_F_fq.gz
  cp $RawDatDir/Run\ 1/HE1_S4_L001_R2_001.fastq.gz  raw_data/paired/S.subartica/HE1/R/HE1_run1_R_fq.gz
  cp $RawDatDir/Run\ 2/HE1_S4_L001_R1_001.fastq.gz raw_data/paired/S.subartica/HE1/F/HE1_run2_F_fq.gz
  cp $RawDatDir/Run\ 2/HE1_S4_L001_R2_001.fastq.gz  raw_data/paired/S.subartica/HE1/R/HE1_run2_R_fq.gz
  cp $RawDatDir/Run\ 3/HE1_S4_L001_R1_001.fastq.gz raw_data/paired/S.subartica/HE1/F/HE1_run3_F_fq.gz
  cp $RawDatDir/Run\ 3/HE1_S4_L001_R2_001.fastq.gz  raw_data/paired/S.subartica/HE1/R/HE1_run3_R_fq.gz

  # R316
  mkdir -p raw_data/paired/S.trifoliorum/R316/F
  mkdir -p raw_data/paired/S.trifoliorum/R316/R
  cp $RawDatDir/Run\ 1/R316_S3_L001_R1_001.fastq.gz  raw_data/paired/S.trifoliorum/R316/F/R316_run1_F_fq.gz
  cp $RawDatDir/Run\ 1/R316_S3_L001_R2_001.fastq.gz  raw_data/paired/S.trifoliorum/R316/R/R316_run1_R_fq.gz
  cp $RawDatDir/Run\ 2/R316_S3_L001_R1_001.fastq.gz  raw_data/paired/S.trifoliorum/R316/F/R316_run2_F_fq.gz
  cp $RawDatDir/Run\ 2/R316_S3_L001_R2_001.fastq.gz  raw_data/paired/S.trifoliorum/R316/R/R316_run2_R_fq.gz
  cp $RawDatDir/Run\ 3/R316_S3_L001_R1_001.fastq.gz  raw_data/paired/S.trifoliorum/R316/F/R316_run3_F_fq.gz
  cp $RawDatDir/Run\ 3/R316_S3_L001_R2_001.fastq.gz  raw_data/paired/S.trifoliorum/R316/R/R316_run3_R_fq.gz


  # S5
  mkdir -p raw_data/paired/S.minor/S5/F
  mkdir -p raw_data/paired/S.minor/S5/R
  cp $RawDatDir/Run\ 1/Sminor_S5_L001_R1_001.fastq.gz  raw_data/paired/S.minor/S5/F/S5_run1_F_fq.gz
  cp $RawDatDir/Run\ 1/Sminor_S5_L001_R1_001.fastq.gz  raw_data/paired/S.minor/S5/R/S5_run1_R_fq.gz
  cp $RawDatDir/Run\ 2/Sminor_S5_L001_R1_001.fastq.gz  raw_data/paired/S.minor/S5/F/S5_run2_F_fq.gz
  cp $RawDatDir/Run\ 2/Sminor_S5_L001_R1_001.fastq.gz  raw_data/paired/S.minor/S5/R/S5_run2_R_fq.gz
  cp $RawDatDir/Run\ 3/Sminor_S5_L001_R1_001.fastq.gz  raw_data/paired/S.minor/S5/F/S5_run3_F_fq.gz
  cp $RawDatDir/Run\ 3/Sminor_S5_L001_R1_001.fastq.gz  raw_data/paired/S.minor/S5/R/S5_run3_R_fq.gz
```


#Data qc

programs: fastqc fastq-mcf kmc

Data quality was visualised using fastqc:

```bash
for RawData in $(ls raw_data/paired/S.*/*/*/*_fq.gz); do
echo $RawData;
ProgDir=~/git_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RawData;
done
```


Trimming was performed on data to trim adapters from sequences and remove poor quality data.
This was done with fastq-mcf.


```bash
for StrainPath in $(ls -d raw_data/paired/S.*/*); do
echo $StrainPath
IluminaAdapters=/home/ransoe/git_repos/tools/seq_tools/ncbi_adapters.fa
ProgDir=/home/ransoe/git_repos/tools/seq_tools/rna_qc
Read_F=$(ls $StrainPath/F/*_fq.gz | grep 'run1')
Read_R=$(ls $StrainPath/R/*_fq.gz | grep 'run1')
echo $Read_F
echo $Read_R
qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA

Read_F=$(ls $StrainPath/F/*_fq.gz | grep 'run2')
Read_R=$(ls $StrainPath/R/*_fq.gz | grep 'run2')
echo $Read_F
echo $Read_R
qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA

Read_F=$(ls $StrainPath/F/*_fq.gz | grep 'run3')
Read_R=$(ls $StrainPath/R/*_fq.gz | grep 'run3')
echo $Read_F
echo $Read_R
qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
done
```

Data quality was visualised once again following trimming:

```bash
for TrimData in $(ls qc_dna/paired/S.*/*/*/*.fq.gz); do
echo $TrimData;
ProgDir=~/git_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $TrimData;
done
```


kmer counting was performed using kmc.
This allowed estimation of sequencing depth and total genome size:

```bash
for StrainPath in $(ls -d qc_dna/paired/S.*/*); do 
echo $StrainPath
ProgDir=/home/ransoe/git_repos/tools/seq_tools/dna_qc
TrimF1_Read=$(ls $StrainPath/F/*.fq.gz | grep 'run1'); 
TrimR1_Read=$(ls $StrainPath/R/*.fq.gz | grep 'run1');
TrimF2_Read=$(ls $StrainPath/F/*.fq.gz | grep 'run2'); 
TrimR2_Read=$(ls $StrainPath/R/*.fq.gz | grep 'run2'); 
TrimF3_Read=$(ls $StrainPath/F/*.fq.gz | grep 'run3'); 
TrimR3_Read=$(ls $StrainPath/R/*.fq.gz | grep 'run3'); 
echo $TrimF1_Read 
echo $TrimR1_Read
echo $TrimF2_Read
echo $TrimR2_Read
echo $TrimF3_Read
echo $TrimR3_Read
qsub $ProgDir/kmc_kmer_counting.sh $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read
done
```

** Estimated Genome Size is:

** Esimated Coverage is:

#Assembly
Assembly was performed using: Spades

A range of hash lengths were used and the best assembly selected for subsequent analysis


```bash
  for StrainPath in $(ls -d qc_dna/paired/S.*/*); do
  echo $StrainPath
	ProgDir=/home/ransoe/git_repos/tools/seq_tools/assemblers/spades/multiple_libraries
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    echo $Strain
    echo $Organism
    TrimF1_Read=$(ls $StrainPath/F/*.fq.gz | grep 'run1');
    TrimR1_Read=$(ls $StrainPath/R/*.fq.gz | grep 'run1');
    TrimF2_Read=$(ls $StrainPath/F/*.fq.gz | grep 'run2');
    TrimR2_Read=$(ls $StrainPath/R/*.fq.gz | grep 'run2');
    TrimF3_Read=$(ls $StrainPath/F/*.fq.gz | grep 'run3');
    TrimR3_Read=$(ls $StrainPath/R/*.fq.gz | grep 'run3');
    echo $TrimF1_Read
    echo $TrimR1_Read
    echo $TrimF2_Read
    echo $TrimR2_Read
    echo $TrimF3_Read
    echo $TrimR3_Read
    OutDir=assembly/spades/$Organism/$Strain
    qsub $ProgDir/subSpades_3lib.sh $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir correct 10
  done
```

#Re-run assembly for Sclerotinia sclerotiorum P7
 for StrainPath in $(ls -d qc_dna/paired/S.sclerotiorum/P7); do
  echo $StrainPath
	ProgDir=/home/ransoe/git_repos/tools/seq_tools/assemblers/spades/multiple_libraries
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    echo $Strain
    echo $Organism
    TrimF1_Read=$(ls $StrainPath/F/*.fq.gz | grep 'run1');
    TrimR1_Read=$(ls $StrainPath/R/*.fq.gz | grep 'run1');
    TrimF2_Read=$(ls $StrainPath/F/*.fq.gz | grep 'run2');
    TrimR2_Read=$(ls $StrainPath/R/*.fq.gz | grep 'run2');
    TrimF3_Read=$(ls $StrainPath/F/*.fq.gz | grep 'run3');
    TrimR3_Read=$(ls $StrainPath/R/*.fq.gz | grep 'run3');
    echo $TrimF1_Read
    echo $TrimR1_Read
    echo $TrimF2_Read
    echo $TrimR2_Read
    echo $TrimF3_Read
    echo $TrimR3_Read
    OutDir=assembly/spades/$Organism/$Strain
    qsub $ProgDir/subSpades_3lib.sh $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $TrimF3_Read $TrimR3_Read $OutDir correct 10
  done 
  ```
  
Assemblies were summarised to allow the best assembly to be determined by eye.

** Assembly stats are:
  * Assembly size:
  * N50:
  * N80:
  * N20:
  * Longest contig:
  **


#Quast
#Re name the contigs to contig names to an acceptable format for NCBI

```bash
for OutDir in $(ls -d assembly/spades/S.*/*/filtered_contigs); do
echo $OutDir
    ProgDir=/home/ransoe/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
    AssFiltered=$OutDir/contigs_min_500bp.fasta
    AssRenamed=$OutDir/contigs_min_500bp_renamed.fasta
    echo $AssFiltered
    echo $AssRenamed
    printf '.\t.\t.\t.\n' > editfile.tab
    $ProgDir/remove_contaminants.py --inp $AssFiltered --out $AssRenamed --coord_file editfile.tab
    rm editfile.tab
done

#Quast to summarise statistics

ProgDir=/home/ransoe/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/*_500bp_renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo $Strain
    echo $Organism
    OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done

See report.txt output. 
```


# Repeat masking
Repeat masking was performed and used the following programs: Repeatmasker Repeatmodeler

The best assembly was used to perform repeatmasking

```bash
ProgDir=/home/ransoe/git_repos/tools/seq_tools/repeat_masking
for BestAss in $(ls assembly/spades/*/*/filtered_contigs/*_500bp_renamed.fasta); do
    echo $BestAss
    qsub $ProgDir/rep_modeling.sh $BestAss
    qsub $ProgDir/transposonPSI.sh $BestAss
done
 ```

** % bases masked by repeatmasker:

** % bases masked by transposon psi: **


# Gene Prediction
Gene prediction followed two steps:
Pre-gene prediction - Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
Gene models were used to predict genes in the Neonectria genome. This used results from CEGMA as hints for gene models.

## Pre-gene prediction
Quality of genome assemblies was assessed by looking for the gene space in the assemblies.

```bash

```

** Number of cegma genes present and complete:
** Number of cegma genes present and partial:

##Gene prediction

Gene prediction was performed for the neonectria genome.
CEGMA genes were used as Hints for the location of CDS.

```bash

```

** Number of genes predicted:

#Functional annotation

Interproscan was used to give gene models functional annotations.

```bash

```

```bash

```

#Genomic analysis
The first analysis was based upon BLAST searches for genes known to be involved in toxin production


##Genes with homology to PHIbase
Predicted gene models were searched against the PHIbase database using tBLASTx.

```bash

```

Top BLAST hits were used to annotate gene models.

```bash

```

** Blast results of note: **
  * 'Result A'
