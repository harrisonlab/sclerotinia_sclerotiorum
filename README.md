# sclerotinia_sclerotiorum
Commands used fr the analysis of Sclerotinia spp. genomes

Sclerotinia sclerotiorum 
====================

Commands used during analysis of the Sclerotinia sclerotiorum genome. Note - all this work was performed in the directory: /home/groups/harrisonlab/project_files/neonectria_ditissima

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
  cp $RawDatDir/Run\ 1/DG4_S1_L001_R1_001.fastq.gz raw_data/paired/S.sclerotiorum/DG4/F/.
  cp $RawDatDir/Run\ 1/DG4_S1_L001_R2_001.fastq.gz  raw_data/paired/S.sclerotiorum/DG4/R/.
  # P7
  mkdir -p raw_data/paired/S.sclerotiorum/P7/F
  mkdir -p raw_data/paired/S.sclerotiorum/P7/R
  cp $RawDatDir/Run\ 1/P7_S2_L001_R1_001.fastq.gz raw_data/paired/S.sclerotiorum/P7/F/.
  cp $RawDatDir/Run\ 1/P7_S2_L001_R2_001.fastq.gz raw_data/paired/S.sclerotiorum/P7/R/.
  # HE1
  mkdir -p raw_data/paired/S.subartica/HE1/F
  mkdir -p raw_data/paired/S.subartica/HE1/R
  cp $RawDatDir/Run\ 1/HE1_S4_L001_R1_001.fastq.gz raw_data/paired/S.subartica/HE1/F/.
  cp $RawDatDir/Run\ 1/HE1_S4_L001_R2_001.fastq.gz  raw_data/paired/S.subartica/HE1/R/.
  # R316
  mkdir -p raw_data/paired/S.trifoliorum/R316/F
  mkdir -p raw_data/paired/S.trifoliorum/R316/R
  cp $RawDatDir/Run\ 1/R316_S3_L001_R1_001.fastq.gz  raw_data/paired/S.trifoliorum/R316/F/.
  cp $RawDatDir/Run\ 1/R316_S3_L001_R2_001.fastq.gz  raw_data/paired/S.trifoliorum/R316/R/.
  # S5
  mkdir -p raw_data/paired/S.minor/S5/F
  mkdir -p raw_data/paired/S.minor/S5/R
  cp $RawDatDir/Run\ 1/Sminor_S5_L001_R1_001.fastq.gz  raw_data/paired/S.minor/S5/F/.
  cp $RawDatDir/Run\ 1/Sminor_S5_L001_R1_001.fastq.gz  raw_data/paired/S.minor/S5/R/.

```


#Data qc

programs: fastqc fastq-mcf kmc

Data quality was visualised using fastqc:

```bash

```

Trimming was performed on data to trim adapters from sequences and remove poor quality data.
This was done with fastq-mcf


```bash

```

Data quality was visualised once again following trimming:

```bash

```


kmer counting was performed using kmc.
This allowed estimation of sequencing depth and total genome size:

```bash

```

** Estimated Genome Size is:

** Esimated Coverage is:

#Assembly
Assembly was performed using: Spades

A range of hash lengths were used and the best assembly selected for subsequent analysis


```bash

```

Quast

```bash

```

Assemblies were summarised to allow the best assembly to be determined by eye.

** Assembly stats are:
  * Assembly size:
  * N50:
  * N80:
  * N20:
  * Longest contig:
  **




### Renaming contigs

Contigs were renamed in accordance with ncbi recomendations.

```bash

```

We run Quast again.
```bash

```

# Repeat masking
Repeat masking was performed and used the following programs: Repeatmasker Repeatmodeler

The best assembly was used to perform repeatmasking

```bash

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
