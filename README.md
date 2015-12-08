# sclerotinia_sclerotiorum
Commands used fr the analysis of Sclerotinia spp. genomes

Sclerotinia sclerotiorum
====================


Useful notes

To kill multiple running command:

```bash
for num in $(seq 6359130 6359131); do
echo $num
qdel $num
done
```

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
  cp $RawDatDir/Run\ 4/DG4_S1_L001_R1_001.fastq.gz raw_data/paired/S.sclerotiorum/DG4/F/DG4_run4_F_fq.gz
  cp $RawDatDir/Run\ 4/DG4_S1_L001_R2_001.fastq.gz  raw_data/paired/S.sclerotiorum/DG4/R/DG4_run4_R_fq.gz
  # P7
  mkdir -p raw_data/paired/S.sclerotiorum/P7/F
  mkdir -p raw_data/paired/S.sclerotiorum/P7/R
  cp $RawDatDir/Run\ 1/P7_S2_L001_R1_001.fastq.gz raw_data/paired/S.sclerotiorum/P7/F/P7_run1_F_fq.gz
  cp $RawDatDir/Run\ 1/P7_S2_L001_R2_001.fastq.gz raw_data/paired/S.sclerotiorum/P7/R/P7_run1_R_fq.gz
  cp $RawDatDir/Run\ 2/P7_S2_L001_R1_001.fastq.gz raw_data/paired/S.sclerotiorum/P7/F/P7_run2_F_fq.gz
  cp $RawDatDir/Run\ 2/P7_S2_L001_R2_001.fastq.gz raw_data/paired/S.sclerotiorum/P7/R/P7_run2_R_fq.gz
  cp $RawDatDir/Run\ 3/P7_S2_L001_R1_001.fastq.gz raw_data/paired/S.sclerotiorum/P7/F/P7_run3_F_fq.gz
  cp $RawDatDir/Run\ 3/P7_S2_L001_R2_001.fastq.gz raw_data/paired/S.sclerotiorum/P7/R/P7_run3_R_fq.gz
  cp $RawDatDir/Run\ 4/P7_S2_L001_R1_001.fastq.gz raw_data/paired/S.sclerotiorum/P7/F/P7_run4_F_fq.gz
  cp $RawDatDir/Run\ 4/P7_S2_L001_R2_001.fastq.gz raw_data/paired/S.sclerotiorum/P7/R/P7_run4_R_fq.gz
  # HE1
  mkdir -p raw_data/paired/S.subartica/HE1/F
  mkdir -p raw_data/paired/S.subartica/HE1/R
  cp $RawDatDir/Run\ 1/HE1_S4_L001_R1_001.fastq.gz raw_data/paired/S.subartica/HE1/F/HE1_run1_F_fq.gz
  cp $RawDatDir/Run\ 1/HE1_S4_L001_R2_001.fastq.gz  raw_data/paired/S.subartica/HE1/R/HE1_run1_R_fq.gz
  cp $RawDatDir/Run\ 2/HE1_S4_L001_R1_001.fastq.gz raw_data/paired/S.subartica/HE1/F/HE1_run2_F_fq.gz
  cp $RawDatDir/Run\ 2/HE1_S4_L001_R2_001.fastq.gz  raw_data/paired/S.subartica/HE1/R/HE1_run2_R_fq.gz
  cp $RawDatDir/Run\ 3/HE1_S4_L001_R1_001.fastq.gz raw_data/paired/S.subartica/HE1/F/HE1_run3_F_fq.gz
  cp $RawDatDir/Run\ 3/HE1_S4_L001_R2_001.fastq.gz  raw_data/paired/S.subartica/HE1/R/HE1_run3_R_fq.gz
  cp $RawDatDir/Run\ 4/HE1_S4_L001_R1_001.fastq.gz raw_data/paired/S.subartica/HE1/F/HE1_run4_F_fq.gz
  cp $RawDatDir/Run\ 4/HE1_S4_L001_R2_001.fastq.gz  raw_data/paired/S.subartica/HE1/R/HE1_run4_R_fq.gz
  # R316
  mkdir -p raw_data/paired/S.trifoliorum/R316/F
  mkdir -p raw_data/paired/S.trifoliorum/R316/R
  cp $RawDatDir/Run\ 1/R316_S3_L001_R1_001.fastq.gz  raw_data/paired/S.trifoliorum/R316/F/R316_run1_F_fq.gz
  cp $RawDatDir/Run\ 1/R316_S3_L001_R2_001.fastq.gz  raw_data/paired/S.trifoliorum/R316/R/.
  cp $RawDatDir/Run\ 2/R316_S3_L001_R1_001.fastq.gz  raw_data/paired/S.trifoliorum/R316/F/.
  cp $RawDatDir/Run\ 2/R316_S3_L001_R2_001.fastq.gz  raw_data/paired/S.trifoliorum/R316/R/.
  cp $RawDatDir/Run\ 3/R316_S3_L001_R1_001.fastq.gz  raw_data/paired/S.trifoliorum/R316/F/.
  cp $RawDatDir/Run\ 3/R316_S3_L001_R2_001.fastq.gz  raw_data/paired/S.trifoliorum/R316/R/.
  cp $RawDatDir/Run\ 4/R316_S3_L001_R1_001.fastq.gz  raw_data/paired/S.trifoliorum/R316/F/.
  cp $RawDatDir/Run\ 4/R316_S3_L001_R2_001.fastq.gz  raw_data/paired/S.trifoliorum/R316/R/.

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
```bash
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
```

#Quast to summarise statistics

```bash
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
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/*_500bp_renamed.fasta); do
ProgDir=/home/ransoe/git_repos/tools/gene_prediction/cegma
    qsub $ProgDir/sub_cegma.sh $Assembly dna
    done
```
Results were viewed in completeness report, gff and fa. 

** Number of cegma genes present and complete:
** Number of cegma genes present and partial:

###Gene predictions

##Gene prediction part 1- Gene model predictions using Augustus

Gene prediction was performed for the neonectria genome.
CEGMA genes were used as Hints for the location of CDS.

```bash
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/*_500bp_renamed.fasta); do
 	Strain=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
 	Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo $Strain
    echo $Organism
    OutDir=gene_pred/augustus/$Organism/$Strain
    ProgDir=/home/ransoe/git_repos/tools/gene_prediction/augustus
    GeneModel=botrytis_cinerea
    qsub $ProgDir/submit_augustus.sh $GeneModel $Assembly false $OutDir
done
```
Count the number of genes in each gene prediction output
```bash
for Genes in $(ls gene_pred/augustus/S.*/*/*_EMR_singlestrand_aug_out.aa); do
cat $Genes |grep '>' |wc -l; 
done
```

** Number of genes predicted:

```bash
##Gene prediction part 2- Gene model predictions using ORFs

Open reading frame predictions were made using the run_ORF_finder.sh.

ProgDir=/home/ransoe/git_repos/tools/gene_prediction/ORF_finder
for Genome in $(ls assembly/spades/*/*/filtered_contigs/*_500bp_renamed.fasta); do
    qsub $ProgDir/run_ORF_finder.sh $Genome
done
``` 
 
 The Gff files from the the ORF finder are not in true Gff3 format. These were corrected using the following commands:
 
```bash
for ORF_Gff in $(ls gene_pred/ORF_finder/S.*/*/*_ORF.gff); do
	Strain=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
 	Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    ProgDir=home/ransoe/git_repos/tools/seq_tools/feature_annotation
    ORF_Gff_mod=gene_pred/ORF_finder/$Organism/$Strain/"$Strain"_ORF_corrected.gff3
    $ProgDir/gff_corrector.pl $ORF_Gff > $ORF_Gff_mod
done 
 ```
 
#Functional annotation

##Genes with homology to PHIbase
Predicted gene models were searched against the PHIbase database using tBLASTx.

```bash
for Subject in $(ls assembly/spades/*/*/filtered_contigs/*_500bp_renamed.fasta); do
    ProgDir=/home/ransoe/git_repos/tools/pathogen/blast
    Query=../../phibase/v3.8/PHI_accessions.fa
    qsub $ProgDir/blast_pipe.sh $Query protein $Subject
done
```

##Interproscan
Interproscan was used to give gene models functional annotations.

To run using my interproscan
```bash
for Genes in $(ls gene_pred/augustus/S.*/*/*_EMR_singlestrand_aug_out.aa); do
echo $Genes
ProgDir=/home/ransoe/git_repos/tools/seq_tools/feature_annotation/interproscan
$ProgDir/sub_interproscan.sh $Genes
done
```

To run using Andy's interproscan
```bash
for Genes in $(ls gene_pred/augustus/S.*/*/*_EMR_singlestrand_aug_out.aa); do
echo $Genes
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
$ProgDir/sub_interproscan.sh $Genes
done 
```


#Genomic analysis

#Signal peptide prediction


```bash
for Proteome in $(ls gene_pred/augustus/S.*/*/*_EMR_singlestrand_aug_out.aa); do
SplitfileDir=/home/ransoe/git_repos/tools/seq_tools/feature_annotation/signal_peptides
ProgDir=/home/ransoe/git_repos/tools/seq_tools/feature_annotation/signal_peptides
Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
SplitDir=gene_pred/augustus_split/$Organism/$Strain
mkdir -p $SplitDir
BaseName="$Organism""_$Strain"_augustus_preds
$SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
for File in $(ls $SplitDir/*_augustus_preds_*); do
Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
while [ $Jobs -ge 1 ]; do
sleep 10
printf "."
Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
done
printf "\n"
echo $File
# qsub $ProgDir/pred_sigP.sh $File
qsub $ProgDir/pred_sigP.sh $File signalp-4.1
done
done
```

The batch files of predicted secreted proteins needed to be combined into a
single file for each strain. This was done with the following commands:
```bash
for SplitDir in $(ls -d gene_pred/augustus_preds/augustus/S.*); do
Strain=$(echo $SplitDir | rev | cut -d '/' -f3 | rev)
Organism=$(echo $SplitDir | rev |cut -d '/' -f4 | rev)
InStringAA=''
InStringNeg=''
InStringTab=''
InStringTxt=''
SigpDir=augustus_signalP-4.1
for GRP in $(ls -l $SplitDir/*_augustus_preds_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do  
	InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_augustus_preds_$GRP""_sp.aa";  
	InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_augustus_preds_$GRP""_sp_neg.aa";  
	InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_augustus_preds_$GRP""_sp.tab";
	InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_augustus_preds_$GRP""_sp.txt";  
done
cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.aa
cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_neg_sp.aa
tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.tab
cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.txt
done 


```

The first analysis was based upon BLAST searches for genes known to be involved in toxin production



Top BLAST hits were used to annotate gene models.

```bash

```

** Blast results of note: **
  * 'Result A'
