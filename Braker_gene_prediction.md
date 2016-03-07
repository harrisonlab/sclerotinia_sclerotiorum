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

#New Braker script
```bash
33hrs
scp WTCHG_143994_18_*.fastq.gz ransoe@149.155.34.72:/home/groups/harrisonlab/project_files/Sclerotinia_spp/RNA_Seq/Timepoint_33/Rep_1
scp WTCHG_144185_16_*.fastq.gz ransoe@149.155.34.72:/home/groups/harrisonlab/project_files/Sclerotinia_spp/RNA_Seq/Timepoint_33/Rep_1

scp WTCHG_143996_19_*.fastq.gz ransoe@149.155.34.72:/home/groups/harrisonlab/project_files/Sclerotinia_spp/RNA_Seq/Timepoint_33/Rep_3
scp WTCHG_144187_19_*.fastq.gz ransoe@149.155.34.72:/home/groups/harrisonlab/project_files/Sclerotinia_spp/RNA_Seq/Timepoint_33/Rep_3

36hrs
scp WTCHG_143994_18_*.fastq.gz ransoe@149.155.34.72:/home/groups/harrisonlab/project_files/Sclerotinia_spp/RNA_Seq/Timepoint_36/Rep_1
scp WTCHG_144185_18_*.fastq.gz ransoe@149.155.34.72:/home/groups/harrisonlab/project_files/Sclerotinia_spp/RNA_Seq/Timepoint_36/Rep_1

scp WTCHG_143996_02_*.fastq.gz ransoe@149.155.34.72:/home/groups/harrisonlab/project_files/Sclerotinia_spp/RNA_Seq/Timepoint_36/Rep_3
scp WTCHG_144187_02_*.fastq.gz ransoe@149.155.34.72:/home/groups/harrisonlab/project_files/Sclerotinia_spp/RNA_Seq/Timepoint_36/Rep_3

```

```bash
cat *_1.fastq>33_Rep_1_merged_1.fastq
cat *_2.fastq>33_Rep_1_merged_2.fastq

cat *_1.fastq>33_Rep_3_merged_1.fastq
cat *_2.fastq>33_Rep_3_merged_2.fastq

cat *_1.fastq>36_Rep_1_merged_1.fastq
cat *_2.fastq>36_Rep_1_merged_2.fastq

cat *_1.fastq>36_Rep_3_merged_1.fastq
cat *_2.fastq>36_Rep_3_merged_2.fastq
```

##Quality check

```bash
for RNASeqdata in $(ls RNA_Seq/Timepoint_33/Rep_*/WTCHG_*.fastq); do
echo $RNASeqdata;
ProgDir=~/git_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RNASeqdata;
done
```
```bash
for RNASeqdata in $(ls RNA_Seq/Timepoint_33/Rep_*/*_merged_*.fastq); do
echo $RNASeqdata;
ProgDir=~/git_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RNASeqdata;
done
```

##Trim

```bash
qsub trimming.sh 33_Rep_1_merged_1.fastq 33_Rep_1_merged_2.fastq  33_Rep_1_merged_1_trim.fastq  33_Rep_1_merged_2_trim.fastq 
qsub trimming.sh 33_Rep_3_merged_1.fastq 33_Rep_3_merged_2.fastq  33_Rep_3_merged_1_trim.fastq  33_Rep_3_merged_2_trim.fastq 
qsub trimming.sh 36_Rep_1_merged_1.fastq 36_Rep_1_merged_2.fastq  36_Rep_1_merged_1_trim.fastq  36_Rep_1_merged_2_trim.fastq 
qsub trimming.sh 36_Rep_3_merged_1.fastq 36_Rep_3_merged_2.fastq  36_Rep_3_merged_1_trim.fastq  36_Rep_3_merged_2_trim.fastq 
```

#Quality check
```bash
for RNASeqdata in $(ls RNA_Seq/Timepoint_3*/Rep_*/*_merged_*_trim.fastq); do
echo $RNASeqdata;
ProgDir=~/git_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RNASeqdata;
done
```

##Align to published genome
```bash
Sclerotiniagenome=/Genomes/Sclerotinia/Ssclerotiorum_v2.fasta
ProgDir=/home/ransoe/git_repos/tools/seq_tools/RNAseq

for Filepath in $(ls -d RNA_Seq/Timepoint_3*/Rep_*); do
	echo $Filepath
	Rep=$(echo $Filepath| rev | cut -d '/' -f1 | rev)
	echo $Rep
	Timepoint=$(echo $Filepath | rev | cut -d '/' -f2 | rev)
	echo $Timepoint
	FileF=$(ls $Filepath/*_merged_1_trim.fastq)
	FileR=$(ls $Filepath/*_merged_2_trim.fastq)
	echo $FileF
	echo $FileR
	OutDir=alignment/trimmed/$Timepoint/$Rep
	echo $OutDir
	qsub $ProgDir/tophat_alignment_edit.sh $Sclerotiniagenome $FileF $FileR $OutDir
done
```

#Pull out aligned.bam file then sort, turn back into fastq and put reads back to F/R

```bash
samtools sort -n accepted_hits.bam accepted_hits_sorted
samtools view -H accepted_hits_sorted.bam
samtools bam2fq accepted_hits_sorted.bam>accepted_hits_sorted.fastq  
grep '@HISEQ.*/1' -A 3 accepted_hits_sorted.fastq --no-group-separator >accepted_hits_read_1.fastq
grep '@HISEQ.*/2' -A 3 accepted_hits_sorted.fastq --no-group-separator >accepted_hits_read_2.fastq
```

##Align against own assembled genomes
##Align against S.min

```bash
Sclerotiniagenome=assembly/spades/S.minor/S5/filtered_contigs/contigs_min_500bp_renamed.fasta 
ProgDir=/home/ransoe/git_repos/tools/seq_tools/RNAseq

for Filepath in $(ls -d alignment/trimmed/Timepoint_3*/Rep_*); do
echo $Filepath
Rep=$(echo $Filepath| rev | cut -d '/' -f2 | rev)
echo $Rep
Tech=$(echo $Filepath | rev | cut -d '/' -f1 | rev)
echo $Tech
FileF=$(ls $Filepath/accepted_hits_read_1.fastq)
FileR=$(ls $Filepath/accepted_hits_read_2.fastq)
echo $FileF
echo $FileR
OutDir=alignment/trimmed/S.minor/"$Rep"_"$Tech"
echo $OutDir
qsub $ProgDir/tophat_alignment_edit.sh $Sclerotiniagenome $FileF $FileR $OutDir
done
```

##Align against S.scl DG4

```bash
Sclerotiniagenome=assembly/spades/S.sclerotiorum/DG4/filtered_contigs/contigs_min_500bp_renamed.fasta 
ProgDir=/home/ransoe/git_repos/tools/seq_tools/RNAseq

for Filepath in $(ls -d alignment/trimmed/Timepoint_3*/Rep_*); do
echo $Filepath
Rep=$(echo $Filepath| rev | cut -d '/' -f2 | rev)
echo $Rep
Tech=$(echo $Filepath | rev | cut -d '/' -f1 | rev)
echo $Tech
FileF=$(ls $Filepath/accepted_hits_read_1.fastq)
FileR=$(ls $Filepath/accepted_hits_read_2.fastq)
echo $FileF
echo $FileR
OutDir=alignment/trimmed/S.sclerotiorum/DG4/"$Rep"_"$Tech"
echo $OutDir
qsub $ProgDir/tophat_alignment_edit.sh $Sclerotiniagenome $FileF $FileR $OutDir
done
```

##Align against S.tri

```bash
Sclerotiniagenome=assembly/spades/S.trifoliorum/R316/filtered_contigs/contigs_min_500bp_renamed.fasta 
ProgDir=/home/ransoe/git_repos/tools/seq_tools/RNAseq

for Filepath in $(ls -d alignment/trimmed/Timepoint_3*/Rep_*); do
echo $Filepath
Rep=$(echo $Filepath| rev | cut -d '/' -f2 | rev)
echo $Rep
Tech=$(echo $Filepath | rev | cut -d '/' -f1 | rev)
echo $Tech
FileF=$(ls $Filepath/accepted_hits_read_1.fastq)
FileR=$(ls $Filepath/accepted_hits_read_2.fastq)
echo $FileF
echo $FileR
OutDir=alignment/trimmed/S.trifoliorum/R316/"$Rep"_"$Tech"
echo $OutDir
qsub $ProgDir/tophat_alignment_edit.sh $Sclerotiniagenome $FileF $FileR $OutDir
done
```

##Align against S.sub

```bash
Sclerotiniagenome=assembly/spades/S.subartica/HE1/filtered_contigs/contigs_min_500bp_renamed.fasta 
ProgDir=/home/ransoe/git_repos/tools/seq_tools/RNAseq

for Filepath in $(ls -d alignment/trimmed/Timepoint_3*/Rep_*); do
echo $Filepath
Rep=$(echo $Filepath| rev | cut -d '/' -f2 | rev)
echo $Rep
Tech=$(echo $Filepath | rev | cut -d '/' -f1 | rev)
echo $Tech
FileF=$(ls $Filepath/accepted_hits_read_1.fastq)
FileR=$(ls $Filepath/accepted_hits_read_2.fastq)
echo $FileF
echo $FileR
OutDir=alignment/trimmed/S.subartica/HE1/"$Rep"_"$Tech"
echo $OutDir
qsub $ProgDir/tophat_alignment_edit.sh $Sclerotiniagenome $FileF $FileR $OutDir
done
```

#Align to S.scl P7

```bash
Sclerotiniagenome=assembly/spades/HiMem/S.sclerotiorum/P7/filtered_contigs/contigs_min_500bp_renamed.fasta 
ProgDir=/home/ransoe/git_repos/tools/seq_tools/RNAseq

for Filepath in $(ls -d alignment/trimmed/Timepoint_3*/Rep_*); do 
echo $Filepath 
Rep=$(echo $Filepath| rev | cut -d '/' -f2 | rev) 
echo $Rep 
Tech=$(echo $Filepath | rev | cut -d '/' -f1 | rev) 
echo $Tech 
FileF=$(ls $Filepath/accepted_hits_read_1.fastq) 
FileR=$(ls $Filepath/accepted_hits_read_2.fastq) 
echo $FileF 
echo $FileR 
OutDir=alignment/trimmed/S.sclerlotiorum/P7/"$Rep"_"$Tech" 
echo $OutDir 
qsub $ProgDir/tophat_alignment_edit.sh $Sclerotiniagenome $FileF $FileR $OutDir 
done
```
####NB. S.sclerotiorum spelt wrong

##Merge bam files

```bash
cp accepted_hits.bam ../accepted_hits_T33Rep1.bam
bamtools merge -in accepted_hits_T33Rep1.bam -in accepted_hits_T33Rep3.bam -in accepted_hits_T36Rep1.bam -in accepted_hits_T36Rep3.bam -out S.minor_accepted_hits.bam
```

#Braker

```bash
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp_renamed.fasta);do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo $Strain
echo $Organism
OutDir=gene_pred/braker/$Organism/$Strain
AcceptedHits=alignment/trimmed/$Organism/$Strain/"$Organism"_accepted_hits.bam
GeneModelName="$Organism"_"$Strain"_braker
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
```

####NB. S.subartica didn't work. S.minor just suddenly appeared overnight - check this.

##Re-run braker for S.subartica only
```bash
for Assembly in $(ls assembly/spades/S.subartica/*/filtered_contigs/contigs_min_500bp_renamed.fasta);do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo $Strain
echo $Organism
OutDir=gene_pred/braker/$Organism/$Strain
AcceptedHits=alignment/trimmed/$Organism/$Strain/"$Organism"_accepted_hits.bam
GeneModelName="$Organism"_"$Strain"_braker
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
```

#Braker for S.sclerotiorum P7 only

```bash
for Assembly in $(ls assembly/spades/HiMem/*/*/filtered_contigs/contigs_min_500bp_renamed.fasta);do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo $Strain
echo $Organism
OutDir=gene_pred/braker/$Organism/$Strain
AcceptedHits=alignment/trimmed/S.sclerlotiorum/$Strain/"$Organism"_P7_accepted_hits.bam
GeneModelName="$Organism"_"$Strain"_braker
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done

```