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
This was done with fastq-mcf on each pair using the script trimming.sh as below.

qsub trimming.sh fastqForwards fastqReverse OutputForwards OutputReverse

An example below.

```bash
qsub trimming.sh WTCHG_144185_02_1.fastq WTCHG_144185_02_2.fastq WTCHG_144185_02_1_trim.fq WTCHG_144185_02_2_trim.fq
```

###Quality check after trimmingq. 
```bash
ProgDir=~/git_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh RNA_Seq/Rep_1/tech_2/WTCHG_144185_02_1_trim.fq;
```

###Continued with un-trimmed data
NB: Test to see if trimmed data affects the alignment. If so this may need re-running.


##Align data
Align reads to published genome using tophat. As the reads are for both Sclerotinia and lettuce align to Sclerotinia 1980 genome.

```bash
Sclerotiniagenome=/Genomes/Sclerotinia/Ssclerotiorum_v2.fasta
ProgDir=/home/ransoe/git_repos/tools/seq_tools/RNAseq
	
for Filepath in $(ls -d RNA_Seq/Rep_1/tech_2); do
	echo $Filepath
	Rep=$(echo $Filepath| rev | cut -d '/' -f2 | rev)
	echo $Rep
	Tech=$(echo $Filepath | rev | cut -d '/' -f1 | rev)
	echo $Tech
	FileF=$(ls $Filepath/*_1.fq)
	FileR=$(ls $Filepath/*_2.fq)
	echo $FileF
	echo $FileR
	OutDir=alignment/$Rep/$Tech
	qsub $ProgDir/tophat_alignment_edit.sh $Sclerotiniagenome $FileF $FileR $OutDir
done
```

##Merge the bam files from technical reps together
e.g. Rep_1, tech_1 and tech_2 concatenated. 
```bash
tech1=Rep_1/tech_1/accepted_hits.bam
tech2=Rep_1/tech_2/accepted_hits.bam
bamtools merge -in accepted_hits_tech1.bam -in accepted_hits_tech2.bam -out Rep1_accepted_hits.bam
```


#Run BRAKER

cp /home/armita/.gm_key ~/.gm_key

ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
Assembly=/Genomes/Sclerotinia.Ssclerotiorum_v2.fasta
OutDir=gene_pred/braker/Sclerotinia_1980
AcceptedHits=alignment/Rep_1/Rep1_accepted_hits.bam
GeneModelName=Sclerotinia1980_braker
qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName


