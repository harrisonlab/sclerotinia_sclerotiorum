# Processing MinION and Illumina Genomes
# RNA-Seq data processing
### Data quality was visualised using fastqc:

```bash
for RawData in raw_rna/S.sclerotiorum/*/*.fastq.gz; do
	ProgDir=/home/ransoe/git_repos/tools/seq_tools/dna_qc
	echo $RawData;
	qsub $ProgDir/run_fastqc.sh $RawData
done
```
### Trimming was performed on data to trim adapters from sequences and remove poor quality data.This was done with fastq-mcf

```bash
for StrainPath in raw_rna/S.sclerotiorum/; do
		ProgDir=/home/ransoe/git_repos/tools/seq_tools/rna_qc
		IlluminaAdapters=/home/ransoe/git_repos/tools/seq_tools/ncbi_adapters.fa
		ReadsF=$(ls $StrainPath/F/*.fastq.gz)
		ReadsR=$(ls $StrainPath/R/*.fastq.gz)
		echo $ReadsF
		echo $ReadsR
		qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters RNA
	done
```

### Data was visualised again following trimming

```bash
	Data quality was visualised once again following trimming:

	for TrimData in qc_rna/raw_rna/*/*/*_trim.fq.gz; do
		ProgDir=/home/ransoe/git_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $TrimData
	done
```

# STAR alignment
## Edited Andy's star_sub.sh script to take previously made STAR index transferred from the York server

## New index stored in JointGenome/old_index/Index

## Copies of the fasta and gtf files used to create the jointgenome index are in JointGenome folder too.

```bash
ProgDir=/home/ransoe
FileF=qc_rna/raw_rna/S.sclerotiorum/F/all_forward_trim.fq.gz
FileR=qc_rna/raw_rna/S.sclerotiorum/R/all_reverse_trim.fq.gz
OutDir=alignment/star/S.sclerotiorum/
Index=JointGenome/old_index/index

qsub $ProgDir/star_sub_edit.sh $FileF $FileR $OutDir $Index
```
### This produced files in /home/groups/harrisonlab/project_files/Sclerotinia_spp

### These were copied to alignment

# Bam file processing
## View and index original file
```bash
samtools view star_aligmentAligned.sortedByCoord.out.bam | less
samtools index star_aligmentAligned.sortedByCoord.out.bam
```
## Extract out bam file of just Sclerotinia aligned reads
```bash
samtools view -b -o sclerotinia.bam star_aligmentAligned.sortedByCoord.out.bam CP017814.1 CP017815.1 CP017816.1 CP017817.1 CP017818.1 CP017819.1 CP017820.1 CP017821.1 CP017822.1 CP017823.1 CP017824.1 CP017825.1 CP017826.1 CP017827.1 CP017828.1 CP017829.1
```
## View output and index
```bash
samtools view sclerotinia.bam | less
```

## Stats of original file
```bash
samtools idxstats star_aligmentAligned.sortedByCoord.out.bam
```

### Output of stats of original file
```bash
sequence name, sequence length, # mapped reads and # unmapped reads

CP017814.1	3951982	934715	0
CP017815.1	3683506	966468	0
CP017816.1	3351453	1053583	0
CP017817.1	2873318	533706	0
CP017818.1	2822964	624443	0
CP017819.1	2483831	586180	0
CP017820.1	2434682	2838530	0
CP017821.1	2299506	756386	0
CP017822.1	2122865	561645	0
CP017823.1	2105496	412773	0
CP017824.1	2052242	469004	0
CP017825.1	1878461	436056	0
CP017826.1	1845946	393864	0
CP017827.1	1815632	414389	0
CP017828.1	1765292	323463	0
CP017829.1	1419421	364431	0
Dovetail_09Sept_Map_inspected_12-07-2015_1_v8_lg_1	214790941	10839001	0
Dovetail_09Sept_Map_inspected_12-07-2015_1_v8_lg_2	217164977	13316611	0
Dovetail_09Sept_Map_inspected_12-07-2015_1_v8_lg_3	257910700	11460886	0
Dovetail_09Sept_Map_inspected_12-07-2015_1_v8_lg_4	377529803	15619757	0
Dovetail_09Sept_Map_inspected_12-07-2015_1_v8_lg_5	339640866	26932002	0
Dovetail_09Sept_Map_inspected_12-07-2015_1_v8_lg_6	193110136	9674884	0
Dovetail_09Sept_Map_inspected_12-07-2015_1_v8_lg_7	195685357	12503888	0
Dovetail_09Sept_Map_inspected_12-07-2015_1_v8_lg_8	309698552	14091152	0
Dovetail_09Sept_Map_inspected_12-07-2015_1_v8_lg_9	204289203	11673322	0
```

## Index and view stats of new Sclerotinia only files
```bash
samtools index sclerotinia.bam
samtools idxstats sclerotinia.bam
```

```bash
CP017814.1	3951982	934715	0
CP017815.1	3683506	966468	0
CP017816.1	3351453	1053583	0
CP017817.1	2873318	533706	0
CP017818.1	2822964	624443	0
CP017819.1	2483831	586180	0
CP017820.1	2434682	2838530	0
CP017821.1	2299506	756386	0
CP017822.1	2122865	561645	0
CP017823.1	2105496	412773	0
CP017824.1	2052242	469004	0
CP017825.1	1878461	436056	0
CP017826.1	1845946	393864	0
CP017827.1	1815632	414389	0
CP017828.1	1765292	323463	0
CP017829.1	1419421	364431	0
Dovetail_09Sept_Map_inspected_12-07-2015_1_v8_lg_1	214790941	0	0
Dovetail_09Sept_Map_inspected_12-07-2015_1_v8_lg_2	217164977	0	0
Dovetail_09Sept_Map_inspected_12-07-2015_1_v8_lg_3	257910700	0	0
Dovetail_09Sept_Map_inspected_12-07-2015_1_v8_lg_4	377529803	0	0
Dovetail_09Sept_Map_inspected_12-07-2015_1_v8_lg_5	339640866	0	0
Dovetail_09Sept_Map_inspected_12-07-2015_1_v8_lg_6	193110136	0	0
Dovetail_09Sept_Map_inspected_12-07-2015_1_v8_lg_7	195685357	0	0
Dovetail_09Sept_Map_inspected_12-07-2015_1_v8_lg_8	309698552	0	0
Dovetail_09Sept_Map_inspected_12-07-2015_1_v8_lg_9	204289203	0	0
```
### Nothing aligned to lettuce genome in this bam file which is good. Not sure why the stats still show the lettuce chromosomes.

# Take bam file, back to fastq then re-run against my genomes?
## Script in alignment/star/S.sclerotiorum folder with $1 bam file $2 output1.fq $3 output2.fq
```bash
qsub bamtofastq.sh sclerotinia.bam sclerotinia_1.fq sclerotinia_2.fq
```
## Copied genome files and fastq files to each of the MinION folders in alignment

## Made indexes from within the folder with specific star_index.sh script in each
### Not the most elegant method but it works

## S.sclerotiorum P7
```bash
cp /home/groups/harrisonlab/project_files/Sclerotinia_spp/assembly/MinION/S.sclerotiorum/P7/S_scl_min_500bp_renamed_mtfree.fasta /home/groups/harrisonlab/project_files/Sclerotinia_spp/alignment/star/MinION_genomes/S.sclerotiorum

cp /home/groups/harrisonlab/project_files/Sclerotinia_spp/alignment/star/S.sclerotiorum/sclerotinia_*.fq /home/groups/harrisonlab/project_files/Sclerotinia_spp/alignment/star/MinION_genomes/S.sclerotiorum
```
### Make index
```bash
qsub star_index.sh S_scl_min_500bp_renamed_mtfree.fasta sclerotinia_1.fq sclerotinia_2.fq
```

## S.minor S5
```bash
cp /home/groups/harrisonlab/project_files/Sclerotinia_spp/assembly/MinION/S.minor/P7/S_minor_min_500bp_renamed.fasta /home/groups/harrisonlab/project_files/Sclerotinia_spp/alignment/star/MinION_genomes/S.sclerotiorum

cp /home/groups/harrisonlab/project_files/Sclerotinia_spp/alignment/star/S.sclerotiorum/sclerotinia_*.fq /home/groups/harrisonlab/project_files/Sclerotinia_spp/alignment/star/MinION_genomes/S.minor
```
### Make index
```bash
qsub star_index.sh S_minor_min_500bp_renamed.fasta sclerotinia_1.fq sclerotinia_2.fq
```

## S.subarctica HE1
```bash
cp /home/groups/harrisonlab/project_files/Sclerotinia_spp/assembly/MinION/S.subarctica/HE1/S_sub_min_500bp_renamed.fasta /home/groups/harrisonlab/project_files/Sclerotinia_spp/alignment/star/MinION_genomes/S.subarctica

cp /home/groups/harrisonlab/project_files/Sclerotinia_spp/alignment/star/S.sclerotiorum/sclerotinia_*.fq /home/groups/harrisonlab/project_files/Sclerotinia_spp/alignment/star/MinION_genomes/S.subarctica
```
### Make index
```bash
qsub star_index.sh S_sub_min_500bp_renamed.fasta sclerotinia_1.fq sclerotinia_2.fq
```

```bash
gzip sclerotinia_*.fq
```

# Running star using star_running.sh script in each folder
```bash
qsub star_running.sh sclerotinia_1.fq.gz sclerotinia_2.fq.gz aligned\ index
```

# Gene prediction

### Before braker predictiction was performed,double checked that I had the genemark key in my user area and copied it over from the genemark install
directory:

```bash
ls ~/.gm_key
cp /home/armita/prog/genemark/2018/gm_key_64 ~/.gm_key
```

# Gene training
# Braker
```bash
for Assembly in $(ls assembly/MinION/*/*/*_min_500bp_*.fasta);do
echo $Assembly
Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/braker/MinION_genomes/$Organism/"$Strain"_braker
AcceptedHits=$(ls alignment/star/MinION_genomes/$Organism/star_aligmentAligned.sortedByCoord.out.bam)
GeneModelName="$Organism"_"$Strain"_braker
rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
```

### Additional genes were added to Braker gene predictions, using CodingQuary in pathogen mode to predict additional regions.

### Firstly, aligned RNAseq data was assembled into transcripts using Cufflinks.

### Note - cufflinks doesn't always predict direction of a transcript and therefore features can not be restricted by strand when they are intersected.

```bash
for Assembly in $(ls assembly/MinION/*/*/*_min_500bp_*.fasta);do
Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/cufflinks/$Organism/$Strain/
mkdir -p $OutDir
AcceptedHits=$(ls alignment/star/MinION_genomes/$Organism/star_aligmentAligned.sortedByCoord.out.bam)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
done
```
# Coding quary

```bash
for Assembly in $(ls assembly/MinION/*/*/*_min_500bp_*.fasta);do
	Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
	Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
	echo "$Organism - $Strain"
	OutDir=gene_pred/codingquary/MinION_genomes/new/$Organism/$Strain
	CufflinksGTF=$(ls gene_pred/cufflinks/$Organism/$Strain/transcripts.gtf)
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
	qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
done
```
# Try again with reads aligned directly to my Genomes
## Moved bam files produced from subset of sclerotinia reads aligned to my genomes into "original_run" file in each of the star/MinION_genomes/* folders

```bash
FileF=/home/groups/harrisonlab/project_files/Sclerotinia_spp/qc_rna/raw_rna/S.sclerotiorum/F/all_forward_trim.fq.gz
FileR=/home/groups/harrisonlab/project_files/Sclerotinia_spp/qc_rna/raw_rna/S.sclerotiorum/R/all_reverse_trim.fq.gz
qsub star_running.sh $FileF $FileR aligned index
```
## Repeat masking
```bash
for Assembly in $(ls assembly/MinION/*/*/*_min_500bp_*.fasta); do
    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=repeat_masked/MinION_genomes/$Organism/"$Strain"/filtered_contigs
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
    qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
    qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
done
	```

###Running directly on blacklace10
	```bash
	for Assembly in $(ls assembly/MinION/*/*/*_min_500bp_*.fasta | tail -n2); do
	    Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
	    Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
	    echo "$Organism - $Strain"
	    OutDir=repeat_masked/MinION_genomes/$Organism/"$Strain"/filtered_contigs
	    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
	    $ProgDir/rep_modeling.sh $Assembly $OutDir
	    # qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
	done
		```

## Busco and QUAST
```bash
for Assembly in $(ls assembly/MinION/*/*/*_min_500bp_*.fasta); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

### Summary of quast/Busco
```bash
for File in $(ls assembly/MinION/*/*/run*/short_summary_*.txt); do
  Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Complete\t$Fragmented\t$Missing\t$Total"
  done
```

```bash
MinION	S.minor	1293	6	16	1315
MinION	S.sclerotiorum	1299	3	13	1315
MinION	S.subarctica	1300	1	14	1315
```
