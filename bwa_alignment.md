#Running BWA

##BWA index of reference genome

```bash
bwa index Ssclerotiorum_v2_sorted.fasta

```


##Script for running bwa
###run_bwa.sh

```bash
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l virtual_free=1G

SEQ=$1
F_RUN=$2
F_RUN=$3
OutDir=$4
CurDir=$PWD

bwa mem $SEQ $F_RUN $R_RUN > $OutDir
```

##Calling bwa script

```bash
for StrainPath in $(ls -d /home/groups/harrisonlab/project_files/Sclerotinia_spp/qc_dna/paired/S.subartica/*); do 
echo $StrainPath
ProgDir=/home/ransoe/git_repos/tools/seq_tools/assemblers
TrimF1_Read=$(ls $StrainPath/F/*.fq.gz | grep 'run1'); 
TrimR1_Read=$(ls $StrainPath/R/*.fq.gz | grep 'run1');
TrimF2_Read=$(ls $StrainPath/F/*.fq.gz | grep 'run2'); 
TrimR2_Read=$(ls $StrainPath/R/*.fq.gz | grep 'run2'); 
TrimF3_Read=$(ls $StrainPath/F/*.fq.gz | grep 'run3'); 
TrimR3_Read=$(ls $StrainPath/R/*.fq.gz | grep 'run3'); 
FASTA=Ssclerotiorum_v2_sorted.fasta
echo $TrimF1_Read 
echo $TrimR1_Read
echo $TrimF2_Read
echo $TrimR2_Read
echo $TrimF3_Read
echo $TrimR3_Read
qsub $ProgDir/run_bwa.sh $FASTA $TrimF1_Read $TrimR1_Read Run1.sam
qsub $ProgDir/run_bwa.sh $FASTA $TrimF2_Read $TrimR2_Read Run2.sam
qsub $ProgDir/run_bwa.sh $FASTA $TrimF3_Read $TrimR3_Read Run3.sam
done
```

##Sam tools bwa sam —> bam 
###Sort and index for viewing in IGV

```bash
samtools view -bS Run1.sam > Run1.bam
samtools sort Run1.bam Run1_sorted
samtools index Run1_sorted.bam 

samtools view -bS Run2.sam > Run2.bam
samtools sort Run2.bam Run2_sorted
samtools index Run2_sorted.bam 

samtools view -bS Run3.sam > Run3.bam
samtools sort Run3.bam Run3_sorted
samtools index Run3_sorted.bam 


samtools mpileup -uf Ssclerotiorum_v2_sorted.fasta Run1_sorted.bam Run2_sorted.bam Run3_sorted.bam | bcftools view -cgb - > allruns.raw.bcf

bcftools view allruns.raw.bcf | /home/armita/all_idris/idris/prog/samtools-0.1.18/bcftools/vcfutils.pl vcf2fq -d10 -D500 > sorted.mpileup.filtered_d10_D500.fq

```
