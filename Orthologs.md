#Orthology analysis

```bash
My Sclerotinia species
S.minor S5
S.subartica HE1
S.sclerotiorum DG4
S.sclerotiorum P7
S.trifoliorum R316
```

#Previously assembled species
##S.sclerotiorum 1980

Make directory to perform analysis in.
#First comparison with the four proteomes I currently have and 1980 (published Sclerotinia sclerotiorum genome) 

  ProjDir=/home/groups/harrisonlab/project_files/Sclerotinia_spp
  cd $ProjDir
  IsolateAbrv=S5_HE1_Ssc1_R316_Ssc2
  WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  
  
  
#Format fasta files

NB: Taxon code must not be more than 4 characters/numbers long and must start with a character

## for S.minor S5
```bash
  Taxon_code=Smin
  Fasta_file=gene_pred/augustus/S.minor/S5/S5_EMR_singlestrand_aug_out.aa
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

## for S.subartica HE1
```bash
  Taxon_code=Ssub
  Fasta_file=gene_pred/augustus/S.subartica/HE1/HE1_EMR_singlestrand_aug_out.aa
  Id_field=1
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

## for S.sclerotiorum DG4

```bash
Taxon_code=Ssc1
Fasta_file=gene_pred/augustus/S.sclerotiorum/DG4/DG4_EMR_singlestrand_aug_out.aa
Id_field=1
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```
## for S.trifoliorum R316

```bash
Taxon_code=Stri
Fasta_file=gene_pred/augustus/S.trifoliorum/R316/R316_EMR_singlestrand_aug_out.aa
Id_field=1
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

## for S.sclerotiorum 1980

```bash
  Taxon_code=Ssc2
  Fasta_file=Sclerotinia_genome/Sclsc1_GeneCatalog_proteins_20110903.aa.fasta
  Id_field=4
  orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
  mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```  

##Filter proteins into good and poor sets.

  Input_dir=$WorkDir/formatted
  Min_length=10
  Max_percent_stops=20
  Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
  Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
  orthomclFilterFasta $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file


##Perform an all-vs-all blast of the proteins

  BlastDB=$WorkDir/blastall/$IsolateAbrv.db

  makeblastdb -in $Good_proteins_file -dbtype prot -out $BlastDB
  BlastOut=$WorkDir/all-vs-all_results.tsv
  mkdir -p $WorkDir/splitfiles

  SplitDir=/home/ransoe/git_repos/tools/seq_tools/feature_annotation/signal_peptides
  $SplitDir/splitfile_500.py --inp_fasta $Good_proteins_file --out_dir $WorkDir/splitfiles --out_base goodProteins

  ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/orthology  
  for File in $(find $WorkDir/splitfiles); do
    Jobs=$(qstat | grep 'blast_500' | grep 'qw' | wc -l)
    while [ $Jobs -gt 1 ]; do
      sleep 10
      printf "."
      Jobs=$(qstat | grep 'blast_500' | grep 'qw' | wc -l)
    done
    printf "\n"
    echo $File
    BlastOut=$(echo $File | sed 's/.fa/.tab/g')
    qsub $ProgDir/blast_500.sh $BlastDB $File $BlastOut
  done
  
###BELOW THIS POINT NOT COMPLETED!
##Merge the all-vs-all blast results

  MergeHits="$IsolateAbrv"_blast.tab
  printf "" > $MergeHits
  for Num in $(ls $WorkDir/splitfiles/*.tab | rev | cut -f1 -d '_' | rev | sort -n); do
    File=$(ls $WorkDir/splitfiles/*_$Num)
    cat $File
  done > $MergeHits

###Perform ortholog identification

  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/orthoMCL
  MergeHits="$IsolateAbrv"_blast.tab
  GoodProts=$WorkDir/goodProteins/goodProteins.fasta
  qsub $ProgDir/qsub_orthomcl.sh $MergeHits $GoodProts

###Plot venn diagrams:

  ProgDir=~/git_repos/emr_repos/tools/pathogen/orthology/venn_diagrams
  $ProgDir/ven_diag_5_way.R --inp $WorkDir/"$IsolateAbrv"_orthogroups.tab --out $WorkDir/"$IsolateAbrv"_orthogroups.pdf
Output was a pdf file of the venn diagram.

The following additional information was also provided. The format of the following lines is as follows:

Isolate name (total number of orthogroups) number of unique singleton genes number of unique groups of inparalogs

  [1] "Pcac (8494)"
  [1] 586
  [1] 126
  [1] "Pcap (7393)"
  [1] 348
  [1] 59
  [1] "Pinf (8079)"
  [1] 601
  [1] 107
  [1] "Ppar (8687)"
  [1] 732
  [1] 95
  [1] "Psoj (7592)"
  [1] 642
  [1] 153
  NULL


See downstream analysis in original P.cactorum file for specific analysis after this point
