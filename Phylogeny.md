# Phylogeny

### FTP all other Sclerotinia genomes from ncbi for use within the phylogeny
```bash
Species='Sclerotinia_glacialis'
Strain='LMK_743'
OutDir=assembly/Scl_database/${Species}_${Strain}
mkdir -p $OutDir
wget -P $OutDir ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/162/575/GCA_002162575.1_ASM216257v1/GCA_002162575.1_ASM216257v1_genomic.fna.gz
gunzip $OutDir/*.gz
cp -s $PWD/$OutDir/*.fna $OutDir/genome.ctg.fa
```

```bash
Species='Sclerotinia_borealis'
Strain='F-4128'
OutDir=assembly/Scl_database/${Species}_${Strain}
mkdir -p $OutDir
wget -P $OutDir ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/503/235/GCA_000503235.1_SBOR_1/GCA_000503235.1_SBOR_1_genomic.fna.gz
gunzip $OutDir/*.gz
cp -s $PWD/$OutDir/*.fna $OutDir/genome.ctg.fa
```

```bash
Species='Sclerotinia_sclerotiorum'
Strain='1980'
OutDir=assembly/Scl_database/${Species}_${Strain}
mkdir -p $OutDir
wget -P $OutDir ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/945/GCF_000146945.2_ASM14694v2/GCF_000146945.2_ASM14694v2_cds_from_genomic.fna.gz
gunzip $OutDir/*.gz
cp -s $PWD/$OutDir/*.fna $OutDir/genome.ctg.fa
```

## Run Busco on all phylogenies
```bash
for Assembly in $(ls assembly/Scl_database/*/genome.ctg.fa); do
Organism=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
echo $Organism
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=gene_pred/busco/$Organism/
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

## NB: realised ran on old Sclerotinia ref genome. 
Re-run busco on new genome

## Run Busco on all phylogenies
```bash
for Assembly in $(ls assembly/Scl_database/Sclerotinia_sclerotiorum_v2_1980/genome.ctg.fa); do
Organism=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
echo $Organism
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=gene_pred/busco/$Organism/
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

## Copy over all busco outputs from "assembly" of my genomes to gene_pred/busco

## Create a list of all BUSCO IDs
```bash
OutDir=analysis/MinION/popgen/busco_phylogeny
mkdir -p $OutDir
BuscoDb="ascomycota_odb9"
ls -1 /home/groups/harrisonlab/dbBusco/$BuscoDb/hmms/*hmm | rev | cut -f1 -d '/' | rev | sed -e 's/.hmm//' > $OutDir/all_buscos_"$BuscoDb".txt
```

## For each busco gene create a folder and move all single copy busco hits from each assembly to the folder. Then create a fasta file containing all the aligned reads for each busco gene for alignment later.

```bash
printf "" > analysis/MinION/popgen/busco_phylogeny/single_hits.txt
for Busco in $(cat analysis/MinION/popgen/busco_phylogeny/all_buscos_*.txt); do
echo $Busco
OutDir=analysis/MinION/popgen/busco_phylogeny/$Busco
mkdir -p $OutDir
for Fasta in $(ls gene_pred/busco/*/*/single_copy_busco_sequences/$Busco*.fna); do
Organism=$(echo $Fasta | rev | cut -f4 -d '/' | rev)
FileName=$(basename $Fasta)
cat $Fasta | sed "s/:.*.fasta:/:"$Organism":/g" | sed "s/:.*.fa:/:"$Organism":/g" > $OutDir/"$Organism"_"$Busco".fasta
done
cat $OutDir/*_*_"$Busco".fasta > $OutDir/"$Busco"_appended.fasta
SingleBuscoNum=$(cat $OutDir/"$Busco"_appended.fasta | grep '>' | cut -f2 -d ':' | sort | uniq | wc -l)
printf "$Busco\t$SingleBuscoNum\n" >> analysis/MinION/popgen/busco_phylogeny/single_hits.txt
done
```

## If all isolates have a single copy of a busco gene, move the appended fasta to a new folder
```bash
OutDir=analysis/MinION/popgen/busco_phylogeny/alignments
mkdir -p $OutDir
OrganismNum=$(cat analysis/MinION/popgen/busco_phylogeny/single_hits.txt | cut -f2 | sort -nr | head -n1)
for Busco in $(cat analysis/MinION/popgen/busco_phylogeny/all_buscos_*.txt); do
echo $Busco
HitNum=$(cat analysis/MinION/popgen/busco_phylogeny/single_hits.txt | grep "$Busco" | cut -f2)
if [ $HitNum == $OrganismNum ]; then
# cp analysis/popgen/busco_phylogeny/$Busco/"$Busco"_appended.fasta $OutDir/.
    cat analysis/MinION/popgen/busco_phylogeny/$Busco/"$Busco"_appended.fasta \
    | sed "s/$Busco://g" \
    | sed "s/genome.ctg.fa://g" \
    | sed "s/_contigs_unmasked.fa//g" \
    | sed -E "s/:.*//g" \
    | tr '.,:' '_' \
    > $OutDir/"$Busco"_appended.fasta
  fi
done
```

## Submit alignment for single copy busco genes with a hit in each organism
```bash
AlignDir=analysis/MinION/popgen/busco_phylogeny/alignments
CurDir=$PWD
cd $AlignDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/phylogenetics
qsub $ProgDir/sub_mafft_alignment.sh
cd $CurDir
```

## Trimming sequence alignments using Trim-Al
### Note - automated1 mode is optimised for ML tree reconstruction

```bash
OutDir=analysis/MinION/popgen/busco_phylogeny/trimmed_alignments
mkdir -p $OutDir
for Alignment in $(ls analysis/MinION/popgen/busco_phylogeny/alignments/*_appended_aligned.fasta); do
TrimmedName=$(basename $Alignment .fasta)"_trimmed.fasta"
echo $Alignment
trimal -in $Alignment -out $OutDir/$TrimmedName -automated1
done
```

```bash
for Alignment in $(ls analysis/MinION/popgen/busco_phylogeny/trimmed_alignments/*aligned_trimmed.fasta); do
Jobs=$(qstat | grep 'sub_RAxML' | grep 'qw' | wc -l)
while [ $Jobs -gt 2 ]; do
sleep 2s
printf "."
Jobs=$(qstat | grep 'sub_RAxML' | grep 'qw' | wc -l)
done		
printf "\n"
echo $Prefix
Prefix=$(basename $Alignment | cut -f1 -d '_')
OutDir=analysis/MinION/popgen/busco_phylogeny/RAxML/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/phylogenetics
qsub $ProgDir/sub_RAxML.sh $Alignment $Prefix $OutDir
done
```

### Run Astral to build a consensus phylogeny from a collective set of "best phylogenies" from each BUSCO locus

## Note - "Recent versions of ASTRAL output a branch support value even without bootstrapping. Our analyses have revealed that this form of support is more reliable than bootstrapping (under the conditions we explored). Nevertheless, you may want to run bootstrapping as well."
## Tutorial tips: https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md#running-with-unresolved-gene-trees

```bash
OutDir=analysis/MinION/popgen/busco_phylogeny/ASTRAL
mkdir -p $OutDir
# cat analysis/popgen/busco_phylogeny/RAxML/*/RAxML_bestTree.* > $OutDir/Pcac_phylogeny.appended.tre
# Andy-Taxa names were noted to be incorrect at this point and were corrected
# Me- no harm in keeping the Sed command in
cat analysis/MinION/popgen/busco_phylogeny/RAxML/*/RAxML_bestTree.*  | sed -r "s/CTG.\w+:/:/g" | sed 's/__/_/g' > $OutDir/Scl_phylogeny.appended.tre

# InTree=$(ls /home/armita/prog/Astral/Astral/test_data/song_primates.424.gene.tre)
# -
# Trimm back branches that have less than 10% bootstrap support for each tree
# in the given file
# -
/home/armita/prog/newick_utilities/newick_utils/src/nw_ed $OutDir/Scl_phylogeny.appended.tre 'i & b<=10' o > $OutDir/Scl_phylogeny.appended.trimmed.tre
# -
# Calculate combined tree
# -
ProgDir=/home/armita/prog/Astral/Astral
# java -Xmx1000M -jar $ProgDir/astral.5.6.1.jar -i $OutDir/Pcac_phylogeny.appended.trimmed.tre -o $OutDir/Pcac_phylogeny.consensus.tre | tee 2> $OutDir/Pcac_phylogeny.consensus.log
java -Xmx1000M -jar $ProgDir/astral.5.6.1.jar -i $OutDir/Scl_phylogeny.appended.tre -o $OutDir/Scl_phylogeny.consensus.tre | tee 2> $OutDir/Scl_phylogeny.consensus.log
java -Xmx1000M -jar $ProgDir/astral.5.6.1.jar -q $OutDir/Scl_phylogeny.consensus.tre -i $OutDir/Scl_phylogeny.appended.tre -o $OutDir/Scl_phylogeny.consensus.scored.tre 2> $OutDir/Scl_phylogeny.consensus.scored.log
```

```bash
cat Scl_phylogeny.consensus.scored.geneious.tre | sed 's/:2/:1/g' > Scl_phylogeny.consensus.scored.geneious2.tre
```

## Finish plotting in R.

```bash

setwd("/Volumes/Backup/Sclerotinia_genomes/Phylogeny")

library(ape)
library(ggplot2)

library(ggtree)
library(phangorn)
library(treeio)

tree <- read.tree("/Volumes/Backup/Sclerotinia_genomes/Phylogeny/Scl_phylogeny.consensus.scored.geneious2.tre")

mydata <- read.csv("/Volumes/Backup/Sclerotinia_genomes/Phylogeny/traits.csv", stringsAsFactors=FALSE)
rownames(mydata) <- mydata$Isolate
mydata <- mydata[match(tree$tip.label,rownames(mydata)),]

t <- ggtree(tree, aes(linetype=nodes$support)) # Core tree
# Adjust terminal branch lengths:
branches <- t$data
tree$edge.length[branches$isTip] <- 1.0
#Tree <- branches$branch.length
#rescale_tree(t, branches$branch.length)

t <- t + geom_treescale(offset=-1.0, fontsize = 3) # Add scalebar
# t <- t + xlim(0, 0.025) # Add more space for labels

# Colouring labels by values in another df
t <- t %<+% mydata # Allow colouring of nodes by another df
#t <- t + geom_tiplab(aes(color=Source), size=3, hjust=0) +
scale_color_manual(values=c("gray39","black")) # colours as defined by col2rgb

tips <- data.frame(t$data)
tips$label <- tips$ID
t <- t + geom_tiplab(data=tips, aes(color=Source), size=3, hjust=0, offset = +0.1) +
  scale_color_manual(values=c("gray39","black")) # colours as defined by col2rgb

# Format nodes by values
nodes <- data.frame(t$data)
#nodes <- nodes[!nodes$isTip,]
nodes$label <- as.numeric(nodes[nodes$label,])
as.numeric(nodes$label)
#nodes$label[nodes$label < 0.80] <- ''
nodes$support[nodes$isTip] <- 'supported'
nodes$support[(!nodes$isTip) & (nodes$label > 0.80)] <- 'supported'
nodes$support[(!nodes$isTip) & (nodes$label < 0.80)] <- 'unsupported'
nodes$support[(!nodes$isTip) & (nodes$label == '')] <- 'supported'
t <- t + aes(linetype=nodes$support)
nodes$label[nodes$label > 0.80] <- ''
t <- t + geom_nodelab(data=nodes, size=2, hjust=-0.05) # colours as defined by col2rgb


# Annotate a clade with a bar line
# t <- t + geom_cladelabel(node=42, label='sect. Alternaria', align=T, colour='black', offset=-1.5)
# t <- t + geom_cladelabel(node=70, label='gaisen clade', align=T, colour='black', offset=-4.5)
# t <- t + geom_cladelabel(node=51, label='tenuissima clade', align=T, colour='black', offset=-4.5)
# t <- t + geom_cladelabel(node=45, label='arborescens clade', align=T, colour='black', offset=-4.5)
# t <- t + geom_cladelabel(node=43, label='sect. Alternaria', align=T, colour='black', offset=-0.0)
# t <- t + geom_cladelabel(node=70, label='gaisen clade', align=T, colour='black', offset=-2.0)
# t <- t + geom_cladelabel(node=46, label='tenuissima clade', align=T, colour='black', offset=-2.0)
# t <- t + geom_cladelabel(node=65, label='arborescens clade', align=T, colour='black', offset=-2.0)

# Save as PDF and force a 'huge' size plot
t <- ggsave("tree5.pdf", width =80, height = 30, units = "cm", limitsize = FALSE)
```
