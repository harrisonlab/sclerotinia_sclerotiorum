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

## Copy over all busco outputs?
