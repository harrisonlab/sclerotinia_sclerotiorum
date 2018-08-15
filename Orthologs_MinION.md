# Orthology analysis between 12 Alternaria spp. isolates

```bash
  IsolateAbrv=Scl_all_isolates
  WorkDir=analysis/orthology/orthomcl/MinION_genomes/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  
```

# 1.a) Format fasta files

##for S. sclerotiorum isolate P7
```bash
Taxon_code=Sscl_1
Fasta_file=$(ls gene_pred/final/MinION_genomes/S.sclerotiorum/P7/final/final_genes_appended_renamed.pep.fasta)
Id_field=1
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

##for S. minor isolate S5
```bash
Taxon_code=Smin_1
Fasta_file=$(ls gene_pred/final/MinION_genomes/S.minor/S5/final/final_genes_appended_renamed.pep.fasta)
Id_field=1
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

##for S. subarctica isolate HE1
```bash
Taxon_code=Ssub_1
Fasta_file=$(ls gene_pred/final/MinION_genomes/S.subarctica/HE1/final/final_genes_appended_renamed.pep.fasta)
Id_field=1
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

##for S. sclerotiorum publlished
```bash
Taxon_code=Sscl_2
Fasta_file=$(ls analysis/orthology/orthomcl/MinION_genomes/S_sclerotiorum_1980_edit.fasta)
Id_field=1
orthomclAdjustFasta $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```

#2) Filter proteins into good and poor sets.
```bash
Input_dir=$WorkDir/formatted
Min_length=10
Max_percent_stops=20
Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
orthomclFilterFasta $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file
```

```bash 
IsolateAbrv=Scl_all_isolates
WorkDir=analysis/orthology/orthomcl/MinION_genomes/$IsolateAbrv
orthofinder -f $WorkDir/formatted -t 3 -a 3
```
