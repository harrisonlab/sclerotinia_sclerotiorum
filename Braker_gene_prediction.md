###Gene prediction using Braker

##Gene prediction is trained using RNA-Seq data 

##Downloading data from ncbi
# Data used for this test analysis SRA code: SRR1915981 was downloaded using fastq-dump


```bash
mkdir -p raw_rna/external/S_sclerotiorum_test_RNA_data
fastq-dump -O raw_rna/external/S_sclerotiorum_test_RNA_data SRR1915981

```