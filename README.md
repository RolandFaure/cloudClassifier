# cloudClassifier
Improves the taxonomic assignment of reads using barcode information of linked-reads data. The formats are built to be compatible with [Kraken2](github.com/DerrickWood/kraken2), but any classifier can be used with a little reformatting.

## Dependancies

1. A working pyhton3 installation
2. Standard unix commands `awk` and `paste`
3. [ETE toolkit](http://etetoolkit.org/download/), intallable through conda
4. The first time you run cloudClassifier, you'll need an internet connexion to download the NCBI taxonomy

## Input

cloudClassifier needs two files as input: 
1. The fastq of the linked reads, with barcodes marked with "BX:Z:" or "BC:Z:" tag. This fastq MUST be sorted by barcode. You can use the `sort` command on unix systems (with option `-S` for limiting memory usage)
2. The taxonomic classification of those reads, presented in the same order as the fastq file. Use the default [Kraken2](github.com/DerrickWood/kraken2) output

## Output

cloudClassifier outputs a file in [Kraken2](github.com/DerrickWood/kraken2) default output format, with improved taxonomic assignment.
