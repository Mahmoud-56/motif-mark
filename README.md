# Motif Mark: A Visualization Tool

This tool visualizes gene sequences with exons, introns, and motifs. It parses a FASTA file, detects exons (uppercase regions), and highlights motif locations
using IUPAC dictionary. The output is a PNG visualization of each gene. 

### Input Files
 - FASTA file: Contains gene sequences (maximum 1000 bases) where exons are uppercase. 
 - Motifs file: List of motifs of 10 bases of less. One per line in a text file 

### Output
- Figure_1.png: Visualization of genes, exons, and motifs.

### Command 

```
./motif-mark-oop.py -f <fasta file> -m <motif file>
```


