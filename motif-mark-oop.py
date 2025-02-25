#!/usr/bin/env python 

import argparse 
import cairo 
import os
import re

def get_args():
    parser = argparse.ArgumentParser(description="Identify and count specific motifs given a FASTA file < 1000 bp")
    parser.add_argument("-f", help="FASTA file containing motifs, exons are capitalized", required=True, type=str)
    parser.add_argument("-m", help="Motif file containing motifs to search for, one motif per line", required=True, type=str)
    return parser.parse_args()


 # Keys are nucleotides codes and values are possible bases

iupac_dict = {"A":'[Aa]',"C":'[Cc]',"G":'[Gg]',"T":'[TUtu]',"U":'[TUtu]',"R":'[AGag]',"Y":'[CTct]',"S":'[GCgc]',"W":'[ATat]',"K":'[GTgt]'
            ,"M":'[ACac]',"B":'[CGTcgt]',"D":'[AGTagt]',"H":'[ACTact]',"V":'[ACGacg]',"N":'[AaCcTtUuGg]',"a":'[Aa]',"c":'[Cc]',"g":'[Gg]'
            ,"t":'[TUtu]',"u":'[TUtu]',"r":'[AGag]',"y":'[CTct]',"s":'[GCgc]',"w":'[ATat]',"k":'[GTgt]',"m":'[ACac]',"b":'[CGTcgt]'
            ,"d":'[AGTagt]',"h":'[ACTact]',"v":'[ACGacg]',"n":'[AaCcTtUuGg]'}



class ParseFasta():

    def __init__(self, file_path):
        with open(file_path, 'r') as file:
            lines = file.readlines()
            self.gene = lines[0].split(" ")[0].strip(">")
            self.seq = "".join([line.strip('\n') for line in lines if not line.startswith(">")])
            self.seq_length = len(self.seq) 

    def find_exon_positions(self):
        exon_st = None
        exons = []  # List to store exon (start, end) tuples
        in_exon = False

        for i, ch in enumerate(self.seq):
            if ch.isupper():
                if not in_exon:  # If we're not already in an exon, this is a new exon start
                    exon_st = i + 1  # 1-based index
                    in_exon = True
                exon_end = i + 1
            else:  # Found a lowercase character (intron)
                if in_exon:  # The exon ends here
                    exons.append((exon_st, exon_end))
                    in_exon = False  # Reset the flag

        # If sequence ends while still in an exon, append the last exon
        if in_exon:
            exons.append((exon_st, exon_end))

        return exons

    def find_motifs(self, motifs):
    motif_positions = {}
    for motif in motifs:
        matches = []
        # Loop through the sequence with a sliding window of motif length
        for i in range(len(self.seq) - len(motif) + 1):
            if re.match(motif, self.seq[i:i+len(motif)], re.IGNORECASE):
                matches.append((i + 1, i + len(motif) + 1))
        motif_positions[motif] = matches
    return motif_positions
                    
class ParseMotif: 
    def __init__(self, file_path):
        with open(file_path, 'r') as file:
            self.motifs = [line.strip() for line in file.readlines()]

    def translate_motifs(self):
        translated_motifs = []
        for motif in self.motifs:
            translated_motif = ""
            for ch in motif:
                if ch in iupac_dict:
                    translated_motif += iupac_dict[ch]
                else:
                    raise ValueError(f"Invalid IUPAC code: {ch}")
            translated_motifs.append(translated_motif)
        return translated_motifs





def main():
    args = get_args()

    # Parse motifs
    motif_parser = ParseMotif(args.m)
    translated_motifs = motif_parser.translate_motifs()

    # Parse FASTA
    fasta_parser = ParseFasta(args.f)
    exons = fasta_parser.find_exon_positions()
    motif_positions = fasta_parser.find_motifs(translated_motifs)

    # Visualization (pseudo-code)
    # 1. Initialize cairo surface and context
    # 2. Draw the sequence as a horizontal line
    # 3. Highlight exons
    # 4. Overlay motifs with unique colors
    # 5. Save or display the image

if __name__ == "__main__":
    main()