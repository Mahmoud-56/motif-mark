#!/usr/bin/env python 

import argparse as ap
import cairo 
import math 

def get_args():
    p = ap.ArgumentParser(description = "Indentify and count specefic motifs given a FASTA file < 1000 bp")
    p.add_argument("-f", help = "FASTA file containing motids, exons are capitalized", required = True)
    p.add_argument("-m", help = "Motif to search for in the FASTA file, one motif per line", required = False)
    return p.parse_args()


class ParseSeq():

    def __init__(self, lines):
        self.gene = lines[0].split(" ")[0].strip(">")
        self.seq = "".join([line.strip('\n') for line in lines if not line.startswith(">")])

    def seq_length(self):
        return int(len(self.seq))

    def iupac_name(motifs_files):

        # Keys are nucleotides codes and values are possible bases

        iupac = {"A":'[Aa]',"C":'[Cc]',"G":'[Gg]',"T":'[TUtu]',"U":'[TUtu]',"R":'[AGag]',"Y":'[CTct]',"S":'[GCgc]',"W":'[ATat]',"K":'[GTgt]',"M":'[ACac]',"B":'[CGTcgt]',"D":'[AGTagt]',"H":'[ACTact]',"V":'[ACGacg]',"N":'[A-Za-z]',
                              "a":'[Aa]',"c":'[Cc]',"g":'[Gg]',"t":'[TUtu]',"u":'[TUtu]',"r":'[AGag]',"y":'[CTct]',"s":'[GCgc]',"w":'[ATat]',"k":'[GTgt]',"m":'[ACac]',"b":'[CGTcgt]',"d":'[AGTagt]',"h":'[ACTact]',"v":'[ACGacg]',"n":'[A-Za-z]'}



args = get_args()

Fasta = open(args.f,'r')

lines = Fasta.readlines()

seq = ParseSeq(lines)



