#!/usr/bin/env python 

import argparse as ap
import cairo 
import math 

def get_args():
    p = ap.ArgumentParser(description = "Indentify and count specefic motifs given a FASTA file < 1000 bp")
    p.add_argument("-f", help = "FASTA file containing motids, exons are capitalized", required = True)
    p.add_argument("-m", help = "Motif to search for in the FASTA file, one motif per line", required = True)
    return p.parse_args()


