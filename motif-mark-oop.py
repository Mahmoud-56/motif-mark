#!/usr/bin/env python 

import argparse
import cairo
import os
import re
from collections import defaultdict


def get_args():
    parser = argparse.ArgumentParser(description="Identify and visualize motifs in a FASTA files with multiple genes")
    parser.add_argument("-f", help="FASTA file containing exons capitalizes and introns as lowercase", required=True, type=str)
    parser.add_argument("-m", help="Motif file with one motif per line", required=True, type=str)
    return parser.parse_args()




# IUPAC codes for motif expansion
iupac_dict = {
    "A": '[Aa]', "C": '[Cc]', "G": '[Gg]', "T": '[TUtu]', "U": '[TUtu]',
    "R": '[AGag]', "Y": '[CTct]', "S": '[GCgc]', "W": '[ATat]', "K": '[GTgt]',
    "M": '[ACac]', "B": '[CGTcgt]', "D": '[AGTagt]', "H": '[ACTact]', "V": '[ACGacg]', "N": '[AaCcTtUuGg]'
}

class ParseFasta:
    def __init__(self, file_path):
        with open(file_path, 'r') as file:
            lines = file.readlines()
            self.gene = lines[0].split(" ")[0].strip(">") #extract the gene 
            self.seq = "".join(line.strip() for line in lines if not line.startswith(">")) #extract the seq
            self.seq_length = len(self.seq) #extract the seq length 

    def find_exon_positions(self):
        """Identify exon start and end positions based on uppercase sequences."""
        exons = []
        exon_start = None
        in_exon = False #starting in an intron

        for i, ch in enumerate(self.seq):
            if ch.isupper():
                if not in_exon: #new exon start 
                    exon_start = i + 1  # 1-based index
                    in_exon = True
                exon_end = i + 1
            else: #lower case letter 
                if in_exon:
                    exons.append((exon_start, exon_end))
                    in_exon = False #reset the flag
        if in_exon:
            exons.append((exon_start, exon_end))

        return exons

    def find_motifs(self, motifs_dict):
        """
        Find motif positions in the sequence using regex from translated IUPAC motifs.

        Input:
            original motif dict where keys are motifs and values are respective regex expressions

        Output:
            a dictionary where keys are original motifs, and values are start and end positions.

        """

        motif_positions = defaultdict(list) # the default value for any new key is an empty list 
        for orig_motif, regex in motifs_dict.items():
            for match in re.finditer(regex, self.seq, re.IGNORECASE):
                motif_positions[orig_motif].append((match.start() + 1, match.end()))

        return motif_positions


class ParseMotif:
    def __init__(self, file_path):
        with open(file_path, 'r') as file:
            self.motifs = [line.strip() for line in file.readlines()]

    def translate_motifs(self):
        """Convert IUPAC motifs to regex while keeping original motifs for legends."""
        translated_motifs = {}
        for motif in self.motifs:
            for ch in motif:
                regex = "".join(iupac_dict.get(ch,ch)) #if value not found add the original motif
                translated_motifs[motif] = regex  # Store original motif as key, regex as value
        return translated_motifs

def draw_visualization(fasta_parser, motifs_dict, output_file="Figure_1.png"):
    """Draws a single image with multiple genes, labeled, with motifs and exons."""
    
    # Constants for drawing
    sequence_width = 800
    gene_spacing = 100
    legend_width = 200
    image_width = sequence_width + legend_width
    image_height = len(fasta_parser) * gene_spacing + 100

    # Create surface and context
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, image_width, image_height)
    ctx = cairo.Context(surface)

    # Set background to white
    ctx.set_source_rgb(1.0, 1.0, 1.0)
    ctx.paint()

    # Assign colors to motifs
    motif_colors = {motif: (i / len(motifs_dict), 0, 1 - i / len(motifs_dict))
                    for i, motif in enumerate(motifs_dict.keys())}

    # Draw each gene
    for i, fasta_parser in enumerate(fasta_parser):
        y_offset = i * gene_spacing + 50
        scale = sequence_width / fasta_parser.seq_length

        # Draw gene name
        ctx.set_source_rgb(0, 0, 0)
        ctx.move_to(10, y_offset + 5)
        ctx.show_text(fasta_parser.gene)
        ctx.stroke()

    # Draw intron line
        ctx.set_source_rgb(0, 0, 0)
        ctx.set_line_width(1)
        ctx.move_to(50, y_offset)
        ctx.line_to(50 + sequence_width, y_offset)
        ctx.stroke()

    # Draw exons
    for start, end in fasta_parser.find_exon_positions():
        exon_start = 50 + (start - 1) * scale
        exon_length = (end - start + 1) * scale
        ctx.rectangle(exon_start, y_offset - 10, exon_length, 20)
        ctx.fill()

    # Draw motifs
        ctx.set_line_width(2)
        for motif, positions in fasta_parser.find_motifs(motifs_dict).items():
            ctx.set_source_rgb(*motif_colors[motif])
            for start, _ in positions:
                tick_x = 50 + (start - 1) * scale
                ctx.move_to(tick_x, y_offset - 20)
                ctx.line_to(tick_x, y_offset + 20)
                ctx.stroke()


    # Draw legend
    legend_x = sequence_width + 60
    legend_y = 50
    ctx.set_font_size(12)

    for i, (motif, color) in enumerate(motif_colors.items()):
        ctx.set_source_rgb(*color)
        ctx.rectangle(legend_x, legend_y + i * 20, 10, 10)
        ctx.fill()

        ctx.set_source_rgb(0, 0, 0)
        ctx.move_to(legend_x + 15, legend_y + i * 20 + 10)
        ctx.show_text(motif)
        ctx.stroke()


    # Save image
    surface.write_to_png(output_file)
    print(f"Visualization saved as {output_file}")


def main():
    args = get_args()
    motif_parser = ParseMotif(args.m)
    motifs_dict = motif_parser.translate_motifs()  # Dictionary {original motif: regex}

    fasta_parser =ParseFasta(args.f)

    draw_visualization(fasta_parser, motifs_dict, output_file="Figure_1.png")

if __name__ == "__main__":
    main()


