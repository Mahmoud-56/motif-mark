#!/usr/bin/env python 

import argparse 
import cairo 
import os
import re




# class ParseFasta():

#     def __init__(self, file_path):
#         with open(file_path, 'r') as file:
#             lines = file.readlines()
#             self.gene = lines[0].split(" ")[0].strip(">")
#             self.seq = "".join([line.strip('\n') for line in lines if not line.startswith(">")])
#             self.seq_length = len(self.seq) 

#     def find_exon_positions(self):
#         exon_st = None
#         exons = []  # List to store exon (start, end) tuples
#         in_exon = False

#         for i, ch in enumerate(self.seq):
#             if ch.isupper():
#                 if not in_exon:  # If we're not already in an exon, this is a new exon start
#                     exon_st = i + 1  # 1-based index
#                     in_exon = True
#                 exon_end = i + 1
#             else:  # Found a lowercase character (intron)
#                 if in_exon:  # The exon ends here
#                     exons.append((exon_st, exon_end))
#                     in_exon = False  # Reset the flag

#         # If sequence ends while still in an exon, append the last exon
#         if in_exon:
#             exons.append((exon_st, exon_end))

#         return exons

#     def find_motifs(self, motifs):
#         motif_positions = {}
#         for motif in motifs:
#             matches = []
#         # Loop through the sequence with a sliding window of motif length
#         for i in range(len(self.seq) - len(motif) + 1):
#             if re.match(motif, self.seq[i:i+len(motif)], re.IGNORECASE):
#                 matches.append((i + 1, i + len(motif) + 1))
#         motif_positions[motif] = matches
#         return motif_positions
                    
# class ParseMotif: 
#     def __init__(self, file_path):
#         with open(file_path, 'r') as file:
#             self.motifs = [line.strip() for line in file.readlines()]

#     def translate_motifs(self):
#         translated_motifs = []
#         for motif in self.motifs:
#             translated_motif = ""
#             for ch in motif:
#                 if ch in iupac_dict:
#                     translated_motif += iupac_dict[ch]
#                 else:
#                     raise ValueError(f"Invalid IUPAC code: {ch}")
#             translated_motifs.append(translated_motif)
#         return translated_motifs

#     def draw_visualization(fasta_parser, exons, motif_positions, output_file="Figure_1.png"):
#     # Canvas setup
#         width, height = 800, 200
#         surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
#         ctx = cairo.Context(surface)

#     # Scaling factor
#         scale = width / fasta_parser.seq_length

#     # Draw the full sequence as a thin line (introns)
#         ctx.set_source_rgb(0.5, 0.5, 0.5)  # Gray for introns
#         ctx.set_line_width(1)
#         ctx.move_to(0, height // 2)
#         ctx.line_to(width, height // 2)
#         ctx.stroke()

#     # Draw exons as black filled rectangles
#         ctx.set_source_rgb(0.0, 0.0, 0.0)  # Black for exons
#         for start, end in exons:
#             exon_start = (start - 1) * scale  # Convert to 0-based pixel coordinates
#             exon_length = (end - start + 1) * scale
#             ctx.rectangle(exon_start, height // 2 - 10, exon_length, 20)  # Height of 20 pixels
#             ctx.fill()

#     # Draw motifs (unique colors, overlaid on sequence/exons)
#         motif_colors = {
#             motif: (i / len(motif_positions), 0, 1 - i / len(motif_positions))  # Simple color gradient
#             for i, motif in enumerate(motif_positions.keys())
#         }
#         for motif, positions in motif_positions.items():
#             ctx.set_source_rgb(*motif_colors[motif])
#             for start, end in positions:
#                 motif_start = (start - 1) * scale
#                 motif_length = (end - start + 1) * scale
#                 ctx.rectangle(motif_start, height // 2 - 5, motif_length, 10)  # Slightly smaller height
#                 ctx.fill()

#         # Save to file
#         surface.write_to_png(output_file)
#         print(f"Visualization saved as {output_file}")



# def main():
#     args = get_args()
#     motif_parser = ParseMotif(args.m)
#     translated_motifs = motif_parser.translate_motifs()
#     fasta_parser = ParseFasta(args.f)
#     exons = fasta_parser.find_exon_positions()
#     motif_positions = fasta_parser.find_motifs(translated_motifs)
#     draw_visualization(fasta_parser, exons, motif_positions)

# if __name__ == "__main__":
#     main()






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

class ParseFasta:
    def __init__(self, file_path):
        with open(file_path, 'r') as file:
            lines = file.readlines()
            self.gene = lines[0].split(" ")[0].strip(">")
            self.seq = "".join([line.strip('\n') for line in lines if not line.startswith(">")])
            self.seq_length = len(self.seq)

    def find_exon_positions(self):
        exon_st = None
        exons = []
        in_exon = False
        for i, ch in enumerate(self.seq):
            if ch.isupper():
                if not in_exon:
                    exon_st = i + 1
                    in_exon = True
                exon_end = i + 1
            else:
                if in_exon:
                    exons.append((exon_st, exon_end))
                    in_exon = False
        if in_exon:
            exons.append((exon_st, exon_end))
        return exons

    def find_motifs(self, motifs):
        motif_positions = {}
        for motif in motifs:
            matches = []
            for i in range(len(self.seq) - len(motif) + 1):
                if re.match(motif, self.seq[i:i+len(motif)]):
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

def draw_visualization(fasta_parser, exons, motif_positions, output_file="Figure_1.png"):
    # Canvas setup: Increased width for legend
    total_width, height = 1000, 200  # 800 for sequence, 200 for legend
    sequence_width = 800
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, total_width, height)
    ctx = cairo.Context(surface)
    scale = sequence_width / fasta_parser.seq_length

    # Draw the full sequence as a gray line (introns)
    ctx.set_source_rgb(0.5, 0.5, 0.5)  # Gray for introns
    ctx.set_line_width(1)
    ctx.move_to(0, height // 2)
    ctx.line_to(sequence_width, height // 2)
    ctx.stroke()

    # Draw exons as black filled rectangles
    ctx.set_source_rgb(0.0, 0.0, 0.0)  # Black for exons
    for start, end in exons:
        exon_start = (start - 1) * scale
        exon_length = (end - start + 1) * scale
        ctx.rectangle(exon_start, height // 2 - 10, exon_length, 20)
        ctx.fill()

    # Draw motifs as vertical ticks above the sequence
    motif_colors = {motif: (i / len(motif_positions), 0, 1 - i / len(motif_positions)) 
                    for i, motif in enumerate(motif_positions.keys())}
    ctx.set_line_width(2)  # Thicker lines for visibility
    for motif, positions in motif_positions.items():
        ctx.set_source_rgb(*motif_colors[motif])
        for start, _ in positions:  # Only use start position for tick
            tick_x = (start - 1) * scale
            ctx.move_to(tick_x, height // 2 - 20)  # Start above sequence
            ctx.line_to(tick_x, height // 2 - 40)  # Extend upward
            ctx.stroke()

    # Draw legend on the right
    legend_x = sequence_width + 20  # Start legend after sequence
    legend_y = 20  # Start near top
    ctx.set_font_size(12)
    for i, (motif, color) in enumerate(motif_colors.items()):
        # Draw color square
        ctx.set_source_rgb(*color)
        ctx.rectangle(legend_x, legend_y + i * 20, 10, 10)
        ctx.fill()
        # Draw motif text
        ctx.set_source_rgb(0, 0, 0)  # Black text
        ctx.move_to(legend_x + 15, legend_y + i * 20 + 10)
        ctx.show_text(motif)  # Display the translated motif regex
        ctx.stroke()

    # Save to file
    surface.write_to_png(output_file)
    print(f"Visualization saved as {output_file}")

def main():
    args = get_args()
    motif_parser = ParseMotif(args.m)
    translated_motifs = motif_parser.translate_motifs()
    fasta_parser = ParseFasta(args.f)
    exons = fasta_parser.find_exon_positions()
    motif_positions = fasta_parser.find_motifs(translated_motifs)
    draw_visualization(fasta_parser, exons, motif_positions)

if __name__ == "__main__":
    main()