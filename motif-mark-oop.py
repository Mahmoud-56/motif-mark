#!/usr/bin/env python3
import argparse
import cairo
import re
from collections import defaultdict


def get_args():
    parser = argparse.ArgumentParser(description="Identify and visualize motifs in a FASTA file with multiple genes")
    parser.add_argument("-f", help="FASTA file containing exons capitalized and introns as lowercase", required=True, type=str)
    parser.add_argument("-m", help="Motif file with one motif per line", required=True, type=str)
    return parser.parse_args()


iupac_dict = {'A': '[Aa]', 'a': '[Aa]',
               'C': '[Cc]', 'c': '[Cc]',
               'T': '[Tt]', 't': '[Tt]', 
               'G': '[Gg]', 'g': '[Gg]',
               'U': '[UuTt]', 'u': '[UuTt]',
               'W': '[AaTt]', 'w': '[AaTt]',
               'S': '[CcGg]', 's': '[CcGg]',
               'M': '[AaCc]', 'm': '[AaCc]',
               'K': '[GgTt]', 'k': '[GgTt]',
               'R': '[AaGg]', 'r': '[AaGg]',
               'Y': '[CcTt]', 'y': '[CcTt]',
               'B': '[CcGgTt]', 'b': '[CcGgTt]',
               'D': '[AaGgTt]', 'd': '[AaGgTt]',
               'H': '[AaCcTt]', 'h': '[AaCcTt]',
               'V': '[AaCcGg]', 'v': '[AaCcGg]',
               'N': '[AaCcTtUuGg]', 'n': '[AaCcTtUuGg]'}



# class ParseMotif:
#     def __init__(self, file_path):
#         with open(file_path, 'r') as file:
#             self.motifs = [line.strip() for line in file.readlines()]

#     def translate_motifs(self):
#         """Convert IUPAC motifs to regex."""
#         for ch in self.motifs:
#             regex = ""
#             regex += iupac_dict[ch]
#             motifs_dict[chr] = regex
#         return motifs_dict

# class ParseFasta:
#     def __init__(self, file_path):
#         self.genes = []
#         with open(file_path, 'r') as file:
#             gene_name = None
#             gene_seq = ""
#             for line in file:
#                 line = line.strip()
#                 if line.startswith(">"):
#                     if gene_name:
#                         self.genes.append((gene_name, gene_seq))
#                     gene_name = line[1:].split(" ")[0]
#                     gene_seq = ""
#                 else:
#                     gene_seq += line
#             if gene_name:
#                 self.genes.append((gene_name, gene_seq))

#     def get_genes(self):
#         return self.genes

#     def find_exon_positions(self, seq):
#         """Identify exon start and end positions based on uppercase letters."""
#         exons = []
#         exon_start = None
#         in_exon = False

#         for i, ch in enumerate(seq):
#             if ch.isupper():
#                 if not in_exon:
#                     exon_start = i + 1
#                     in_exon = True
#                 exon_end = i + 1
#             else:
#                 if in_exon:
#                     exons.append((exon_start, exon_end))
#                     in_exon = False
#         if in_exon:
#             exons.append((exon_start, exon_end))
#         return exons

#     def find_motifs(self, seq, motifs_dict):
#         """Find motif positions in the sequence using regex."""
#         motif_positions = defaultdict(list)
#         for translated_motifs, regex in motifs_dict.items():
#             for match in re.finditer(regex, seq, re.IGNORECASE):
#                 motif_positions[translated_motifs].append((match.start() + 1, match.end()))
#         return motif_positions

class ParseMotif:
    def __init__(self, file_path):
        with open(file_path, 'r') as file:
            self.motifs = [line.strip() for line in file.readlines()]

    def translate_motifs(self):
        """Convert IUPAC motifs to regex."""
        return {motif: "".join(iupac_dict.get(ch, ch) for ch in motif) for motif in self.motifs}

class ParseFasta:
    def __init__(self, file_path):
        self.genes = []
        with open(file_path, 'r') as file:
            gene_name = None
            gene_seq = ""
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    if gene_name:
                        self.genes.append((gene_name, gene_seq))
                    gene_name = line[1:].split(" ")[0]
                    gene_seq = ""
                else:
                    gene_seq += line
            if gene_name:
                self.genes.append((gene_name, gene_seq))

    def get_genes(self):
        return self.genes

    def find_exon_positions(self, seq):
        """Identify exon start and end positions based on uppercase letters."""
        exons = []
        exon_start = None
        in_exon = False

        for i, ch in enumerate(seq):
            if ch.isupper():
                if not in_exon:
                    exon_start = i + 1
                    in_exon = True
                exon_end = i + 1
            else:
                if in_exon:
                    exons.append((exon_start, exon_end))
                    in_exon = False
        if in_exon:
            exons.append((exon_start, exon_end))
        return exons

    def find_motifs(self, seq, motifs_dict):
        """Find motif positions in the sequence using regex."""
        motif_positions = defaultdict(list)
        for motif, regex in motifs_dict.items():
            for match in re.finditer(regex, seq):
                motif_positions[motif].append((match.start() + 1, match.end()))
        return motif_positions

colors = [
    (1, 0, 0),    # Red
    (0, 0, 1),    # Blue
    (0, 1, 0),    # Green
    (1, 0.5, 0),  # Orange
    (0.5, 0, 0.5),# Purple
    (0, 0.5, 0.5),# Teal
    (1, 0, 1),    # Magenta
    (0, 1, 1),    # Cyan
    (0.5, 0.5, 0) # Olive
]

def draw_visualization(fasta_parser, motifs_dict, output_file="Figure_1.png"):
    """Draw a visualization of genes with exons, introns, and motifs."""
    sequence_width = 800
    gene_spacing = 100
    legend_height = 100  # Reduce height for better spacing
    image_width = sequence_width + 200
    image_height = len(fasta_parser.genes) * gene_spacing + legend_height + 150  # Add extra padding

    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, image_width, image_height)
    ctx = cairo.Context(surface)
    ctx.set_source_rgb(1.0, 1.0, 1.0)
    ctx.paint()

    # Draw plot title
    ctx.set_source_rgb(0, 0, 0)
    ctx.set_font_size(36)  
    ctx.move_to(image_width / 2 - 100, 50)  # Center the title
    ctx.show_text("Motif Plot")

    # Assign unique colors to motifs
    motif_colors = {motif: colors[i % len(colors)] for i, motif in enumerate(motifs_dict.keys())}

    for i, (gene_name, seq) in enumerate(fasta_parser.get_genes()):
        y_offset = i * gene_spacing + 150
        scale = sequence_width / len(seq)

        # Adjust gene name position
        ctx.set_source_rgb(0, 0, 0)
        ctx.set_font_size(20)  
        ctx.move_to(40, y_offset + 10)  # Move text slightly left and down
        ctx.show_text(gene_name)

        # Draw intron line
        ctx.set_line_width(4)
        ctx.move_to(150, y_offset)  
        ctx.line_to(150 + sequence_width, y_offset)
        ctx.stroke()

        # Draw exons
        for start, end in fasta_parser.find_exon_positions(seq):
            exon_start = 150 + (start - 1) * scale
            exon_length = (end - start + 1) * scale
            ctx.rectangle(exon_start, y_offset - 10, exon_length, 20)
            ctx.fill()

        # Draw motifs
        ctx.set_line_width(4)
        for motif, positions in fasta_parser.find_motifs(seq, motifs_dict).items():
            ctx.set_source_rgb(*motif_colors[motif])
            for start, _ in positions:
                tick_x = 150 + (start - 1) * scale
                ctx.move_to(tick_x, y_offset - 15)
                ctx.line_to(tick_x, y_offset + 15)
                ctx.stroke()

    # Draw legend at the bottom
    legend_y = image_height - 60  
    legend_x_start = 100

    ctx.set_font_size(18)
    ctx.set_source_rgb(0, 0, 0)
    ctx.move_to(legend_x_start, legend_y - 10)
    ctx.show_text("Legend:")

    x_offset = legend_x_start + 80
    for motif, color in motif_colors.items():
        ctx.set_source_rgb(*color)
        ctx.set_line_width(5)
        ctx.move_to(x_offset, legend_y)
        ctx.line_to(x_offset + 30, legend_y)
        ctx.stroke()

        ctx.set_source_rgb(0, 0, 0)
        ctx.move_to(x_offset + 40, legend_y + 5)
        ctx.show_text(motif)

        x_offset += 150  # Space out legend items

    # Save image
    surface.write_to_png(output_file)
    print(f"Visualization saved as {output_file}")

def main():
    args = get_args()
    motif_parser = ParseMotif(args.m)
    motifs_dict = motif_parser.translate_motifs()
    fasta_parser = ParseFasta(args.f)
    draw_visualization(fasta_parser, motifs_dict, output_file="Figure_1.png")

if __name__ == "__main__":
    main()
