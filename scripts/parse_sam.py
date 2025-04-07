import re
import pysam
import argparse

parser = argparse.ArgumentParser(description="Parse a SAM file and generate a GTF for a custom gene.")
parser.add_argument("--gene_name", required=True, help="Name of the custom gene (e.g., GFP_abc)")
parser.add_argument("--transcript_id", required=True, help="Transcript ID for the gene")
args = parser.parse_args()

# input and output file names
sam_file = "custom_aligned.sam"
gtf_file = "custom_genes.gtf"

gene_name = args.gene_name
transcript_id = args.transcript_id

# Function to parse the CIGAR string and extract exon positions
def parse_cigar(start, cigar):
    exons = []
    pos = start
    exon_start = None
    cigar_tuples = re.findall(r'(\d+)([MIDNSHP=X])', cigar)

    for length, op in cigar_tuples:
        length = int(length)
        if op == 'M':
            if exon_start is None:
                exon_start = pos
            pos += length
        elif op == 'N':
            if exon_start is not None:
                exons.append((exon_start, pos - 1))
                exon_start = None
            pos += length
        elif op in ['D', 'I']:
            pass
    
    if exon_start is not None:
        exons.append((exon_start, pos - 1))
    
    return exons

# Read the SAM file and write the GTF
samfile = pysam.AlignmentFile(sam_file, "r")
with open(gtf_file, "w") as gtf_out:
    for read in samfile.fetch():
        chrom = read.reference_name
        start = read.reference_start + 1
        cigar = read.cigarstring
        strand = "-" if read.flag & 16 else "+"

        exons = parse_cigar(start, cigar)
        end = exons[-1][1] if exons else start

        gtf_out.write(f"{chrom}\tCustom\tgene\t{start}\t{end}\t.\t{strand}\t.\tgene_id \"{gene_name}\"; gene_name \"{gene_name}\"; gene_version \"1\"; gene_type \"protein_coding\";\n")
        gtf_out.write(f"{chrom}\tCustom\ttranscript\t{start}\t{end}\t.\t{strand}\t.\tgene_id \"{gene_name}\"; transcript_id \"{transcript_id}\"; transcript_version \"1\"; transcript_type \"protein_coding\"; transcript_name \"{gene_name}\";\n")

        for i, (exon_start, exon_end) in enumerate(exons, start=1):
            gtf_out.write(f"{chrom}\tCustom\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\tgene_id \"{gene_name}\"; transcript_id \"{transcript_id}\"; exon_number \"{i}\"; exon_id \"{gene_name}_exon{i}\"; exon_version \"1\";\n")

print(f"GTF file generated: {gtf_file}")
