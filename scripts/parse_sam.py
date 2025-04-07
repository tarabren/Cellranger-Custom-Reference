import re
import pysam

# Input and output file paths
sam_file = "custom_aligned.sam"
gtf_file = "custom_genes.gtf"
gene_name = "Custom Gene"
transcript_id = ""

# Function to parse the CIGAR string and extract exon positions
def parse_cigar(start, cigar):
    exons = []
    pos = start
    exon_start = None

    # Regular expression to extract CIGAR operations
    cigar_tuples = re.findall(r'(\d+)([MIDNSHP=X])', cigar)

    for length, op in cigar_tuples:
        length = int(length)
        
        if op == 'M':  # Matches (exons)
            if exon_start is None:
                exon_start = pos
            pos += length
        elif op == 'N':  # Skipped region (introns)
            if exon_start is not None:
                exons.append((exon_start, pos - 1))  # Store previous exon
                exon_start = None  # Reset exon start
            pos += length  # Move past intron
        elif op in ['D', 'I']:  # Deletions/Insertions (small changes, usually ignored)
            pass  
    
    # Capture last exon if applicable
    if exon_start is not None:
        exons.append((exon_start, pos - 1))

    return exons

# Read the SAM file
samfile = pysam.AlignmentFile(sam_file, "r")
with open(gtf_file, "w") as gtf_out:
    for read in samfile.fetch():
        chrom = read.reference_name
        start = read.reference_start + 1  # SAM format is 0-based, GTF is 1-based
        cigar = read.cigarstring
        strand = "-" if read.flag & 16 else "+"
        
        # Get exon positions from CIGAR
        exons = parse_cigar(start, cigar)
        end = exons[-1][1] if exons else start  # Last exon end position
        
        # Write GTF entries
        gtf_out.write(f"{chrom}\tCustom\tgene\t{start}\t{end}\t.\t{strand}\t.\tgene_id \"{gene_name}\"; gene_name \"{gene_name}\"; gene_version \"1\"; gene_type \"protein_coding\";\n")
        gtf_out.write(f"{chrom}\tCustom\ttranscript\t{start}\t{end}\t.\t{strand}\t.\tgene_id \"{gene_name}\"; transcript_id \"{transcript_id}\"; transcript_version \"1\"; transcript_type \"protein_coding\"; transcript_name \"{gene_name}\";\n")
        
        # Write exon entries
        for i, (exon_start, exon_end) in enumerate(exons, start=1):
            gtf_out.write(f"{chrom}\tCustom\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\tgene_id \"{gene_name}\"; transcript_id \"{transcript_id}\"; exon_number \"{i}\"; exon_id \"{gene_name}_exon{i}\"; exon_version \"1\";\n")

print(f"GTF file generated: {gtf_file}")
