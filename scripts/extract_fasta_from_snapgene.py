from Bio import SeqIO

# Define file paths
input_file = "your_construct.dna"  # Your SnapGene file
output_fasta = "custom_sequence.fasta"  # Output FASTA file

# Read the SnapGene file and extract the sequence
with open(input_file, "rb") as file:
    record = SeqIO.read(file, "snapgene")  # Read as SnapGene format

# Convert the sequence to uppercase
record.seq = record.seq.upper()

# Save as FASTA
with open(output_fasta, "w") as fasta_file:
    fasta_file.write(f">{record.id}\n{record.seq}\n")

print(f"Extracted sequence saved to {output_fasta} (uppercase format)")
