#!/bin/bash
#
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100g

module load cellranger

cellranger mkref --genome=Custom_Ref \
                 --fasta=/reference_genome/fasta/genome.fa \
                 --genes=/reference_genome/genes/genes.gtf

echo "✅ Cell Ranger reference creation completed!"
