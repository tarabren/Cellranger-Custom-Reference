# Cellranger Custom Reference
This repo provides a workflow for adding a **custom gene** (e.g., a transgene or fusion gene) to a single-cell RNA-seq reference for use with [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome).

---

## Overview

This workflow guides you through:

1. Extracting a custom sequence from a SnapGene `.dna` file
2. Aligning the custom sequence to the reference genome using `minimap2`
3. Parsing the alignment to generate GTF annotations
4. Adding custom FASTA and GTF entries to a Cell Ranger reference
5. Rebuilding the reference using `cellranger mkref`
6. Running alignment and quantification with `cellranger count`

---

## 🛠 Requirements

- Python ≥ 3.7
- [minimap2](https://github.com/lh3/minimap2)
- [samtools](http://www.htslib.org/)
- [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation) (tested with v6+)
- SnapGene `.dna` file (optional, can substitute with FASTA)

---

## File Descriptions

### 🧬 Input Files
- `your_construct.dna`: SnapGene format for the custom gene
- `reference_genome/`: Directory for reference genome
    - `reference_genome/fasta/genome.fa`: The reference FASTA (from 10x Genomics)
    - `reference_genome/genes/genes.gtf`: The reference GTF file (from 10x Genomics)

### ⚙️ Scripts

| Script | Description | Output |
|--------|-------------|--------|
| `scripts/extract_fasta_from_snapgene.py` | Extracts a DNA sequence from a SnapGene `.dna` file into FASTA format | custom_sequence.fasta |
| `scripts/parse_sam.py` | Parses a SAM file to identify exons and creates a GTF annotation | custom_genes.gtf |
| `scripts/run_mkref.sh` | Builds a custom Cell Ranger reference using modified FASTA and GTF | Custom_Ref/ |
---

## Workflow

1. **Extract custom sequence from snap file:**

    ```bash
    module load python/3.12.2
    python scripts/extract_fasta_from_snapgene.py
    ```

2. **Add Custom Sequence to genome.fa File:**

    ```bash
    cat custom_sequence.fasta >> reference_genome/fasta/genome.fa
    ```

3. **Align with minimap2:**

    ```bash
    minimap2 -ax splice reference_genome/fasta/genome.fa custom_sequence.fasta > custom_aligned.sam
    ```

    Check output:
     ```bash
    module load samtools
    samtools view custom_aligned.sam | head -n 10
    ```
4. **Parse alignment and add to GTF:**

    ```python
    python scripts/parse_sam.py --gene_name your_gene_name --transcript_id your_transcript_id
    ```

    If your reference file is in .gz format, run the following to decompress it:
    ```bash
    gunzip reference_genome/genes/genes.gtf.gz reference_genome/genes/genes.gtf
    ```

    Now add the custom gtf file to the reference
    ```bash
    cat custom_genes.gtf >> reference_genome/genes/genes.gtf
    ```

5. **Rebuild the reference:**

    ```bash
    bash scripts/run_mkref.sh
    ```

You're custom reference has been added! 


