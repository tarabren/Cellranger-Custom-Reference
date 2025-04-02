# Cellranger-Custom-Reference
This repo provides a workflow for adding a **custom gene** (e.g., a transgene or fusion gene) to a single-cell RNA-seq reference for use with [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome).

---

## 📦 Overview

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

## 🗂 File Descriptions

### 🧬 Input Files
- `your_construct.dna`: SnapGene format for the custom gene
- `reference_genome.fa`: The reference FASTA (from 10x Genomics)
- `genes.gtf`: The reference GTF file (from 10x Genomics)

### ⚙️ Scripts

| Script | Description |
|--------|-------------|
| `scripts/extract_fasta_from_snapgene.py` | Extracts a DNA sequence from a SnapGene `.dna` file into FASTA format |
| `scripts/parse_sam_to_gtf.py` | Parses a SAM file to identify exons and creates a GTF annotation |
| `scripts/run_mkref.sh` | Builds a custom Cell Ranger reference using modified FASTA and GTF |
| `scripts/run_counts.sh` | Runs `cellranger count` with the new custom reference |

---

## 🚀 Quick Start

1. **Extract custom sequence:**

    ```bash
    python scripts/extract_fasta_from_snapgene.py --input your_construct.dna --output custom_sequence.fasta
    ```

2. **Align with minimap2:**

    ```bash
    minimap2 -ax splice reference_genome.fa custom_sequence.fasta > custom_aligned.sam
    ```

3. **Parse alignment to GTF:**

    ```bash
    python scripts/parse_sam_to_gtf.py --input custom_aligned.sam --output custom_genes.gtf
    ```

4. **Update reference files:**

    Append the custom sequence to your reference FASTA:

    ```bash
    cat custom_sequence.fasta >> path_to_cellranger_reference/fasta/genome.fa
    ```

    Decompress and append the custom GTF:

    ```bash
    gunzip path_to_cellranger_reference/genes/genes.gtf.gz
    cat custom_genes.gtf >> path_to_cellranger_reference/genes/genes.gtf
    ```

5. **Rebuild the reference:**

    ```bash
    bash scripts/run_mkref.sh
    ```

6. **Run Cell Ranger count:**

    ```bash
    bash scripts/run_counts.sh
    ```



