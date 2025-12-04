# SimpleSnek: A Simple Snakemake Pipeline for Variant Calling

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.github.io)

This repository contains a Snakemake workflow for conducting variant analysis from DArT fastq. The pipeline takes raw, single-end FASTQ files and performs read trimming, mapping, GATK-based joint-calling, and final filtering to produce high-quality SNP and INDEL variant sets.

This workflow is specifically designed to:
* Process **single-end** sequencing data.
* Handle **polyploid** _or_ **diploid** (`-ploidy 4`) variant calling with GATK.
* Perform **joint-calling** across all samples found in the input directory.
* Use Trimmomatic, BWA, Picard, Samtools, GATK4, and BCFtools.

---

## Table of Contents
* [Pipeline Overview](#pipeline-overview)
* [Required Setup](#required-setup)
* [Configuration](#configuration)
* [Usage](#usage)
* [Output Files](#output-files)

---

## Pipeline Overview

The workflow begins by discovering all samples in the specified FASTQ directory and processes them through the following major stages:

1.  **Read Trimming:** Low-quality bases and adapter remnants are removed from single-end reads using `Trimmomatic`.
2.  **Mapping:** Trimmed reads are aligned to the reference genome using `BWA-MEM`. Unmapped reads are filtered out.
3.  **Sorting & Indexing:** Aligned BAM files are sorted and indexed with `Samtools`.
4.  **Mark & Remove Duplicates:** PCR duplicates are identified with `Picard MarkDuplicates` and **removed** from the BAM files.
5.  **Per-Sample Calling:** `GATK HaplotypeCaller` is run on each sample individually in `GVCF` mode.
6.  **Consolidation:** All individual GVCFs are combined into a single database using `GATK CombineGVCFs`.
7.  **Joint Genotyping:** The combined GVCF is genotyped to create a single, multi-sample VCF using `GATK GenotypeGVCFs`.
8.  **Variant Separation:** The joint VCF is split into two files: one for SNPs and one for INDELs using `GATK SelectVariants`.
9.  **Hard Filtering:** `BCFtools` applies hard filters to the SNP and INDEL files to remove low-confidence variants.

### Workflow Diagram
```

( Raw Single-End FASTQ files )
│
▼
[ 1. Trimmomatic (SE) ]
│
▼
[ 2. BWA MEM → Samtools view ]
│
▼
[ 3. Samtools sort & index ]
│
▼
[ 4. Picard MarkDuplicates (REMOVE\_DUPLICATES=true) ]
│
▼
[ 5. GATK HaplotypeCaller (per-sample GVCF) ]
│
▼
[ 6. GATK CombineGVCFs ]
(Combined GVCF)
│
▼
[ 7. GATK GenotypeGVCFs ]
(Joint-called VCF)
│
├─────────────┐
▼             ▼
[ 8. GATK SelectVariants (SNPs) ]   [ 8. GATK SelectVariants (INDELs) ]
│             │
▼             ▼
[ 9. BCFtools filter ]       [ 9. BCFtools filter ]
│             │
▼             ▼
(final\_snps.vcf.gz)      (final\_indels.vcf.gz)

```

---

## Required Setup

This pipeline does **not** use Conda. It relies on a specific directory structure and manually installed software.

**1. Directory Structure:**
You must organize your project files as shown below. The `Snakefile` uses relative paths (`../../`) that depend on this layout.

```

/your\_main\_project\_directory/
├── analysis/              \<-- Your current working directory
│   ├── Snakefile          \<-- This snakemake file
|   ├── SimpleSnek.smk          \<-- This snakemake file
|   ├── environment          \<-- This environment file
│   └── ... (final output VCFs will be created here)
│
├── bin/                   \<-- Directory for executable software
│   ├── Trimmomatic-0.39/
│   │   └── trimmomatic-0.39.jar
│   ├── gatk-4.3.0.0/
│   │   └── gatk
│   └── picard.jar
│
├── reference\_genome/
│   └── final\_genome.fa    \<-- Must be indexed by BWA, Samtools, and Picard
│
└── trim/                  \<-- Pipeline will create and use these folders
└── mapped/
└── sorted/
└── marked/
└── vcf/

# The location of your input files is configured in the Snakefile

# For example: /workdir/ams866/Alfalfa/validation\_plates/FASTQ/
````

_Note: you may have to manually create the trim, mapped, sorted, marked, and vcf folders manually before running the pipeline_

**2. Software Installation:**
Ensure the following software is installed and that the paths in the `Snakefile` point to the correct locations.
* **Java**
* **Trimmomatic** (v0.39)
* **BWA**
* **Samtools**
* **Picard**
* **GATK** (v4.3.0.0)
* **BCFtools**

**3. Reference Genome Indexing:**
Before running the main pipeline, your reference genome (`final_genome.fa`) must be indexed. You can uncomment and run the `index_genome` rule in the `Snakefile` once, or run these commands manually:
```bash
bwa index /path/to/your/reference_genome/final_genome.fa
samtools faidx /path/to/your/reference_genome/final_genome.fa
java -jar /path/to/your/bin/picard.jar CreateSequenceDictionary R=/path/to/your/reference_genome/final_genome.fa O=/path/to/your/reference_genome/final_genome.dict
````

-----

## Configuration

All configuration is done by **editing the `Snakefile` directly**.

1.  **`fastq_dir`**: This is the most important variable. Change its value to the full path of the directory containing your input `.FASTQ.gz` files.
    ```python
    # FASTQ directory
    fastq_dir = "/workdir/ams866/Alfalfa/validation_plates/FASTQ"
    ```
2.  **Tool Paths and Parameters**: Check all `shell:` commands to ensure the paths to executables (`java -jar ../../bin/...`) and the reference genome (`../../reference_genome/...`) are correct for your setup. Adjust parameters like thread counts or memory (`-Xmx20g`) as needed.

-----

## Usage

Navigate to the directory containing your `Snakefile` to execute the commands below.

**1. Dry-Run (Highly Recommended):**
This will check for errors and show you which jobs will be run without actually executing them.

```bash
snakemake all -n
```

> **Note:** You must explicitly specify the `all` target. Running just `snakemake --cores <N>` without `all` may not work as expected.

**2. Full Pipeline Execution:**
Run the entire pipeline. Snakemake will use the number of cores you specify.

```bash
snakemake all --cores <N>
```

> Replace `<N>` with the total number of CPU cores you want to allocate to the job.

**3. Visualize the Workflow:**
To generate an image of the workflow graph (requires `graphviz`):

```bash
snakemake all --dag | dot -Tpng > workflow_dag.png
```
**4. Run Partial Pipeline:**
To run only up to a specific step:
```bash
snakemake map_all --cores <N>   # Stop after mapping
snakemake mark_all --cores <N>  # Stop after duplicate marking
snakemake call_all --cores <N>  # Stop after GVCF generation
```
-----

## Output Files

The primary final outputs of this pipeline will be created in the `analysis` directory (or wherever you run Snakemake from):

  * `final_combined_filtered_snps.vcf.gz`: The final set of filtered SNPs.
  * `final_combined_filtered_indels.vcf.gz`: The final set of filtered INDELs.

Intermediate files are stored in the sibling directories (`../../trim`, `../../mapped`, etc.) as defined by the directory structure.


-----

## Contributing

Contributions are welcome\! If you have suggestions for improvements or find a bug, please feel free to:

1.  Fork the repository.
2.  Create a new branch (`git checkout -b feature/your-feature`).
3.  Make your changes.
4.  Commit your changes (`git commit -am 'Add some feature'`).
5.  Push to the branch (`git push origin feature/your-feature`).
6.  Create a new Pull Request.

-----

## License

This project is licensed under the MIT License - see the [LICENSE](https://www.google.com/search?q=LICENSE) file for details.
