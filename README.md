```
--------------------------------------------------------
    888      e      Y88b    /        e88~~\  ~~~888~~~
    888     d8b      Y88b  /        d888        888   
    888    /Y88b      Y88b/         8888 __     888   
    888   /  Y88b     /Y88b         8888   |    888   
|   88P  /____Y88b   /  Y88b        Y888   |    888   
 \__8"  /      Y88b /    Y88b        "88__/     888   

                   5baseTAPS 1.0.0
--------------------------------------------------------
```

# 5baseTAPS: Nextflow Pipeline

A Nextflow pipeline for whole-genome cytosine methylation and SNP variant calling from
**Illumina 5-base** and **TAPS Watchmaker** duplex UMI sequencing data.

Built at [The Jackson Laboratory](https://www.jax.org/) Genome Technologies core, extending
[nf-core/fastquorum](https://nf-co.re/fastquorum) with rastair methylation-aware CpG methylation and SNP calling, and GATK
HaplotypeCaller variant calling.

---

## Overview

Both library types (Illumina 5-base and TAPS Watchmaker) produce reads with the same
base-level methylation signature: 5-methylcytosine (5mC) appears as thymine (C→T), while
unmethylated cytosines are read as C. The pipeline exploits this shared read representation
to run a single analysis path regardless of library chemistry.

### Workflow

**Step-by-step:**

1. Raw read QC ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. [Optional] Chunked alignment of raw FASTQs for large samples ([seqkit](https://bioinf.shenwei.me/seqkit/) + [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2))
3. UMI extraction and grouping ([fgbio FastqToBam](http://fulcrumgenomics.github.io/fgbio/), [GroupReadsByUmi](http://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html))
4. Duplex consensus calling and filtering ([fgbio CallDuplexConsensusReads](http://fulcrumgenomics.github.io/fgbio/tools/latest/CallDuplexConsensusReads.html), [FilterConsensusReads](http://fulcrumgenomics.github.io/fgbio/tools/latest/FilterConsensusReads.html))
5. Duplex QC metrics ([fgbio CollectDuplexSeqMetrics](http://fulcrumgenomics.github.io/fgbio/tools/latest/CollectDuplexSeqMetrics.html))
6. Final consensus re-alignment ([bwa-mem2](https://github.com/bwa-mem2/bwa-mem2))
7. Methylation-aware CpG methylation + SNP calling + per-CpG BED/VCF + per-read BED ([rastair](https://www.rastair.com/))
8. M-bias QC report ([rastair](https://www.rastair.com/))
9. Lambda (negative) and pUC19 (positive) methylation control QC
10. CpG site mask generation for GATK (derived from the rastair call BED)
11. DRAGstr model calibration + scatter-gather SNP/INDEL calling ([GATK HaplotypeCaller](https://gatk.broadinstitute.org/))
12. MultiQC report combining all QC metrics ([MultiQC](https://multiqc.info/))

---

## Requirements

- **Nextflow** ≥ 24.04.2
- **Singularity** (all processes run in containers; Docker or Apptainer also supported)
- **SLURM** (configured for JAX clusters; other executors require config changes)

---

## Running the Pipeline

First, clone the repository:

```bash
git clone https://github.com/TheJacksonLaboratory/5baseTAPS.git
cd 5baseTAPS
```

### Running a Quick Test

```bash
module load singularity nextflow
nextflow run . -profile sumner_test --outdir test_results/
```

### Input Samplesheet

A complete guide to running the pipeline is available in [docs/usage.md](docs/usage.md).

Prepare a CSV samplesheet with one row per FASTQ pair (multi-lane samples use multiple rows with the same `sample` name). The `read_structure` column depends on how the sequencing run was demultiplexed:

**Type 1 — Inline UMI** (UMI embedded in read sequence, raw BCL Convert output):

```csv
sample,fastq_1,fastq_2,read_structure
SAMPLE1,/path/to/SAMPLE1_R1.fastq.gz,/path/to/SAMPLE1_R2.fastq.gz,7M1S+T 7M1S+T
SAMPLE2,/path/to/SAMPLE2_L1_R1.fastq.gz,/path/to/SAMPLE2_L1_R2.fastq.gz,7M1S+T 7M1S+T
SAMPLE2,/path/to/SAMPLE2_L2_R1.fastq.gz,/path/to/SAMPLE2_L2_R2.fastq.gz,7M1S+T 7M1S+T
```

**Type 2 — UMI in read header** (BCL Convert already extracted UMI into the read name and trimmed adapters):

```csv
sample,fastq_1,fastq_2,read_structure
SAMPLE1,/path/to/SAMPLE1_R1.fastq.gz,/path/to/SAMPLE1_R2.fastq.gz,+T +T
SAMPLE2,/path/to/SAMPLE2_L1_R1.fastq.gz,/path/to/SAMPLE2_L1_R2.fastq.gz,+T +T
SAMPLE2,/path/to/SAMPLE2_L2_R1.fastq.gz,/path/to/SAMPLE2_L2_R2.fastq.gz,+T +T
```

> `7M1S+T 7M1S+T` — 7 bp UMI + 1 bp spacer trimmed from each read (Type 1)
> `+T +T` — reads are all template; UMI transferred from read name to BAM tag (Type 2)
>
> See [docs/usage.md](docs/usage.md) for full details on both input types.

Multi-lane rows with the same `sample` identifier are merged automatically before UMI grouping.

### Run via SLURM Head Script

Submit a SLURM wrapper script that calls `nextflow run` so the Nextflow process itself is
managed by the scheduler. Nextflow then submits each pipeline task as a separate SLURM job.

```bash
#!/bin/bash
#SBATCH --job-name=nf-5baseTAPS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x-%j.log
#SBATCH --error=logs/%x-%j.log

module load singularity nextflow
nextflow run . \
    -profile sumner2_singularity \
    --input samplesheet.csv \
    --genome CHM13 \
    --outdir results
```

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | — | Path to samplesheet CSV (required) |
| `--genome` | `CHM13` | Reference genome key (`CHM13`,`GRCh38`,`GRCm38`, etc.) |
| `--fasta` | — | Path to reference FASTA (if it differs from the `--genome` fasta, all indexes are rebuilt) |
| `--fasta_fai` | — | Path to FASTA index (`.fai`) |
| `--dict` | — | Path to sequence dictionary (`.dict`) |
| `--bwamem2` | — | Path to bwa-mem2 index directory |
| `--outdir` | — | Output directory (required) |
| `--run_gatk` | `true` | Run GATK HaplotypeCaller variant calling; set to false for methylation only |
| `--align_raw_bam_chunks` | `1` | Number of parallel alignment chunks applied per samplesheet row (FASTQ pair); set > 1 for high-depth FASTQ pairs (4-8) |
| `--filter_min_reads` | `1 0 0` | Minimum reads to retain a duplex consensus (format: `AB SS DS`) |
| `--filter_max_base_error_rate` | `0.1` | Maximum per-base error rate for consensus filtering |
| `--duplex_seq` | `true` | Enable duplex consensus mode (set `false` for single-strand UMI libraries) |
| `--rastair_rscript_dir` | `null` | Override rastair R scripts directory (uses bundled `bin/rastair_scripts/` by default) |

---

## Output Structure

```
<outdir>/
├── <sample>/
│   ├── bam/               — duplex consensus BAM + index
│   ├── methylation/       — per-CpG BED/VCF, per-read BED, M-bias HTML, methylKit, summaries
│   ├── variants/          — GATK VCFs, CpG site mask
│   └── qc/                — per-sample QC (alignment, coverage, duplex, FastQC, variant)
├── report/                — MultiQC HTML report
└── pipeline_info/         — Nextflow execution reports and parameter logs
```

See [docs/output.md](docs/output.md) for a complete description of all output files.

---


---

## Credits

The JAX-GT pipeline was developed by the JAX Genome Technologies bioinformatics team, built on
top of [nf-core/fastquorum](https://nf-co.re/fastquorum) (Nils Homer & Zach Norgaard,
Fulcrum Genomics) and the TAPS methylation conversion subworkflow 
adapted from [nf-core/methylseq](https://nf-co.re/methylseq) (Phil Ewels et al.;
doi:[10.5281/zenodo.1343417](https://doi.org/10.5281/zenodo.1343417)).

Full tool citations are listed in [CITATIONS.md](CITATIONS.md).

For questions or inquiries about the pipeline, contact the JAX Genome Technologies team at [GTdrylab@jax.org](mailto:GTdrylab@jax.org).
