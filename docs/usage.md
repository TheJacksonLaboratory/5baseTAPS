# 5baseTAPS: Usage

## Pipeline parameters

Provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files specified with `-c` can only be used for infrastructure configuration (resource tuning, executor settings), not pipeline parameters — see the [Nextflow configuration docs](https://www.nextflow.io/docs/latest/config.html).

## Samplesheet input

Prepare a CSV samplesheet with one row per FASTQ pair. Multi-lane samples use multiple rows with the same `sample` name. Use `--input` to specify its location:

```bash
--input '[path to samplesheet file]'
```

### Samplesheet format

5baseTAPS supports two FASTQ input types depending on how the sequencing run was demultiplexed. The `read_structure` column in the samplesheet is the key differentiator.

Multi-lane samples use multiple rows with the same `sample` name — lanes are merged automatically before UMI grouping. Single- and multi-lane samples can be combined in one file.

---

#### Input Type 1 — Inline UMI (raw FASTQ)

The UMI is embedded at the 5′ end of both R1 and R2 read sequences, immediately followed by a spacer base. This is the raw output from sequencing when BCL Convert is run **without** UMI extraction.

Read structure: `7M1S+T 7M1S+T` — 7 bp UMI + 1 bp spacer + template, on both R1 and R2.

```
Read 1: [AAACCCG][A][insert sequence ...]
          7bp UMI  sp  template (+T)
```

The pipeline trims the UMI and spacer from the read sequences using `fgbio FastqToBam`, populates the BAM `RX` tag, then proceeds to alignment and UMI grouping. Reads are uniform length (e.g. 143 bp) after trimming.

```csv
sample,fastq_1,fastq_2,read_structure
SAMPLE1,/path/to/SAMPLE1_R1.fastq.gz,/path/to/SAMPLE1_R2.fastq.gz,7M1S+T 7M1S+T
SAMPLE2,/path/to/SAMPLE2_L001_R1.fastq.gz,/path/to/SAMPLE2_L001_R2.fastq.gz,7M1S+T 7M1S+T
SAMPLE2,/path/to/SAMPLE2_L002_R1.fastq.gz,/path/to/SAMPLE2_L002_R2.fastq.gz,7M1S+T 7M1S+T
SAMPLE2,/path/to/SAMPLE2_L003_R1.fastq.gz,/path/to/SAMPLE2_L003_R2.fastq.gz,7M1S+T 7M1S+T
```

---

#### Input Type 2 — UMI in read header (BCL Convert demultiplexed)

The UMI has already been extracted by BCL Convert and placed into the FASTQ read name. Adapter trimming is also applied by BCL Convert, so reads are variable length. The read sequences contain only insert bases — no UMI prefix to trim.

Read name format (Illumina BCL Convert, dual UMI):
```
@INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y:UMI1+UMI2 1:N:0:SAMPLEINDEX
                                       ^^^^^^^^^^
                                       7+7 bp dual UMI in 8th colon field
```

Use `read_structure = "+T +T"` (all bases are template; no UMI in reads). The pipeline passes `--extract-umis-from-read-names` to `fgbio FastqToBam`, which reads the UMI from the 8th colon-delimited field of the read name and populates the BAM `RX` tag directly. No read trimming occurs. UMI grouping and consensus calling then proceed identically to Type 1.

```csv
sample,fastq_1,fastq_2,read_structure
SAMPLE1,/path/to/SAMPLE1_R1.fastq.gz,/path/to/SAMPLE1_R2.fastq.gz,+T +T
SAMPLE2,/path/to/SAMPLE2_L001_R1.fastq.gz,/path/to/SAMPLE2_L001_R2.fastq.gz,+T +T
SAMPLE2,/path/to/SAMPLE2_L002_R1.fastq.gz,/path/to/SAMPLE2_L002_R2.fastq.gz,+T +T
SAMPLE2,/path/to/SAMPLE2_L003_R1.fastq.gz,/path/to/SAMPLE2_L003_R2.fastq.gz,+T +T
```

> **Note:** All lanes of a sample must use the same `read_structure`. Different samples in the same samplesheet may use different types — mixing Type 1 and Type 2 samples in a single run is supported.

---

### Samplesheet column reference

| Column           | Description |
| ---------------- | ----------- |
| `sample`         | Sample name. Identical across rows from the same sample (multi-lane). Spaces are converted to underscores. |
| `fastq_1`        | Full path to R1 FASTQ. Extension must be `.fastq`, `.fq`, `.fastq.gz`, or `.fq.gz`. |
| `fastq_2`        | Full path to R2 FASTQ. Extension must be `.fastq`, `.fq`, `.fastq.gz`, or `.fq.gz`. |
| `read_structure` | fgbio [read structure](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures) string. Use `7M1S+T 7M1S+T` for inline UMI (Type 1) or `+T +T` for UMI-in-header (Type 2). Must be identical across all lanes of the same sample. |

### Main Options

Two modes of running this pipeline are supported:

1. **Research and Development (R&D):** `--mode rd`. Allows branching off from the pipeline to test multiple consensus calling or filtering parameters.
2. **High Throughput (HT):** `--mode ht`. Intended for production environments where performance takes precedence.

For duplex sequencing (reads from the same source molecule may observe either strand), set `--duplex_seq true`. This uses [`fgbio CallDuplexConsensusReads`](https://fulcrumgenomics.github.io/fgbio/tools/latest/CallDuplexConsensusReads.html). Otherwise, [`fgbio CallMolecularConsensusReads`](https://fulcrumgenomics.github.io/fgbio/tools/latest/CallMolecularConsensusReads.html) is used.

### Consensus Calling Options

Options for [`fgbio CallMolecularConsensusReads`](https://fulcrumgenomics.github.io/fgbio/tools/latest/CallMolecularConsensusReads.html) and [`CallDuplexConsensusReads`](https://fulcrumgenomics.github.io/fgbio/tools/latest/CallDuplexConsensusReads.html) are prefixed with `call_`.

- `--call_min_reads` — minimum read count to call a consensus
- `--call_min_baseq` — minimum input base quality to use when calling a consensus

### Consensus Filtering Options

Options for [`fgbio FilterConsensusReads`](https://fulcrumgenomics.github.io/fgbio/tools/latest/FilterConsensusReads.html) are prefixed with `filter_`.

- `--filter_min_reads` — minimum read count to retain a consensus (up to three values for duplex reads, e.g. `'1 0 0'`)
- `--filter_min_baseq` — minimum base quality to retain a consensus base
- `--filter_max_base_error_rate` — maximum error rate for a single consensus base

### Reference Genome Options

**Pre-built genomes** (use `--genome <key>`):

#### JAX-GT pre-built (recommended)

| Key | Species | Assembly | Source |
|-----|---------|----------|--------|
| `CHM13` | *Homo sapiens* | T2T-CHM13v2 | UCSC **(default)** |
| `GRCh38` | *Homo sapiens* | GRCh38 | NCBI |

These genomes have BWAmem2 indices pre-built and validated on the JAX HPC. Other keys below follow the standard iGenomes path convention and require BWAmem2 indices to be built locally.

#### Supported iGenomes

**Homo sapiens**

| Key | Assembly | Source |
|-----|----------|--------|
| `CHM13`  | T2T-CHM13v2 | UCSC |
| `GRCh38` | GRCh38 | NCBI |
| `GRCh37` | GRCh37 | Ensembl |
| `hg38` | hg38 | UCSC |
| `hg19` | hg19 | UCSC |

**Mus musculus**

| Key | Assembly | Source |
|-----|----------|--------|
| `GRCm38` | GRCm38 | Ensembl |
| `mm10` | mm10 | UCSC |

**Other organisms**

| Key | Species | Assembly | Source |
|-----|---------|----------|--------|
| `bosTau8` | *Bos taurus* | bosTau8 | UCSC |
| `UMD3.1` | *Bos taurus* | UMD3.1 | Ensembl |
| `CanFam3.1` | *Canis familiaris* | CanFam3.1 | Ensembl |
| `canFam3` | *Canis familiaris* | canFam3 | UCSC |
| `WBcel235` | *C. elegans* | WBcel235 | Ensembl |
| `ce10` | *C. elegans* | ce10 | UCSC |
| `GRCz10` | *Danio rerio* | GRCz10 | Ensembl |
| `danRer10` | *Danio rerio* | danRer10 | UCSC |
| `BDGP6` | *Drosophila melanogaster* | BDGP6 | Ensembl |
| `dm6` | *Drosophila melanogaster* | dm6 | UCSC |
| `EquCab2` | *Equus caballus* | EquCab2 | Ensembl |
| `equCab2` | *Equus caballus* | equCab2 | UCSC |
| `Galgal4` | *Gallus gallus* | Galgal4 | Ensembl |
| `galGal4` | *Gallus gallus* | galGal4 | UCSC |
| `Mmul_1` | *Macaca mulatta* | Mmul_1 | Ensembl |
| `CHIMP2.1.4` | *Pan troglodytes* | CHIMP2.1.4 | Ensembl |
| `panTro4` | *Pan troglodytes* | panTro4 | UCSC |
| `Rnor_6.0` | *Rattus norvegicus* | Rnor_6.0 | Ensembl |
| `rn6` | *Rattus norvegicus* | rn6 | UCSC |
| `R64-1-1` | *Saccharomyces cerevisiae* | R64-1-1 | Ensembl |
| `sacCer3` | *Saccharomyces cerevisiae* | sacCer3 | UCSC |
| `Sscrofa10.2` | *Sus scrofa* | Sscrofa10.2 | Ensembl |
| `susScr3` | *Sus scrofa* | susScr3 | UCSC |
| `TAIR10` | *Arabidopsis thaliana* | TAIR10 | Ensembl |
| `AGPv3` | *Zea mays* | AGPv3 | Ensembl |
| `Sbi1` | *Sorghum bicolor* | Sbi1 | Ensembl |
| `Gm01` | *Glycine max* | Gm01 | Ensembl |
| `IRGSP-1.0` | *Oryza sativa japonica* | IRGSP-1.0 | Ensembl |
| `EB2` | *Bacillus subtilis 168* | EB2 | Ensembl |
| `EB1` | *E. coli K-12 DH10B* | EB1 | Ensembl |
| `EF2` | *Schizosaccharomyces pombe* | EF2 | Ensembl |

#### Explicit reference file specification

You can provide reference files directly instead of using `--genome`:

- `--fasta` — path to the genome FASTA file
- `--fasta_fai` — path to the FASTA index (`samtools faidx`)
- `--dict` — path to the sequence dictionary (`samtools dict`)
- `--bwamem2` — path to the **directory** containing the bwa-mem2 index

Use `--save_reference` to save generated indices in the results directory for reuse in subsequent runs.

---

## Running the pipeline

### On an HPC with SLURM (recommended)

The recommended way to run the pipeline on an HPC cluster is to submit a SLURM job that runs Nextflow as the workflow manager. Nextflow then submits individual process jobs to the cluster scheduler automatically.

Example SLURM submission script:

```bash
#!/usr/bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -p compute
#SBATCH --mem=24G
#SBATCH -t 2-15:55
#SBATCH --job-name=run-5base
#SBATCH --output="%x.%j.slurm.o"
#SBATCH --error="%x.%j.slurm.e"

module load singularity
module load nextflow/25.04.2

mkdir -p /flashscratch/$USER/5base-run && cd /flashscratch/$USER/5base-run

genRef=GRCh38
FASTA=/flashscratch/nf-JAX-5base/genomes/Homo_sapiens/NCBI/$genRef/Sequence/WholeGenomeFasta/genome.fa
BWAMEM2=/flashscratch/nf-JAX-5base/genomes/Homo_sapiens/NCBI/$genRef/Sequence/BWAmem2Index
csv=/path/to/samplesheet.csv

nextflow run /path/to/5baseTAPS/main.nf \
    --input $csv \
    --genome $genRef \
    --bwamem2   ${BWAMEM2} \
    --fasta     ${FASTA} \
    --fasta_fai ${FASTA}.fai \
    --dict      ${FASTA/fa/dict} \
    -profile sumner2_apptainer \
    --filter_min_reads '1 0 0' \
    --outdir result_$genRef
```

> **Note:** Replace `sumner2_apptainer` with the profile matching your HPC site (see `conf/` for available profiles). The `module load` commands and paths are site-specific.
>
> **Note:** If you are using Nextflow **v26 or later**, add `export NXF_ANSI_LOG=false` before the `nextflow run` command.

### Parallel alignment chunks (`--align_raw_bam_chunks`)

For large samples (≥500M reads), alignment is the pipeline bottleneck. The `--align_raw_bam_chunks INT` parameter splits the raw BAM into N chunks that are aligned in parallel then merged, reducing wall time:

```bash
nextflow run /path/to/5baseTAPS/main.nf \
    ...
    --align_raw_bam_chunks 4
```

Default is `1` (no splitting). A value of 4 gives approximately 1.8× wall-time speedup on 500M+ read samples.

> **Warning:** Each chunk becomes a separate SLURM job submitted concurrently. Most HPC schedulers enforce a per-user job submission limit. Running multiple samples with a high chunk count can exhaust this limit, causing jobs to queue rather than run in parallel and negating the speedup. As a practical guide: `(number of samples) × (chunk count)` alignment jobs will be submitted simultaneously — keep this product within your site's job limit. When in doubt, use `--align_raw_bam_chunks 4` as a conservative default.


### Annotated VCFs (`--annotate_vcf`)

Optional and off by default. With `--annotate_vcf true` the pipeline writes **additional, annotated copies** of the variant files. The original VCFs are not modified.

```bash
nextflow run /path/to/5baseTAPS/main.nf ... --genome GRCh38 --annotate_vcf true
```

| Input VCF | New annotated file | Location |
|-----------|--------------------|----------|
| `<sample>.haplotypecaller.filtered-SNP.vcf.gz` | `<sample>.haplotypecaller.filtered-SNP.annotated.vcf.gz` (+ `.tbi`) | `<sample>/variants/` |
| `<sample>.haplotypecaller.filtered-INDEL.vcf.gz` | `<sample>.haplotypecaller.filtered-INDEL.annotated.vcf.gz` (+ `.tbi`) | `<sample>/variants/` |
| `<sample>.rastair_call.vcf.gz` | `<sample>.rastair_call.annotated.vcf.gz` (+ `.tbi`) | `<sample>/methylation/` |

**What is a subset of the original**

* The GATK annotated files contain the same records as the corresponding original VCFs.
* **The rastair annotated file is a subset of `rastair_call.vcf.gz`.** The original is mostly CpG methylation rows with no variant (`ALT=.`); only records with a real ALT allele are annotated and kept, and only on the chromosomes in `--annotate_chroms` (default `chr1`–`chr22`, `chrX`, `chrY`). Variant records on other contigs (alt, random, unplaced) are not in the annotated file. Add contigs with `--annotate_chroms` if you need them.

**Requirements**

* `--genome GRCh38` only. Any other genome stops the run with an error, because the annotation resources are GRCh38.
* The resources are read from the genome configuration (the `GRCh38` block of `conf/genomes_elion.config` on the JAX Elion cluster): `dbsnp`, `gnomad_af`, `clinvar`, `dbnsfp`, `snpeff_db` and `snpeff_config`. A site that is not Elion must define these six entries for GRCh38 in its own genome config; if any is missing the run stops with an error naming it.
* `dbsnp`, `gnomad_af` and `clinvar` must be bgzipped VCFs with a tabix index (`.tbi` next to the file) and `chr`-prefixed contig names (`chr1`, `chr2`, ...), like the pipeline's VCFs. **A resource with different contig names (for example `1`, `2`) does not cause an error: the step finishes successfully and annotates nothing.** The NCBI ClinVar download uses bare names, see [Preparing the ClinVar file](#preparing-the-clinvar-file).
* `dbnsfp` is the opposite: it must be a bare-chromosome-name file (no `chr` prefix), because SnpSift `dbnsfp` strips `chr` from the VCF contigs before looking up the database; a `chr`-prefixed database annotates nothing. The pipeline stops if the dbNSFP step annotates 0 entries although SnpEff called 20 or more missense variants on that chromosome.
* `snpeff_db` is the SnpEff database name (`hg38` on Elion) and `snpeff_config` the path to the `snpEff.config` whose `data.dir` contains it.

Example entries for a site that defines its own genome config (the paths are placeholders; the Elion values are in `conf/genomes_elion.config`):

```groovy
params {
    genomes {
        GRCh38 {
            // ... fasta, fasta_fai, dict, bwamem2 ...
            dbsnp         = "/path/to/dbsnp_151.vcf.gz"                        // + .tbi, chr-prefixed
            gnomad_af     = "/path/to/af-only-gnomad.hg38.vcf.gz"              // + .tbi, chr-prefixed
            clinvar       = "/path/to/clinvar_<date>.GRCh38.chr.vcf.gz"        // + .tbi, chr-prefixed
            dbnsfp        = "/path/to/dbNSFP4.2a.gatk_formatted.nochr.txt.gz"  // + .tbi + .data_types, no 'chr'
            snpeff_db     = "hg38"
            snpeff_config = "/path/to/snpEff.config"
        }
    }
}
```

**Cost:** the annotation is scattered per chromosome, so it adds about 220 tasks per sample (24 chromosomes × 3 VCFs × 3 steps, plus 4 selection and merge tasks). Keep your scheduler's job-submission limit in mind, as for `--align_raw_bam_chunks`.

The fields added are described in [docs/output.md](output.md#annotated-vcfs).

#### Preparing the ClinVar file

The NCBI file (`https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz`) names contigs `1`, `2`, ... `MT` and includes a few records on non-primary contigs. The recipe below renames the primary contigs to the `chr` names of your reference, drops the other records and indexes the result. It needs `bcftools` and your reference `.fai`:

```bash
FAI=/path/to/genome.fa.fai          # GRCh38 reference with chr-prefixed names (chrM for the mitochondrion)
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.md5
[ "$(cut -d' ' -f1 clinvar.vcf.gz.md5)" = "$(md5sum clinvar.vcf.gz | cut -d' ' -f1)" ] || echo "md5 mismatch"
DATE=$({ zcat clinvar.vcf.gz || true; } | grep -m1 '^##fileDate=' | cut -d= -f2 | tr -d -)

zcat clinvar.vcf.gz | awk -v FAI=$FAI 'BEGIN{FS=OFS="\t"; while((getline l < FAI)>0){split(l,a,"\t"); len[a[1]]=a[2]}}
/^##/ {print; next}
/^#CHROM/ { n=split("chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM",c," "); for(i=1;i<=n;i++) print "##contig=<ID=" c[i] ",length=" len[c[i]] ">"; print; next}
{ k=$1; if(k=="MT") k="M"; if(k ~ /^([1-9]|1[0-9]|2[0-2]|X|Y|M)$/){ $1="chr" k; print } }' \
  | bcftools view -Oz -o clinvar_${DATE}.GRCh38.chr.vcf.gz
bcftools index -t clinvar_${DATE}.GRCh38.chr.vcf.gz
```

Name the output with the ClinVar release date, as above, so a later refresh becomes a new file and the version used is visible in the VCF header's command lines.

### Output files

The pipeline creates the following in your working directory:

```bash
work/           # Nextflow working directory (intermediate files)
<OUTDIR>/       # Final results (defined with --outdir)
.nextflow_log   # Nextflow log file
```

For a full description of all output files and directory structure, see [docs/output.md](output.md).

### Using a params file

A params file is a plain YAML file you create yourself that lists pipeline parameters, replacing the need to type them on the command line each run. Nextflow does not generate this file automatically (the auto-generated `pipeline_info/params_<timestamp>.json` after each run records what was used, but is not directly reusable as `-params-file` input).

Create a file (e.g. `params.yaml`) with your parameters:

```yaml
input: '/path/to/samplesheet.csv'
outdir: './results/'
genome: 'CHM13'
filter_min_reads: '1 0 0'
```

Then pass it with `-params-file`:

```bash
nextflow run /path/to/5baseTAPS/main.nf -profile sumner2_apptainer -params-file params.yaml
```

### Reproducibility

Pin a specific pipeline version with `-r`:

```bash
nextflow run /path/to/5baseTAPS/main.nf -r 1.0.0 ...
```

The version is recorded in all pipeline reports and the MultiQC summary.

---

## Core Nextflow arguments

> **Note:** Core Nextflow options use a _single_ hyphen; pipeline parameters use a _double_ hyphen.

### `-profile`

Selects a configuration profile for your compute environment. Available profiles:

- `singularity` — use Singularity containers (recommended on HPC)
- `docker` — use Docker containers (for local development)
- `apptainer` — use Apptainer containers
- `conda` — use Conda environments (last resort; containers preferred)
- Site-specific profiles (e.g. `sumner2_apptainer`) — see `conf/` directory for all available profiles

Multiple profiles can be combined: `-profile singularity,custom_config`.

### `-resume`

Resume a previous run using cached results where inputs are unchanged:

```bash
nextflow run /path/to/5baseTAPS/main.nf -resume ...
```

Use `nextflow log` to list previous run names and resume a specific one with `-resume [run-name]`.

### `-c`

Specify a custom config file for infrastructure tuning (resource adjustments, executor settings). Do not use `-c` for pipeline parameters — use `--param` flags or `-params-file` instead.

---

## Custom configuration

### Resource requests

Each pipeline step has default CPU/memory/time allocations defined in `conf/base.config`. Steps that fail with retryable errors are automatically resubmitted with increased resources (2×, then 3×). To override defaults, see the [Nextflow configuration docs](https://www.nextflow.io/docs/latest/config.html).

### Custom containers

To override the container used by a specific process, add a `withName:` selector to a custom config file:

```groovy
process {
    withName: 'PROCESS_NAME' {
        container = 'your-custom-image:tag'
    }
}
```

### Custom tool arguments

Pass additional arguments to individual tools via `ext.args` in a custom config:

```groovy
process {
    withName: 'PROCESS_NAME' {
        ext.args = '--your-extra-flag value'
    }
}
```

---

## Nextflow memory requirements

Limit Nextflow's JVM memory usage by adding to your `~/.bashrc`:

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
