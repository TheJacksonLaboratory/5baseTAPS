# 5baseTAPS: Output

## 5-base/TAPS Overview

JAX Genome Technologies offers sequencing and secondary analysis services for whole-genome cytosine methylation studies using two library types:

* **[Illumina 5-base](https://www.illumina.com/science/genomics-research/articles/5-base-solution.html)** — Illumina's commercial duplex UMI library preparation using proprietary enzymatic methylation conversion chemistry
* **[TAPS Watchmaker](https://www.watchmakergenomics.com/taps.html?tab=OVERVIEW)** — Watchmaker Genomics' library preparation using TAPS (TET-assisted pyridine borane sequencing) chemistry

Although the two library types use distinct chemistries, both produce reads with the same base-level methylation signature: 5-methylcytosine (5mC) appears as thymine (C→T) while unmethylated cytosines are read as C. This shared read-level representation is why the same downstream bioinformatics pipeline can process data from either library type without modification.

The JAX-GT Nextflow pipeline integrates frameworks from [nf-core/fastquorum](https://nf-co.re/fastquorum) and [fgbio](http://fulcrumgenomics.github.io/fgbio/) for UMI-based duplex consensus calling, [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) for genome alignment, [rastair](https://www.rastair.com/) for methylation-aware CpG methylation and SNP calling, and [GATK HaplotypeCaller](https://gatk.broadinstitute.org) in DRAGEN mode for germline SNP/INDEL calling.

**Main Steps:**
1. Preprocessing: raw read quality (FastQC), alignment to reference genome (bwa-mem2), multi-lane merging.
2. UMI consensus: UMI grouping (fgbio GroupReadsByUmi), duplex/molecular consensus calling and filtering (fgbio), consensus re-alignment (bwa-mem2).
3. Methylation calling: CpG methylation calls (rastair), M-bias diagnostics, methylKit export, methylation summaries with spike-in controls.
4. Variant calling: GATK HaplotypeCaller (DRAGEN STR model, scattered by chromosome), DRAGENHardQUAL filtering, CpG-site masking to remove TAPS 5mC→T conversion artifacts.
5. QC: coverage (mosdepth), alignment stats (samtools flagstat), integrated MultiQC report (read the [interpretation guide](./interpret_multiqc.md)).

**Supported reference genomes** (use `--genome <key>`):

Pre-built (bwa-mem2 indices locally cached on JAX HPC):
* `CHM13` — Homo sapiens T2T CHM13v2 (default)
* `GRCh38` — Homo sapiens GRCh38 / hg38 [NCBI]

Available via [iGenomes](https://nf-co.re/docs/usage/reference_genomes) (downloaded automatically on first use):
* `GRCh37` — Homo sapiens GRCh37 [Ensembl]
* `hg38` — Homo sapiens hg38 [UCSC]
* `hg19` — Homo sapiens hg19 [UCSC]
* `GRCm38` — Mus musculus GRCm38 [Ensembl]
* `mm10` — Mus musculus mm10 [UCSC]
* and more — see [full iGenomes list](https://nf-co.re/docs/usage/reference_genomes)

## Pipeline Output Structure

```
$OUTDIR/
    ├── report/
    │   ├── 5-baseTAPS_multiqc_report.html                 ← Start here
    │   ├── 5-baseTAPS_multiqc_report_data/
    │   └── 5-baseTAPS_multiqc_report_plots/
    ├── <sample>/
    │   ├── bam/
    │   │   ├── <sample>.cons.filtered.bam                 ← Primary BAM
    │   │   └── <sample>.cons.filtered.bam.bai
    │   ├── methylation/
    │   │   ├── <sample>.rastair_call.bed.gz               ← Primary methylation
    │   │   ├── <sample>.rastair_call.vcf.gz
    │   │   ├── <sample>.rastair_call.annotated.vcf.gz         ← only with --annotate_vcf (ALT-only subset)
    │   │   ├── <sample>.rastair_call.annotated.vcf.gz.tbi
    │   │   ├── <sample>.rastair_perread.bed.gz
    │   │   ├── <sample>.rastair_perread.bed.gz.tbi
    │   │   ├── <sample>.rastair_methylkit.txt.gz
    │   │   ├── <sample>.rastair_mbias.html
    │   │   ├── <sample>.methylation_summary.tsv
    │   │   ├── <sample>.lambda_negCtrl.methylation_summary.tsv
    │   │   └── <sample>.puc19_posCtrl.methylation_summary.tsv
    │   ├── variants/
    │   │   ├── <sample>.haplotypecaller.filtered-SNP.vcf.gz    ← Primary SNPs
    │   │   ├── <sample>.haplotypecaller.filtered-SNP.vcf.gz.tbi
    │   │   ├── <sample>.haplotypecaller.filtered-INDEL.vcf.gz
    │   │   ├── <sample>.haplotypecaller.filtered-INDEL.vcf.gz.tbi
    │   │   ├── <sample>.haplotypecaller.filtered-{SNP,INDEL}.annotated.vcf.gz{,.tbi}   ← only with --annotate_vcf
    │   │   ├── <sample>.haplotypecaller.vcf.gz
    │   │   ├── <sample>.haplotypecaller.vcf.gz.tbi
    │   │   ├── <sample>.rastair_cpg_sites.bed.gz
    │   │   └── <sample>.rastair_cpg_sites.bed.gz.tbi
    │   └── qc/
    │       ├── fastqc/
    │       │   ├── <sample>_1_fastqc.html
    │       │   ├── <sample>_1_fastqc.zip
    │       │   ├── <sample>_2_fastqc.html
    │       │   └── <sample>_2_fastqc.zip
    │       ├── alignment/
    │       │   ├── <sample>.prededup.flagstat
    │       │   ├── <sample>.postdedup.flagstat
    │       │   ├── <sample>.lambda_negCtrl.flagstat
    │       │   ├── <sample>.puc19_posCtrl.flagstat
    │       │   └── <sample>.cons.stats
    │       ├── coverage/
    │       │   ├── <sample>.mosdepth.summary.txt
    │       │   └── <sample>.mosdepth.global.dist.txt
    │       ├── duplex/
    │       │   ├── <sample>.grouped-family-sizes.txt
    │       │   ├── <sample>.duplex_seq_metrics.duplex_yield_metrics.txt
    │       │   ├── <sample>.duplex_seq_metrics.family_sizes.txt
    │       │   ├── <sample>.duplex_seq_metrics.duplex_family_sizes.txt
    │       │   ├── <sample>.duplex_seq_metrics.umi_counts.txt
    │       │   ├── <sample>.duplex_seq_metrics.duplex_umi_counts.txt
    │       │   └── <sample>.duplex_seq_metrics.duplex_qc.pdf
    │       ├── variant_call/
    │       │   ├── <sample>.bcftools_stats.txt
    │       │   └── <sample>.dragstr_model.txt
    │       └── mqc/
    │           ├── <sample>_read_partitioning_mqc.tsv
    │           ├── <sample>_mapping_summary_mqc.tsv
    │           ├── <sample>_duplex_summary_mqc.tsv
    │           ├── <sample>_duplex_familysize_mqc.tsv
    │           ├── <sample>_methyl_summary_mqc.tsv
    │           ├── <sample>_methyl_controls_mqc.tsv
    │           ├── <sample>_coverage_metrics_mqc.tsv
    │           └── <sample>_variant_calling_mqc.tsv
    └── pipeline_info/
        ├── execution_trace_*.txt
        └── params_*.json
```

---

## 📂 Primary Data Outputs

The following files are the "source of truth" for downstream biological analysis.

### 1. Consensus BAM (final)

The UMI consensus pipeline (`consensus_alignment/`) processes reads through four steps; only the final BAM is published:

| Step | Output | Published |
|------|--------|-----------|
| 1. GroupReadsByUmi | `grouped.bam` | No — intermediate |
| 2. CallDuplexConsensusReads | `cons.unmapped.bam` | No — intermediate |
| 3. AlignConsensusBAM | `mapped.bam` | No — intermediate |
| 4. FilterConsensusReads | **`cons.filtered.bam`** | **Yes — final BAM** |

* **`<sample>/bam/<sample>.cons.filtered.bam`**: The final consensus BAM. Reads are aligned to the reference genome with bwa-mem2 and quality-filtered by fgbio FilterConsensusReads. This is the input to rastair (methylation calling), mosdepth (coverage), GATK HaplotypeCaller (variant calling), and samtools flagstat (post-dedup QC). The `grouped-family-sizes.txt` from step 1 is still published separately for duplex QC metrics.

---
#### Example 1 — Single-strand consensus (aD=1, bD=0)
A molecule where only A-strand reads were captured (most common at typical sequencing depths).
```
GT26-11111:17	163	chr1	10000	17	19M1D13M4I36M3I63M1D5M	=	10053	180	CTAACCCTAACCCTAACCC...	EEEEE...HHHII	MC:Z:32M1D61M1D32M18S	MD:Z:0N18^T5T5T100^T5	RG:Z:A	MI:Z:17	NM:i:12	MQ:i:17	AS:i:94	XS:i:92	RX:Z:AACTAAC-GTTGTAT	aD:i:1	bD:i:0	cD:i:1	aE:f:0	bE:f:0	cE:f:0	aM:i:1	bM:i:0	cM:i:1	ac:Z:CTAACC...	ad:B:s,1,1,...	ae:B:s,0,0,...	aq:Z:EEEEE...	ms:i:101
```

#### Example 2 — True duplex consensus (aD=1, bD=1, cD=2)
A molecule with both A-strand and B-strand reads. Note `bc`/`bd`/`be`/`bq` tags are present alongside `ac`/`ad`/`ae`/`aq`.

**Position 15 — discordant CpG:**
The reference base at position 15 is **C** (part of a CpG dinucleotide, confirmed by `MD:Z:14c0g`):
- A-strand (`ac`): **T** — 5mC→T TAPS conversion, indicating methylation
- B-strand (`bc`): **C** — cytosine retained, indicating no methylation (or incomplete conversion)
- Consensus SEQ: **N** — fgbio masks positions where A-strand and B-strand disagree, assigning quality `#` (Phred 2)
```
GT26-1111:22743	99	chr1	10768243	60	143M	=	10768360	260	TAGTCCCAGCTACTNGGAGGCTGAGGTGGGAGG...	EEEEEEE8EEEEEEE8...qqqqq	MC:Z:143M	MD:Z:14c0g20t106	RG:Z:A	MI:Z:22743	NM:i:3	MQ:i:60	AS:i:131	XS:i:63	RX:Z:NCGTTGT-TACTCAT	aD:i:1	bD:i:1	cD:i:2	aE:f:0	bE:f:0	cE:f:0.00699301	aM:i:1	bM:i:1	cM:i:2	ac:Z:TAGTCCCAGCTACTТGGG...	bc:Z:TAGTCCCAGCTACTCAGG...	ad:B:s,1,1,...	bd:B:s,1,1,...	ae:B:s,0,0,...	be:B:s,0,0,...	aq:Z:EEEEEEE8...	bq:Z:E8EEEEEE...	ms:i:138
```

##### Read Name
- fgbio replaces the Illumina instrument read name with **`{sampleName}:{moleculeID}`** — e.g. `GT26-1111:17`
- `GT26-1111` is the sample ID; `17` is the UMI molecule number assigned by GroupReadsByUmi

---

##### Core Alignment Fields

| Field | Value | Meaning |
|-------|-------|---------|
| FLAG | 99 / 163 | Properly paired; first or second in pair |
| RNAME | chr1 | Chromosome |
| POS | 10000 | 1-based alignment start |
| MAPQ | 17 | Mapping quality |
| CIGAR | 19M1D13M4I36M3I63M1D5M | Alignment with indels |
| RNEXT | = | Mate on same chromosome |
| PNEXT | 10053 | Mate position |
| TLEN | 180 | Template length |

---

##### Sequence & Quality

- `SEQ` → Consensus read sequence (methylated CpG sites appear as T due to TAPS conversion; unmethylated C retained)
- `QUAL` → Phred-scaled base qualities from fgbio consensus error model

---

##### Standard Alignment Tags

| Tag | Meaning |
|-----|---------|
| RG:Z | Read group (`A` = duplex strand A) |
| NM:i | Edit distance to reference |
| MD:Z | Mismatch/deletion string |
| AS:i | Alignment score |
| XS:i | Suboptimal alignment score |
| UQ:i | Phred-scaled probability of the alignment being incorrect |
| MC:Z | Mate CIGAR |
| MQ:i | Mate mapping quality |

---

##### UMI and fgbio Consensus Tags ([fgbio CallDuplexConsensusReads documentation](https://fulcrumgenomics.github.io/fgbio/tools/latest/CallDuplexConsensusReads.html))

```bash
RX:Z:AACTAAC-GTTGTAT   # Duplex UMI: strand A - strand B (hyphen separator)
MI:Z:17                 # Molecule ID assigned by GroupReadsByUmi

# Max raw-read depth at any position in the consensus (a=A-strand, b=B-strand, c=total)
aD:i:1   bD:i:0   cD:i:1

# Min raw-read depth at any position in the consensus (0 = some positions had no coverage)
aM:i:1   bM:i:0   cM:i:1

# Fraction of raw-read bases disagreeing with the final consensus call (error rate)
aE:f:0   bE:f:0   cE:f:0

# Per-base arrays — A-strand always present; B-strand tags (bc/bd/be/bq) only present if bD>0
ac:Z:...        # A-strand single-strand consensus bases
bc:Z:...        # B-strand single-strand consensus bases (duplex only)
ad:B:s,...      # A-strand per-base raw read depth (int16 array)
bd:B:s,...      # B-strand per-base raw read depth (duplex only)
ae:B:s,...      # A-strand per-base disagreeing base count (int16 array)
be:B:s,...      # B-strand per-base disagreeing base count (duplex only)
aq:Z:...        # A-strand per-base consensus quality
bq:Z:...        # B-strand per-base consensus quality (duplex only)
ms:i:138        # Mate score
```

> **Note:** Unlike bismark-aligned BAMs used in DRAGEN 5-base, this BAM does **not** carry `XM:Z` methylation string tags. Methylation calls are produced separately by rastair and delivered as BED/VCF files (see Section 2 below).

---

### 2. CpG Methylation Calls

* **`<sample>/methylation/<sample>.rastair_call.bed.gz`**: Per-CpG methylation calls in BED format. Each record is a CpG site covered by at least one duplex consensus read. Only sites with actual read coverage are reported — zero-coverage reference CpG sites are excluded.

The following is an example rastair call BED record (header shown for reference):
```
#chr      start     end       name  beta_est  strand  unmod  mod  no_snp  snp  coverage  genotype  gt_p_score  gt_conf_score  cpg
chr1      12986558  12986559  .     1.00      +       0      3    3       0    6         C/C       99          9              REF
```

#### Column Description

| Column | Name | Description |
|--------|------|-------------|
| 1 | chr | Reference chromosome |
| 2 | start | 0-based genomic start of CpG |
| 3 | end | 1-based genomic end of CpG |
| 4 | name | Always `.` |
| 5 | beta_est | Methylation fraction 0.00–1.00 (mod / (mod + unmod)) |
| 6 | strand | `+` (C on forward strand) or `-` (G on reverse strand) |
| 7 | unmod | Reads with C at this position (unmodified, 5mC not converted) |
| 8 | mod | Reads with T at this position (TAPS: 5mC→T conversion) |
| 9 | no_snp | Reads supporting the reference allele (no SNP) |
| 10 | snp | Reads supporting a SNP allele |
| 11 | coverage | Total reads at this site (unmod + mod + snp) |
| 12 | genotype | Called genotype (e.g., `C/C`, `G/G`) |
| 13 | gt_p_score | Genotype posterior probability score |
| 14 | gt_conf_score | Genotype confidence score |
| 15 | cpg | CpG status: `REF` (reference CpG) or `NEW` (SNP creates a new CpG context not present in the reference) |

---

* **`<sample>/methylation/<sample>.rastair_call.vcf.gz`**: Combined methylation + SNP calls from rastair. Contains two record types:
  - **CpG records** (INFO flag `CPG`, `ALT=.`): per-CpG methylation calls, same data as the BED file. Methylation evidence in INFO field `M5mC_Strands` (4 integers: unmod, mod, no_snp, snp reads) and FORMAT field `M5mC` (methylation level).
  - **SNP records** (no `CPG` flag, `ALT` is an actual allele): genome-wide SNP calls detected simultaneously from the same duplex reads. FORMAT field is `ML` rather than `M5mC`.

  This is a rastair-specific VCF format, not a standard GATK/DeepVariant variant VCF.

  > **Note:** The SNP calls here are a by-product of rastair's methylation calling — they are produced opportunistically from the same duplex reads without the dedicated SNP-calling optimisations (DRAGEN model, DragSTR calibration, scatter-gather per chromosome) applied by the GATK HaplotypeCaller workflow. For variant analysis, prefer the `variants/<sample>.haplotypecaller.*.vcf.gz` outputs.

---

### 3. Per-Read Methylation

* **`<sample>/methylation/<sample>.rastair_perread.bed.gz`** / **`.tbi`**: Per-read CpG methylation assignments used for single-molecule methylation analysis or epiallele phasing. The tabix index enables fast random-access queries by genomic region.

---

### 4. Methylation Summary Statistics

* **`<sample>/methylation/<sample>.methylation_summary.tsv`**: Per-chromosome and genome-wide CpG coverage and methylation summary. The `ALL` row gives genome-wide totals and is the source for the MultiQC methylation metrics section.

| Column | Description |
|--------|-------------|
| sample | Sample ID |
| chr | Chromosome (`ALL` = genome-wide totals) |
| global_meth_pct | Global CpG methylation % |
| cpg_cov_ge1x | CpG sites with ≥ 1× duplex coverage |
| cpg_cov_ge5x | CpG sites with ≥ 5× duplex coverage |
| cpg_cov_ge10x | CpG sites with ≥ 10× duplex coverage |
| cpg_cov_ge20x | CpG sites with ≥ 20× duplex coverage |

---

### 5. Spike-in Control Methylation Summaries

* **`<sample>/methylation/<sample>.lambda_negCtrl.methylation_summary.tsv`**: Methylation summary for the lambda phage (unmethylated negative control) spike-in.
* **`<sample>/methylation/<sample>.puc19_posCtrl.methylation_summary.tsv`**: Methylation summary for the pUC19 plasmid (fully methylated positive control) spike-in.

These files are used to estimate TAPS conversion efficiency and validate assay performance:

| Control | Expected | Interpretation |
|---------|----------|----------------|
| Lambda (negCtrl) | < 2% methylated | Conversion efficiency = 100 − lambda % |
| pUC19 (posCtrl) | > 95% methylated | Confirms 5mC detection is working |

---

### 6. M-bias Report

* **`<sample>/methylation/<sample>.rastair_mbias.html`**: Self-contained HTML report with four diagnostic sections:

| Section | Description | Controlled by |
|---------|-------------|---------------|
| **M-bias** | Per-position methylation level (β) for each read pair and strand across read length. Flat = clean; elevated at ends = trimming artifact or incomplete conversion. | Always produced |
| **V-bias** | Per-fragment heatmap of β vs. position-in-fragment and insert size, split by orientation. Detects position/orientation-specific conversion bias not visible in the marginal M-bias curves. | `--plot_vbias` (default: on) |
| **GC bias** | Methylation level vs. GC content of the local genomic window around each CpG site. | `--plot_gc` (default: on) |
| **CpG density bias** | Methylation level vs. local CpG density. | `--plot_gc` (default: on) |

---

### 7. methylKit Output

* **`<sample>/methylation/<sample>.rastair_methylkit.txt.gz`**: CpG methylation in [methylKit](https://bioconductor.org/packages/release/bioc/html/methylKit.html) format for R-based differential methylation analysis.

| Column | Description |
|--------|-------------|
| chrBase | chr.position composite key |
| chr | Chromosome |
| base | 1-based cytosine position |
| strand | `F` (forward) or `R` (reverse) |
| coverage | Total duplex read depth |
| freqC | % methylated (fraction of reads with T at this TAPS position) |
| freqT | % unmethylated |

---

### 8. Variant Calls

* **`<sample>/variants/<sample>.haplotypecaller.filtered-SNP.vcf.gz`** / **`.tbi`**: PASS germline SNPs after GATK DRAGENHardQUAL filtering and CpG-site masking. The CpG mask excludes genomic CpG positions where TAPS 5mC→T conversion is indistinguishable from a true C>T SNP, substantially reducing false-positive SNP calls in TAPS data.

* **`<sample>/variants/<sample>.haplotypecaller.filtered-INDEL.vcf.gz`** / **`.tbi`**: PASS germline INDELs after hard filtering.

* **`<sample>/variants/<sample>.haplotypecaller.vcf.gz`** / **`.tbi`**: Unfiltered merged VCF from all 24 scattered HaplotypeCaller shards (chr1–22, X, Y). Use the `FILTER` field (`PASS` vs `DRAGENHardQUAL`) to select high-confidence variants.

* **Annotated copies (optional, `--annotate_vcf`)**: `<sample>.haplotypecaller.filtered-SNP.annotated.vcf.gz` and `...-INDEL.annotated.vcf.gz` in `variants/`, `<sample>.rastair_call.annotated.vcf.gz` in `methylation/`. See [Annotated VCFs](#annotated-vcfs-optional) below.

* **`<sample>/variants/<sample>.rastair_cpg_sites.bed.gz`** / **`.tbi`**: 3-column BED of all CpG positions called by rastair (REF CpGs + SNP-created NEW CpGs), with adjacent CpG strand pairs merged. This is the sample-specific CpG exclusion mask passed to GATK HaplotypeCaller via `--exclude-intervals`. Useful as a benchmarking reference for comparing called CpG sites against other callers or genome-wide CpG annotations.

#### Annotated VCFs (optional)

Written only with `--annotate_vcf true` (GRCh38 only, see [docs/usage.md](usage.md#annotated-vcfs---annotate_vcf)):
`<sample>.haplotypecaller.filtered-SNP.annotated.vcf.gz`, `<sample>.haplotypecaller.filtered-INDEL.annotated.vcf.gz` (in `variants/`) and `<sample>.rastair_call.annotated.vcf.gz` (in `methylation/`), each with a `.tbi` index. Each is a copy of its original VCF with annotation added; CHROM, POS, REF, ALT, QUAL, FILTER and the sample columns are unchanged, so `FILTER` values (`PASS`, `DRAGENHardQUAL`, ...) carry over. (`QUAL` keeps its value but may be written as `8.0` instead of `8`.)

> **Which records are in the annotated files:** `--annotate_chroms` (default `chr1`–`chr22`, `chrX`, `chrY`) applies to all three annotated files; records on other chromosomes are not included. The GATK annotated files then contain the same records as their originals (HaplotypeCaller only calls those chromosomes), unless you give a shorter list. **The rastair annotated file is always a subset of `rastair_call.vcf.gz`:** only records with a real ALT allele; CpG methylation rows (`ALT=.`) and variant records on other contigs are not included.

**Sources** (versions as configured on the JAX Elion cluster; other sites set their own in the genome config):

| Source | Used for | Version | Reference |
|--------|----------|---------|-----------|
| dbSNP | rsID in the `ID` column | build 151 (`dbsnp_151.vcf.gz`) | https://www.ncbi.nlm.nih.gov/snp/ |
| gnomAD | population allele frequency (`gnomAD_AF`) | `af-only-gnomad.hg38`, the GATK germline resource derived from gnomAD (the file header does not state the gnomAD release) | https://gnomad.broadinstitute.org/ |
| ClinVar | clinical significance | release 2026-10-04 (named in the `clinvar` path) | https://www.ncbi.nlm.nih.gov/clinvar/ |
| SnpEff | functional consequence on transcripts (`ANN`) | 5.1d, database `hg38` | https://pcingola.github.io/SnpEff/ (ANN field: https://pcingola.github.io/SnpEff/snpeff/inputoutput/) |
| dbNSFP | in-silico scores for missense variants | 4.2a | https://www.dbnsfp.org/ (scores: SIFT, [PolyPhen-2](http://genetics.bwh.harvard.edu/pph2/), [REVEL](https://sites.google.com/site/revelgenomics/), [CADD](https://cadd.gs.washington.edu/)) |

| Added field | Source | Meaning |
|-------------|--------|---------|
| `ID` column | dbSNP 151 | rsID (`.` if the allele is not in dbSNP) |
| `INFO/gnomAD_AF` | gnomAD allele-frequency resource (`af-only-gnomad.hg38`) | Allele frequency per ALT allele. **Absent means the allele was not found, not that the frequency is 0.** |
| `INFO/CLNSIG`, `CLNREVSTAT`, `CLNDN`, `CLNDISDB`, `GENEINFO`, `ALLELEID` | ClinVar (the release is the one named in the `clinvar` path of the genome config, for example `clinvar_20261004...`; the path is recorded in the VCF header's `##bcftools_annotateCommand` line) | Germline classification, review status, disease name and IDs, gene, ClinVar allele ID. A ClinVar match can lack `CLNSIG` when ClinVar holds no germline classification for it (`ALLELEID` is then still present). |
| `INFO/ANN`, `LOF`, `NMD` | SnpEff 5.1d, database `hg38` (RefSeq transcripts) | Functional consequence, one entry per transcript: Allele \| Annotation \| Annotation_Impact \| Gene_Name \| Gene_ID \| Feature_Type \| Feature_ID \| Transcript_BioType \| Rank \| HGVS.c \| HGVS.p \| cDNA.pos/length \| CDS.pos/length \| AA.pos/length \| Distance \| ERRORS/WARNINGS/INFO |
| `INFO/dbNSFP_SIFT_score`, `dbNSFP_SIFT_pred`, `dbNSFP_Polyphen2_HDIV_pred`, `dbNSFP_REVEL_score`, `dbNSFP_CADD_phred`, `dbNSFP_aapos`, `dbNSFP_Ensembl_transcriptid`, `dbNSFP_genename` | dbNSFP 4.2a via SnpSift 5.1d | In-silico scores for missense SNVs. These are **comma-separated lists with one value per Ensembl transcript**, in the same order as `dbNSFP_Ensembl_transcriptid`, although the VCF header declares them `Number=A`. SIFT: score below 0.05 is predicted damaging (`D`), otherwise tolerated (`T`). |

**Reading the values**

* **`gnomAD_AF`:** overall allele frequency in gnomAD (the resource has no per-ancestry breakdown). A rare variant has a small value; no value means the allele is not in the resource.
* **`CLNSIG`:** ClinVar's aggregate germline classification. The standard terms are `Pathogenic`, `Likely_pathogenic`, `Uncertain_significance`, `Likely_benign` and `Benign` (ACMG/AMP terms, see https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/); others such as `risk_factor` also occur. `CLNREVSTAT` says how well supported the classification is (see https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/).
* **`ANN` impact (third field):** `HIGH` = probably disrupts the protein (for example `stop_gained`, `frameshift_variant`); `MODERATE` = non-disruptive but may change protein function (for example `missense_variant`); `LOW` = mostly harmless (for example `synonymous_variant`); `MODIFIER` = mostly non-coding, effect hard to predict. With several transcripts, SnpEff lists the effects sorted by putative impact first. Definitions: https://pcingola.github.io/SnpEff/snpeff/inputoutput/.
* **dbNSFP scores** (from the dbNSFP 4.2a readme):
  * `dbNSFP_SIFT_score`: 0 to 1, **smaller = more likely damaging**. `dbNSFP_SIFT_pred` is `D` (damaging) below 0.05, otherwise `T` (tolerated).
  * `dbNSFP_Polyphen2_HDIV_pred`: `D` probably damaging, `P` possibly damaging, `B` benign.
  * `dbNSFP_REVEL_score`: 0 to 1, an ensemble of 13 scores, **larger = more likely damaging**.
  * `dbNSFP_CADD_phred`: phred-like rank score, **larger = more likely damaging**.
  * `dbNSFP_aapos`: amino-acid position (`-1` for splice-site variants).
  * These are scores for missense (non-synonymous) variants; other variants have no value.
* **Licences:** the dbNSFP readme states that REVEL scores are free for non-commercial use and CADD scores are free for academic, non-commercial applications. Check these licences before commercial use.
* **Use:** these are annotations from public resources at the versions above, for research use. They are lookups and predictions, not a clinical interpretation.

The exact commands and resource paths used are recorded in the VCF header (`##bcftools_annotateCommand`, `##SnpEffCmd`, `##SnpSiftCmd`).

---

### 9. Key QC & Metric Files

| File | Description |
| :--- | :--- |
| **`report/5-baseTAPS_multiqc_report.html`** | Integrated QC report — start here ([interpretation guide](./interpret_multiqc.md)) |
| **`<s>/qc/alignment/<s>.postdedup.flagstat`** | Samtools flagstat on final consensus BAM |
| **`<s>/qc/alignment/<s>.prededup.flagstat`** | Samtools flagstat before UMI deduplication |
| **`<s>/qc/alignment/<s>.lambda_negCtrl.flagstat`** | Flagstat for lambda spike-in reads |
| **`<s>/qc/alignment/<s>.puc19_posCtrl.flagstat`** | Flagstat for pUC19 spike-in reads |
| **`<s>/qc/coverage/<s>.mosdepth.summary.txt`** | Mean, median, and percentile WGS coverage |
| **`<s>/qc/duplex/<s>.duplex_seq_metrics.duplex_yield_metrics.txt`** | Duplex yield: input reads, consensus reads, fold-reduction |
| **`<s>/qc/duplex/<s>.duplex_seq_metrics.duplex_qc.pdf`** | fgbio duplex QC summary PDF |
| **`<s>/qc/variant_call/<s>.bcftools_stats.txt`** | bcftools stats on PASS SNP VCF |
| **`<s>/qc/variant_call/<s>.dragstr_model.txt`** | DRAGEN STR model calibration parameters |
| **`<s>/methylation/<s>.rastair_mbias.html`** | M-bias report (M-bias, V-bias, GC/CpG density bias) |
| **`<s>/qc/fastqc/<s>_1_fastqc.html`** | FastQC quality report for R1 raw reads |

NOTE: `<s>` represents the sample ID wildcard.

---

### `duplex_seq_metrics.duplex_qc.pdf` — Plot Descriptions

Generated by [fgbio CollectDuplexSeqMetrics](https://fulcrumgenomics.github.io/fgbio/tools/latest/CollectDuplexSeqMetrics.html). Contains **8 pages** for WGS production samples.

**Terminology:**

| Term | Definition |
|------|-----------|
| By Coord+Strand (CS) | Families grouped by genomic coordinate + strand, ignoring UMI — the "raw" grouping |
| SS Families | Single-strand families — each UMI orientation (AB, BA) counted separately |
| DS Families | Double-strand families — one entry per source molecule with reads on both strands |
| Duplexes – Ideal | DS families where at least 1 read was observed on each strand (AB ≥ 1 AND BA ≥ 1) |
| Duplexes – Actual | DS families that passed fgbio quality filters and produced a valid duplex consensus read |

---

**Page 1 — Family Size Distributions**

Three curves (CS, DS, SS) plotted as count vs. family size on a log2-scaled x-axis.

At production WGS depth (300–600M reads), all three curves nearly overlap and show very similar decay — the CS curve no longer dominates over DS/SS as it does at shallow depths. At shallow test depths (< 30M reads), the CS curve is much higher than DS/SS.

- **Expected:** Peak at size 1 (e.g., ~4×10^8 for a 583M-read sample), steeply declining. Virtually no families beyond size 8–16.
- **Investigate if:** A secondary bump at high family sizes — indicates PCR jackpotting (a subset of molecules over-amplified into very large families).

---

**Page 2 — Cumulative Family Size Distributions**

Fraction of families with ≥ each family size (same three curves). At production WGS depth, all three curves are nearly identical and drop sharply: typically ~10% of DS families have ≥ 2 reads per strand at 30–50× WGS coverage.

- **Use:** The y-value at size = N on the DS curve tells you the yield fraction that survives a `--min-reads N` filter in FilterConsensusReads. Useful for evaluating trade-offs before tightening filtering thresholds.

---

**Page 3 — Duplex Tag Family Size Distribution (all families)**

2D heatmap: AB-strand reads (x-axis) vs. BA-strand reads (y-axis). Color scale (log) from dark green (millions of families) to light blue (~1 family). Includes SS families in the BA=0 row.

- **AB axis** typically extends into the thousands for large samples (e.g., 0–6000 for the 583M-read sample) due to rare jackpotted families with extreme AB read counts.
- **BA axis** is much shorter (0–5 for the 583M-read sample); BA reads are far sparser than AB reads even at high WGS depth.
- **For shallow samples (< 30M reads):** BA axis is compressed to 0–1 and the entire heatmap appears as uniform-colored vertical bands, with no individual cells distinguishable.
- **Expected:** Nearly all mass as a dense strip at BA=0 (SS families). DS families (BA > 0) appear as very faint thin horizontal lines. The plot looks mostly empty aside from the BA=0 strip.
- **Investigate if:** Multiple dark-green outlier columns at large AB (> 1000) values — indicates jackpotted molecules with thousands of PCR copies from a single source.

---

**Page 4 — Duplex Tag Family Size Distribution (DS families only)**

Same heatmap filtered to families where AB > 0 AND BA > 0. The AB axis contracts significantly (to ~150 for the 583M-read sample) because jackpotted SS families are excluded.

- **Expected:** Darkest cells at small AB and BA (AB=1–5, BA=1–2). The modal family has just 1–2 reads per strand. Sparse outliers at AB=50–150 are normal at 583M reads.
- **For shallow test samples:** This plot will be completely blank — essentially no DS families exist at 10M reads. This is expected, not a failure.
- **Investigate if** (production samples): Distribution skewed with most mass at AB >> BA, indicating systematic strand-biased PCR amplification.

---

**Page 5 — Duplex Yield by Input Read Pairs**

Area chart: DS Families (blue), Duplexes – Actual (green), Duplexes – Ideal (red) as a function of total input read pairs.

- **Expected:** DS Families (blue) grows linearly — the physical library is not exhausted. Duplexes–Ideal (red) and Duplexes–Actual (green) are thin slivers near the bottom of the plot.
- **Reading the three layers:**
  - **Blue (DS Families):** all double-strand families — molecules where reads from both strands were captured
  - **Red (Duplexes – Ideal):** subset with ≥ 1 read per strand — theoretical maximum for duplex consensus calling
  - **Green (Duplexes – Actual):** subset passing fgbio quality filters — the final usable duplex consensus count
  - Large gap between blue and red: most DS families have too few reads on one strand to qualify even in theory
  - Small gap between red and green: a fraction of Ideal duplexes are removed by quality filtering
- **Investigate if:** DS Families curve starts to plateau — additional sequencing will yield diminishing unique molecules, indicating library saturation.

---

**Page 6 — Ratio of Actual vs. Ideal Duplex Yield**

Single curve (Actual / Ideal duplexes) plotted against input read pairs.

- **Expected:** A flat horizontal line, stable across all read depths. For the 583M-read GT26-1111 sample, the ratio is flat at ~0.12 (12%). A flat ratio means quality filtering is depth-independent — adding more reads increases Ideal and Actual proportionally.
- **Investigate if:** Ratio decreases at higher read depths — may indicate quality degradation at late sequencing cycles or library-specific chemistry issues.

---

**Page 7 — UMI Representation**

Scatter plot: x = observations of each UMI sequence in raw reads; y = unique tag families observed for that UMI.

- **Expected:** Tight linear relationship from origin. Each UMI is used proportionally to how often it was sequenced — indicating good UMI library diversity. The 583M-read sample shows a tight line up to ~1.5×10^7 observations with a few higher outliers.
- **Investigate if:** A population of points well above the main diagonal — indicates UMI collision: the same UMI sequence was shared by multiple source molecules (UMI library complexity exhausted). This would inflate apparent DS family counts.

---

**Page 8 — Read Distribution Among Families**

Three curves (CS, DS, SS) plotted as **total reads allocated to families of size N** vs. family size (log2 scaled). Similar shape to Page 1 but y-axis is read count, not family count.

- **Key difference from Page 1:** Large families are weighted more heavily here. A jackpotted family of size 1000 appears as 1 count in Page 1 but contributes 1000 reads to Page 8. If jackpotting is consuming a meaningful fraction of your sequencing, Page 8 will show a tail or secondary bump at high family sizes that is absent in Page 1.
- **Expected:** Dominates at size 1. No secondary peak at large family sizes. At 583M reads, Page 8 looks visually similar to Page 1, confirming jackpotted families are not a significant fraction of total reads.

---

## 5baseTAPS Pipeline Generated Files

| Output File | Folder | Description |
|:------------|:-------|:------------|
| `<s>.cons.filtered.bam` | `<s>/bam/` | **Primary BAM** — duplex/molecular consensus reads aligned to reference genome |
| `<s>.cons.filtered.bam.bai` | `<s>/bam/` | BAM index |
| `<s>.rastair_call.bed.gz` | `<s>/methylation/` | **Primary methylation** — per-CpG beta values and read counts; coverage ≥ 1 |
| `<s>.rastair_call.vcf.gz` | `<s>/methylation/` | Combined methylation + SNP calls (CpG records + genome-wide SNP records) |
| `<s>.rastair_perread.bed.gz` | `<s>/methylation/` | Per-read CpG methylation assignments |
| `<s>.rastair_perread.bed.gz.tbi` | `<s>/methylation/` | Tabix index for per-read BED |
| `<s>.rastair_methylkit.txt.gz` | `<s>/methylation/` | Methylation in methylKit R format for differential analysis |
| `<s>.rastair_mbias.html` | `<s>/methylation/` | M-bias report — M-bias, V-bias, GC/CpG density bias |
| `<s>.methylation_summary.tsv` | `<s>/methylation/` | Per-chromosome and genome-wide CpG coverage + global methylation % |
| `<s>.lambda_negCtrl.methylation_summary.tsv` | `<s>/methylation/` | Lambda spike-in methylation summary (conversion efficiency QC) |
| `<s>.puc19_posCtrl.methylation_summary.tsv` | `<s>/methylation/` | pUC19 spike-in methylation summary (positive control QC) |
| `<s>.haplotypecaller.filtered-SNP.vcf.gz{,.tbi}` | `<s>/variants/` | **Primary SNPs** — PASS germline SNPs, CpG-masked, GATK hard-filtered |
| `<s>.haplotypecaller.filtered-INDEL.vcf.gz{,.tbi}` | `<s>/variants/` | PASS germline INDELs after hard filtering |
| `<s>.haplotypecaller.vcf.gz{,.tbi}` | `<s>/variants/` | Unfiltered merged GATK HC VCF (all variants, use FILTER field) |
| `<s>.rastair_cpg_sites.bed.gz{,.tbi}` | `<s>/variants/` | Sample-specific CpG site mask (REF + NEW CpGs called by rastair) used for GATK HC exclusion |
| `<s>_1_fastqc.html{,.zip}` | `<s>/qc/fastqc/` | FastQC quality report for R1 raw reads |
| `<s>_2_fastqc.html{,.zip}` | `<s>/qc/fastqc/` | FastQC quality report for R2 raw reads |
| `<s>.prededup.flagstat` | `<s>/qc/alignment/` | Samtools flagstat before UMI deduplication |
| `<s>.postdedup.flagstat` | `<s>/qc/alignment/` | Samtools flagstat on final consensus BAM |
| `<s>.lambda_negCtrl.flagstat` | `<s>/qc/alignment/` | Flagstat for lambda spike-in reads |
| `<s>.puc19_posCtrl.flagstat` | `<s>/qc/alignment/` | Flagstat for pUC19 spike-in reads |
| `<s>.mosdepth.summary.txt` | `<s>/qc/coverage/` | Mean, median, and percentile WGS coverage by chromosome |
| `<s>.mosdepth.global.dist.txt` | `<s>/qc/coverage/` | Cumulative coverage distribution (used in MultiQC mosdepth plot) |
| `<s>.grouped-family-sizes.txt` | `<s>/qc/duplex/` | UMI family size frequency table |
| `<s>.duplex_seq_metrics.duplex_yield_metrics.txt` | `<s>/qc/duplex/` | Duplex yield: reads in, consensus reads out, fold-reduction, duplex yield % |
| `<s>.duplex_seq_metrics.family_sizes.txt` | `<s>/qc/duplex/` | Raw UMI family size histogram (all reads) |
| `<s>.duplex_seq_metrics.duplex_family_sizes.txt` | `<s>/qc/duplex/` | Duplex-only family size histogram |
| `<s>.duplex_seq_metrics.umi_counts.txt` | `<s>/qc/duplex/` | UMI sequence frequency table |
| `<s>.duplex_seq_metrics.duplex_umi_counts.txt` | `<s>/qc/duplex/` | Duplex UMI pair frequency table |
| `<s>.duplex_seq_metrics.duplex_qc.pdf` | `<s>/qc/duplex/` | fgbio duplex QC summary PDF |
| `<s>.bcftools_stats.txt` | `<s>/qc/variant_call/` | bcftools stats on PASS SNP VCF (Ti/Tv, SNP/INDEL counts, quality) |
| `<s>.dragstr_model.txt` | `<s>/qc/variant_call/` | DRAGEN STR model calibration file (CalibrateDragstrModel) |
| `<s>.cons.stats` | `<s>/qc/alignment/` | Samtools stats on final consensus BAM (MAPQ ≥ 20; insert size distribution for MultiQC) |
| `5-baseTAPS_multiqc_report.html` | `report/` | **Integrated QC report** — mapping, duplex, methylation, coverage, variants |

NOTE: `<s>` is a placeholder for the sample ID.

---

