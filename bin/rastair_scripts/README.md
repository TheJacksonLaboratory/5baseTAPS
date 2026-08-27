# JAX Rastair Mbias Scripts

These scripts override the built-in mbias R scripts shipped inside the rastair container image.
The Nextflow pipeline passes `--r-script-dir` to `rastair mbias` pointing here, so `mbias.R`
and `QC_report_jax.Rmd` are used in place of the container defaults.

## Files

- **`mbias.R`** — Entry point called by rastair. Parses CLI arguments and calls
  `rmarkdown::render()` on `QC_report_jax.Rmd`.
- **`QC_report_jax.Rmd`** — RMarkdown report generating M-bias, V-bias, GC/CpG bias plots
  and per-chromosome cutoff files.

## JAX Modifications vs. Container Original

### 1. Vectorized position parsing (`parse_positions_vec`)

**Problem:** The original `parse_positions()` function was called row-by-row inside a
`data.table` `by=` group, resulting in one `strsplit` call per read (~300 million R function
invocations for a 570M-read sample). On a 583M-read whole-genome sample this step ran for
over 20 hours without completing.

**Fix:** Replaced with `parse_positions_vec()`, which operates on full columns using
vectorized `strsplit` + `unlist` — a single C-level pass per column type. Performance on
the same 570M-read dataset: **~34 minutes** (vs. 20+ hours).

### 2. Memory-efficient column handling

Inside `parse_positions_vec`, after all `strsplit` passes are complete:
- Position columns (`mod_cpgs`, `unmod_cpgs`, `snp_cpgs`, `mod_denovos`, `unmod_denovos`)
  are removed from `dt` in-place (`dt[, (pos_cols) := NULL]`) before expanding to the
  billion-row result table, freeing ~30GB of string data.
- `rm(parts); gc()` is called before the final `cbind` to release intermediate data.tables.
- `meta_cols` is reduced to only columns needed downstream:
  - Without vbias (`--no-vbias`): `c("chr", "flag")` — saves ~80GB at the cbind step
    (1.3B rows × 2 cols vs. 9 cols).
  - With vbias: `c("chr", "start", "end", "read_id", "orientation", "insert_size", "flag")`.
  - `mapq` and `read_length` are never used after this point and are dropped in both cases.

### 3. Deduplicated binomial CI computation (`binom_ci_dedup`)

**Problem:** The original called `binom.test()` once per row of `mbias_table` (~91K rows
for a whole-genome sample), many of which share the same `(mod, n)` integer pair.

**Fix:** `binom_ci_dedup()` computes `binom.test` only on unique `(mod, n)` pairs and
joins results back, significantly reducing redundant computation.

### 4. Multi-threaded data.table aggregation

`setDTthreads(max(1L, as.integer(params$threads)))` is called after library load so the
data.table aggregation step uses all allocated CPUs (typically 2 in the Nextflow config).

### 5. Stderr timing messages

Key steps emit timestamped progress to stderr (visible in SLURM `.e` logs, not in HTML):
- `fread` completion: reads loaded, chromosome count, elapsed seconds
- `parse_positions_vec` completion: CpG observations, elapsed seconds
- `binom CIs` completion: elapsed seconds
- `plot_mbias`: per-chromosome progress

`cat(..., file=stderr())` is used instead of `message()` to bypass knitr's output capture
and prevent timing messages from appearing in the HTML report.

### 6. OT-first strand ordering in cutoff files

The original `write.table()` call had data-dependent row ordering (OT/OB order varied
by chromosome). Fixed with `out_table[order(strand, decreasing=TRUE), ]` to enforce
OT before OB consistently across all chromosomes.

### 7. Report title

Changed from `"QC Report"` to `"Rastair Methylation Bias Report"`.

## Performance (583M reads, GRCh38, 2 CPUs, 184GB RAM)

| Step | Time |
|------|------|
| fread | 18 min |
| CpG + flag filter | 7.5 min |
| Dedup | 5.5 min |
| parse_positions_vec | 34 min |
| mbias_table aggregation | <1 min |
| binom CIs (deduped) | 2.3 min |
| plot_mbias (191 chromosomes) | ~2 min |
| **Total** | **~70 min** |

The original script did not complete after 20+ hours on the same dataset.
