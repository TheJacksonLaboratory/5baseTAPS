#!/usr/bin/env python3
"""
Standalone 5-base TAPS QC report generator.

Discovers all samples in a completed (or partially completed) pipeline output
directory, runs all per-sample metric scripts that have sufficient inputs, and
then runs MultiQC on the results.

Usage:
    run_taps_qc.py <outdir> [options]

    <outdir>      Pipeline output directory (contains <sample>/ subdirectories
                  each with methylation/, variants/, qc/ subdirs)

Options:
    --mqc-outdir  DIR   Where to write *_mqc.tsv files (default: <outdir>/report/mqc)
    --multiqc     PATH  Path to multiqc executable (default: auto-detect from PATH)
    --no-multiqc        Generate *_mqc.tsv files only; skip MultiQC
    --samples     S...  Only process these sample IDs (space-separated, default: all)
    --dry-run           Print what would be run without executing

Example:
    run_taps_qc.py /path/to/outdir
    run_taps_qc.py /path/to/outdir --samples SAMPLE1 --no-multiqc
"""

import sys
import os
import argparse
import subprocess
import shutil
from pathlib import Path
from glob import glob


# ── File discovery ────────────────────────────────────────────────────────────

def find_samples(outdir):
    """Discover sample IDs from <outdir>/<sample>/qc/ subdirectories."""
    samples = sorted(d.name for d in outdir.iterdir() if d.is_dir() and (d / "qc").is_dir())
    if not samples:
        sys.exit(f"ERROR: no sample directories (with qc/ subdirectory) found in {outdir}")
    return samples


def first(pattern):
    """Return first glob match or None."""
    matches = glob(str(pattern))
    return matches[0] if matches else None


def find_files(outdir, sample):
    """
    Locate all possible QC input files for a sample.
    Returns a dict of file-role -> Path (or None if not found).
    """
    o = outdir
    s = sample
    q = o / s / "qc"
    al = q / "alignment"
    co = q / "coverage"
    du = q / "duplex"
    vc = q / "variant_call"
    me = o / s / "methylation"   # methylation summaries co-located with call outputs

    def g(pattern):
        return first(pattern)

    return {
        # Mapping / spike-in composition
        "prededup_flagstat":   g(al / f"{s}.prededup.flagstat"),
        "lambda_flagstat":     g(al / f"{s}.lambda_negCtrl.flagstat"),
        "puc19_flagstat":      g(al / f"{s}.puc19_posCtrl.flagstat"),
        # Duplex metrics
        "postdedup_flagstat":  (
            g(al / f"{s}.postdedup.flagstat") or
            g(al / f"{s}.flagstat")
        ),
        "duplex_yield":        g(du / f"{s}.duplex_seq_metrics.duplex_yield_metrics.txt"),
        "family_sizes":        g(du / f"*.family_sizes.txt"),
        # WGS coverage
        "mosdepth_summary":    g(co / f"{s}.mosdepth.summary.txt"),
        "mosdepth_dist":       g(co / f"{s}.mosdepth.global.dist.txt"),
        # Methylation summaries (in methylation/ alongside call outputs)
        "main_methyl":         g(me / f"{s}.methylation_summary.tsv"),
        "lambda_methyl":       g(me / f"{s}.lambda_negCtrl.methylation_summary.tsv"),
        "puc19_methyl":        g(me / f"{s}.puc19_posCtrl.methylation_summary.tsv"),
        # Variant calling
        "bcftools_stats":      g(vc / f"{s}.bcftools_stats.txt"),
    }


# ── Script locations ──────────────────────────────────────────────────────────

SCRIPT_DIR = Path(__file__).parent.resolve()


def script(name):
    return str(SCRIPT_DIR / name)


# ── Per-sample metric runners ─────────────────────────────────────────────────

def run_mapping(sample, files, mqc_outdir, dry_run):
    required = ["prededup_flagstat", "lambda_flagstat", "puc19_flagstat"]
    if not all(files.get(k) for k in required):
        missing = [k for k in required if not files.get(k)]
        print(f"  [SKIP mapping] {sample}: missing {missing}")
        return False
    cmd = [
        sys.executable, script("mapping_metrics.py"),
        sample, files["prededup_flagstat"], files["lambda_flagstat"], files["puc19_flagstat"],
    ]
    return _run(cmd, cwd=mqc_outdir, dry_run=dry_run, label=f"mapping {sample}")


def run_methyl(sample, files, mqc_outdir, dry_run, genome="Sample"):
    required = ["main_methyl", "lambda_methyl", "puc19_methyl"]
    if not all(files.get(k) for k in required):
        missing = [k for k in required if not files.get(k)]
        print(f"  [SKIP methyl] {sample}: missing {missing}")
        return False
    cmd = [
        sys.executable, script("methyl_metrics.py"),
        sample, files["main_methyl"], files["lambda_methyl"], files["puc19_methyl"], genome,
    ]
    return _run(cmd, cwd=mqc_outdir, dry_run=dry_run, label=f"methyl {sample}")


def run_duplex(sample, files, mqc_outdir, dry_run):
    required = ["prededup_flagstat", "postdedup_flagstat", "duplex_yield", "family_sizes"]
    if not all(files.get(k) for k in required):
        missing = [k for k in required if not files.get(k)]
        print(f"  [SKIP duplex] {sample}: missing {missing}")
        return False
    duplex_dir = Path(files["duplex_yield"]).parent
    duplex_files = list(duplex_dir.glob("*duplex_seq_metrics*.txt"))
    cmd = [
        sys.executable, script("duplex_metrics.py"),
        sample, files["prededup_flagstat"], files["postdedup_flagstat"],
    ] + [str(f) for f in duplex_files]
    return _run(cmd, cwd=mqc_outdir, dry_run=dry_run, label=f"duplex {sample}")


def run_coverage(sample, files, mqc_outdir, dry_run):
    required = ["mosdepth_summary", "mosdepth_dist", "main_methyl"]
    if not all(files.get(k) for k in required):
        missing = [k for k in required if not files.get(k)]
        print(f"  [SKIP coverage] {sample}: missing {missing}")
        return False
    cmd = [
        sys.executable, script("wgs_coverage_metrics.py"),
        sample, files["mosdepth_summary"], files["mosdepth_dist"], files["main_methyl"],
    ]
    return _run(cmd, cwd=mqc_outdir, dry_run=dry_run, label=f"coverage {sample}")


def run_vc(sample, files, mqc_outdir, dry_run):
    if not files.get("bcftools_stats"):
        print(f"  [SKIP vc] {sample}: missing bcftools_stats")
        return False
    cmd = [
        sys.executable, script("vc_metrics.py"),
        sample, files["bcftools_stats"],
    ]
    return _run(cmd, cwd=mqc_outdir, dry_run=dry_run, label=f"vc {sample}")


def _run(cmd, cwd, dry_run, label):
    if dry_run:
        print(f"  [DRY-RUN] {label}: {' '.join(cmd)}")
        return True
    try:
        result = subprocess.run(cmd, cwd=str(cwd), check=True, capture_output=True, text=True)
        print(f"  [OK] {label}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"  [FAIL] {label}:\n{e.stderr.strip()}")
        return False


# ── MultiQC ───────────────────────────────────────────────────────────────────

def run_multiqc(outdir, mqc_outdir, multiqc_bin, config, logo, dry_run):
    if not multiqc_bin:
        # Try known SIF locations, then fall back to PATH
        sif_candidates = [
            "/flashscratch/nf-JAX-5base/singularity/cache/depot.galaxyproject.org-singularity-multiqc-1.35--pyhdfd78af_1.img",
            "/gt/research_development/singularity/cache/depot.galaxyproject.org-singularity-multiqc-1.35--pyhdfd78af_1.img",
        ]
        sif = next((Path(p) for p in sif_candidates if Path(p).exists()), None)
        if sif:
            singularity = shutil.which("singularity") or next(
                iter(sorted(Path("/cm/local/apps/apptainer").glob("*/bin/singularity"))), None
            ) or "singularity"
            multiqc_bin = f"{singularity} exec {sif} multiqc"
        else:
            multiqc_bin = shutil.which("multiqc") or "multiqc"

    # Split multiqc_bin in case it's a compound command (e.g. "singularity exec ... multiqc")
    import shlex
    multiqc_cmd_prefix = shlex.split(multiqc_bin)

    config_path = Path(config).resolve() if config else None
    logo_path = Path(logo).resolve() if logo else None

    if config_path and not config_path.is_file():
        sys.exit(f"ERROR: MultiQC config not found: {config_path}")
    if logo_path and not logo_path.is_file():
        sys.exit(f"ERROR: MultiQC logo not found: {logo_path}")

    config_args = ["-c", str(config_path)] if config_path else []
    logo_args = ["--cl-config", f'custom_logo: "{logo_path}"'] if logo_path else []

    multiqc_outdir = outdir / "report"

    cmd = multiqc_cmd_prefix + [
        str(outdir),       # single dir — searches recursively for *_mqc.tsv, flagstat, mosdepth, etc.
        "--outdir", str(multiqc_outdir),
        "--filename", "5-baseTAPS_multiqc_report",
        "--force",
    ] + config_args + logo_args

    if dry_run:
        print(f"\n[DRY-RUN] MultiQC: {' '.join(cmd)}")
        return
    print(f"\nRunning MultiQC...")
    try:
        subprocess.run(cmd, check=True)
        report = multiqc_outdir / "5-baseTAPS_multiqc_report.html"
        print(f"\nReport: {report}")
    except subprocess.CalledProcessError as e:
        print(f"[FAIL] MultiQC exited with {e.returncode}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Standalone TAPS QC report generator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("outdir", help="Pipeline output directory")
    parser.add_argument("--mqc-outdir", default=None,
                        help="Where to write *_mqc.tsv files (default: <outdir>/report/mqc)")
    parser.add_argument("--multiqc", default=None,
                        help="Path to multiqc executable (default: auto-detect)")
    parser.add_argument("--config", default=None,
                        help="Path to MultiQC config YAML")
    parser.add_argument("--logo", default=None,
                        help="Path to custom MultiQC logo")
    parser.add_argument("--no-multiqc", action="store_true",
                        help="Generate *_mqc.tsv files only, skip MultiQC")
    parser.add_argument("--samples", nargs="+", default=None,
                        help="Only process these sample IDs")
    parser.add_argument("--genome", default="Sample",
                        help="Reference genome name for methylation column label (e.g. CHM13, GRCh38)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print commands without executing")
    args = parser.parse_args()

    outdir = Path(args.outdir).resolve()
    if not outdir.is_dir():
        sys.exit(f"ERROR: {outdir} does not exist")

    mqc_outdir = Path(args.mqc_outdir).resolve() if args.mqc_outdir else outdir / "report" / "mqc"
    if not args.dry_run:
        mqc_outdir.mkdir(parents=True, exist_ok=True)

    samples = args.samples or find_samples(outdir)
    print(f"Output dir : {outdir}")
    print(f"MQC outdir : {mqc_outdir}")
    print(f"Samples    : {', '.join(samples)}\n")

    n_ok = n_skip = n_fail = 0
    for sample in samples:
        print(f"── {sample}")
        files = find_files(outdir, sample)
        runners = [
            lambda s, f, o, d: run_mapping(s, f, o, d),
            lambda s, f, o, d: run_methyl(s, f, o, d, genome=args.genome),
            lambda s, f, o, d: run_duplex(s, f, o, d),
            lambda s, f, o, d: run_coverage(s, f, o, d),
            lambda s, f, o, d: run_vc(s, f, o, d),
        ]
        for runner in runners:
            result = runner(sample, files, mqc_outdir, args.dry_run)
            if result is True:
                n_ok += 1
            elif result is False:
                n_skip += 1

    print(f"\nDone: {n_ok} OK, {n_skip} skipped (missing inputs)")

    if not args.no_multiqc:
        run_multiqc(outdir, mqc_outdir, args.multiqc, args.config, args.logo, args.dry_run)


if __name__ == "__main__":
    main()
