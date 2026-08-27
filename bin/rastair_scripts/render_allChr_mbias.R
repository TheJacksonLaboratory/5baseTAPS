#!/usr/bin/env Rscript
# Stage 2: gather per-chromosome RDS files and render the full HTML report.
# Mirrors the interface of mbias.R but takes --rds-dir instead of --bed.

suppressMessages({
  library(rmarkdown)
  library(argparser)
})

rmd_file <- "allChr_mbias_report_jax.Rmd"

get_script_dir <- function(command_line_args) {
  command_line = paste(command_line_args, collapse=" ")
  dir = gsub(".*--file=(.*) --args.*", "\\1", as.character(command_line), perl = TRUE)
  return(dir)
}

parser <- arg_parser("Gather per-chr RDS files and render M-bias HTML report")
parser <- add_argument(parser, "--rds-dir",
                       help="Directory containing chrN_mbias_table.rds files from stage 1")
parser <- add_argument(parser, "--output-prefix",
                       help="Output directory",
                       default=".")
parser <- add_argument(parser, "--threads",
                       help="Number of threads for data.table",
                       default=1)
parser <- add_argument(parser, "--vbias",
                       help="Render V-bias plots (requires --vbias in stage 1 chr jobs)",
                       flag=TRUE)
parser <- add_argument(parser, "--gc",
                       help="Render GC/CpG bias plots",
                       flag=TRUE)
parser <- add_argument(parser, "--vcf",
                       help="Path to rastair_call.vcf.gz for GC bias",
                       default="")
parser <- add_argument(parser, "--reference",
                       help="Path to reference FASTA for GC bias",
                       default="")
parser <- add_argument(parser, "--bcftools",
                       help="Path to bcftools binary",
                       default="bcftools")
parser <- add_argument(parser, "--sample-name",
                       help="Sample name shown in the report title",
                       default="")

args <- parse_args(parser)

if (is.na(args$rds_dir) || is.null(args$rds_dir) || !dir.exists(args$rds_dir)) {
  cat("Error: --rds-dir must point to a directory containing chrN_mbias_table.rds files\n")
  quit(status=1)
}

params_list <- list(
  rds_dir     = normalizePath(args$rds_dir),
  output_dir  = normalizePath(args$output_prefix),
  threads     = args$threads,
  plot_vbias  = args$vbias,
  plot_gc     = args$gc,
  input_vcf   = if (nchar(args$vcf) > 0) args$vcf else NA,
  reference   = if (nchar(args$reference) > 0) args$reference else NA,
  bcftools    = args$bcftools,
  sample_name = if (!is.na(args$sample_name)) args$sample_name else ""
)

script_dir <- get_script_dir(commandArgs())
script_dir <- dirname(normalizePath(script_dir))
rmd_file   <- paste0(script_dir, "/", rmd_file)

cat("Rmd_file:", rmd_file, "\n")
if (!file.exists(rmd_file)) {
  cat("Error: RMarkdown file not found:", rmd_file, "\n")
  quit(status=1)
}

output_file <- paste0(normalizePath(args$output_prefix), "/qc_report.html")
str(params_list)
cat("Rendering RMarkdown document...\n")
tryCatch({
  rmarkdown::render(
    input           = rmd_file,
    output_file     = output_file,
    params          = params_list,
    intermediates_dir = getwd(),
    clean           = TRUE,
    quiet           = FALSE
  )
  cat("Report successfully generated:", output_file, "\n")
}, error = function(e) {
  cat("Error rendering RMarkdown document:\n")
  cat(conditionMessage(e), "\n")
  quit(status=1)
})

cat("Done!\n")
