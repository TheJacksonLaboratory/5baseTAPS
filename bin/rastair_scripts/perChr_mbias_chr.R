#!/usr/bin/env Rscript
# Stage 1: per-chromosome data transform.
# Loads one chromosome via tabix, runs parse_positions_vec, saves mbias_table RDS.

suppressMessages({
  library(data.table)
  library(argparser)
})

parser <- arg_parser("Per-chromosome M-bias data transform")
parser <- add_argument(parser, "--bed",          help="Input per-read bed.gz (tabix indexed)")
parser <- add_argument(parser, "--chr",          help="Chromosome to process, e.g. chr1")
parser <- add_argument(parser, "--output",       help="Output RDS file path", default="mbias_table.rds")
parser <- add_argument(parser, "--tabix-path",   help="Path to tabix",        default="tabix")
parser <- add_argument(parser, "--include-flag", help="Include bitflag",      type="integer", default=3)
parser <- add_argument(parser, "--exclude-flag", help="Exclude bitflag",      type="integer", default=3852)
parser <- add_argument(parser, "--threads",      help="data.table threads",   type="integer", default=1)
parser <- add_argument(parser, "--wgbs",         help="Flip mod/unmod for WGBS", flag=TRUE)
parser <- add_argument(parser, "--vbias",        help="Compute V-bias aggregates for stage 2", flag=TRUE)

args <- parse_args(parser)
setDTthreads(max(1L, args$threads))

.ts <- function(msg) cat(format(Sys.time(), "[%H:%M:%S]"), " [", args$chr, "] ", msg, "\n", file=stderr(), sep="")

parse_positions_vec <- function(dt, plot_vbias=FALSE) {
  if (plot_vbias) {
    meta_cols <- c("chr", "start", "end", "read_id", "orientation", "insert_size", "flag")
  } else {
    meta_cols <- c("chr", "flag")
  }
  col_types  <- list(mod_cpgs="M", unmod_cpgs="U", snp_cpgs="S", mod_denovos="m", unmod_denovos="u")
  n <- nrow(dt)
  parts <- lapply(names(col_types), function(col) {
    splits <- strsplit(as.character(dt[[col]]), ",", fixed=TRUE)
    lens   <- lengths(splits)
    pos    <- suppressWarnings(as.integer(unlist(splits, use.names=FALSE)))
    ridx   <- rep(seq_len(n), lens)
    keep   <- !is.na(pos)
    if (!any(keep)) return(NULL)
    data.table(row_idx=ridx[keep], pos_in_read=pos[keep], type=col_types[[col]])
  })
  dt[, (names(col_types)) := NULL]
  gc()
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0) return(dt[integer(0), ..meta_cols])
  result <- rbindlist(parts)
  rm(parts); gc()
  cbind(dt[result$row_idx, ..meta_cols], result[, .(pos_in_read, type)])
}

.ts("loading...")
.t0 <- proc.time()
data <- fread(cmd=paste(args$tabix_path, "-h", normalizePath(args$bed), args$chr),
              header=TRUE, stringsAsFactors=FALSE)
setnames(data, 1, "chr")
data[, chr         := as.factor(chr)]
data[, orientation := as.factor(orientation)]
.ts(paste0("loaded: ", nrow(data), " reads (", round((proc.time()-.t0)["elapsed"], 1), "s)"))

if (args$wgbs) {
  data[, num_mod := num_cpg - num_mod]
  setnames(data, c("mod_cpgs","unmod_cpgs","mod_denovos","unmod_denovos"),
                 c("unmod_cpgs","mod_cpgs","unmod_denovos","mod_denovos"))
}

data <- data[num_cpg > 0]
data <- data[bitwAnd(flag, args$include_flag) == args$include_flag]
data <- data[bitwAnd(flag, args$exclude_flag) == 0]
.ts(paste0("after CpG + flag filter: ", nrow(data), " reads"))

read_len <- max(data$read_length, na.rm=TRUE)

data <- data[, .N, by=.(chr, start, end, read_id, mapq, orientation, insert_size,
                         read_length, flag, mod_cpgs, unmod_cpgs, snp_cpgs,
                         mod_denovos, unmod_denovos)]
.ts(paste0("after dedup: ", nrow(data), " reads"))

.t1 <- proc.time()
data <- parse_positions_vec(data, plot_vbias=args$vbias)
.ts(paste0("parse_positions done: ", nrow(data), " CpG obs (",
           round((proc.time()-.t1)["elapsed"], 1), "s)"))

# V-bias aggregates: compute per-chr pre-aggregated tiles before dropping extra cols
vdata <- NULL
if (args$vbias && nrow(data) > 0) {
  .ts("computing vbias aggregates...")
  vd <- copy(data)
  vd[, Orientation := fifelse(bitwAnd(flag, 96) == 96 | bitwAnd(flag, 144) == 144, "OT", "OB")]
  vd[, abs_pos := start + pos_in_read]
  frag_bounds <- vd[, .(frag_start = min(start), frag_end = max(end)), by = .(read_id, Orientation)]
  vd <- vd[frag_bounds, on = .(read_id, Orientation)]
  vd[Orientation == "OT", pos_in_fragment := abs_pos - frag_start]
  vd[Orientation == "OB", pos_in_fragment := frag_end - abs_pos]
  vdata <- vd[pos_in_fragment >= 0L & pos_in_fragment < insert_size, .(
    n_mod   = sum(type == "M"),
    n_total = .N
  ), by = .(chr, pos_in_fragment, insert_size, Orientation)]
  rm(vd, frag_bounds); gc()
  .ts(paste0("vbias aggregates: ", nrow(vdata), " tiles"))
}

mbias_table <- data[, .(
  mod   = sum(type %in% c("M", "m")),
  unmod = sum(type %in% c("U", "u"))
), by = .(
  chr,
  pos       = pos_in_read,
  read_pair = as.factor(fifelse(bitwAnd(flag, 64) > 0, "First", "Second")),
  strand    = as.factor(fifelse(bitwAnd(flag, 96) == 96 | bitwAnd(flag, 144) == 144, "OT", "OB"))
)]

# Apply 5'->3' flip here so stage 2 does not need per-chr read_len
mbias_table[strand=="OT" & read_pair=="Second" | strand=="OB" & read_pair=="First",
            pos := read_len - pos]

rm(data); gc()
.ts(paste0("mbias_table: ", nrow(mbias_table), " rows"))

# Save table + read_len + vdata (NULL if --vbias not set) so stage 2 can use them
saveRDS(list(mbias_table=mbias_table, read_len=read_len, vdata=vdata), args$output)
.ts(paste0("saved: ", args$output))
