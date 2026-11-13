#!/usr/bin/env Rscript
# Step 13: GENESPACE riparian plot of the Cactus alignment, from the blk CSV of 13a_hal_to_blk.py.
#
# The drawing is GENESPACE's own riparian_engine(); this script only picks colours and calls it.
#   * braids: one colour per HelSus chromosome (GENESPACE palette, gs_colors), traced by 13a;
#     blocks that 13a could not trace back to HelSus are --other (grey).
#   * the HelSus SLC4A8-AQP5 interval (13a highlight rows) is drawn last, in --highlight.
#   * --mono: every non-interval braid grey, so only the interval is coloured
#     (the look of plot_riparian(highlightBed = ...)).
# Chromosomes are scaled by physical position (useOrder = FALSE): there are no genes here.
# reorderBySynteny = FALSE, because GENESPACE's reordering needs its own gene-based bed file;
# non-reference contigs are ordered by the number in their name.
#
# Usage
#   Rscript 13b_plot_riparian.R --blk blk_genome.csv --genomes genomes.tsv --out genespace_genome \
#       [--width 10] [--height 7] [--min-chr-len 1e6] [--highlight "#7B3294"] [--mono] \
#       [--inversion-colour "#A6611A"] [--braid-alpha 0.6]
# Writes <out>.pdf and <out>.png.

suppressPackageStartupMessages({
  library(GENESPACE)
  library(data.table)
  library(ggplot2)
})

# ---- arguments: --key value, or --flag ----
a <- commandArgs(trailingOnly = TRUE)
opt <- list(blk = NULL, genomes = NULL, out = NULL, width = 10, height = 7, min_chr_len = 1e5,
            highlight = "#7B3294", other = "grey55", mono = FALSE, inversion_colour = NULL,
            braid_alpha = 0.6, font = 5, dpi = 600)
i <- 1
while (i <= length(a)) {
  k <- gsub("-", "_", sub("^--", "", a[i]))
  if (!k %in% names(opt)) stop(sprintf("unknown option %s", a[i]))
  if (is.logical(opt[[k]])) { opt[[k]] <- TRUE; i <- i + 1; next }
  v <- a[i + 1]
  opt[[k]] <- if (is.numeric(opt[[k]])) as.numeric(v) else v
  i <- i + 2
}
for (k in c("blk", "genomes", "out")) if (is.null(opt[[k]])) stop(sprintf("--%s is required", k))

blk <- fread(opt$blk, colClasses = list(character = c("chr1", "chr2", "refChr")))
gids <- fread(opt$genomes)$code
gids <- gids[gids %in% c(blk$genome1, blk$genome2)]
ref <- blk$refGenome[1]
message(sprintf("[info] %d blocks, %d interval rows, %d genomes: %s",
                sum(blk$highlight == 0), sum(blk$highlight == 1), length(gids),
                paste(gids, collapse = ", ")))

# ---- drop blocks on contigs riparian_engine will not draw (it leaves them with NA positions) ----
ext <- rbind(blk[, .(genome = genome1, chr = chr1, e = pmax(startBp1, endBp1))],
             blk[, .(genome = genome2, chr = chr2, e = pmax(startBp2, endBp2))])[
               , .(len = max(e)), by = .(genome, chr)]
keep <- ext[len >= opt$min_chr_len, paste(genome, chr)]
n0 <- nrow(blk)
blk <- blk[paste(genome1, chr1) %in% keep & paste(genome2, chr2) %in% keep]
if (nrow(blk) < n0)
  message(sprintf("[info] --min-chr-len %g: %d contig(s), %d block(s) left out",
                  opt$min_chr_len, sum(ext$len < opt$min_chr_len), n0 - nrow(blk)))

# ---- colours ----
blk[is.na(refChr), refChr := ""]
rc <- unique(blk$refChr[blk$refChr != ""])
rc <- rc[order(suppressWarnings(as.numeric(gsub("\\D+", "", rc))), rc)]
pal <- setNames(gs_colors(length(rc)), rc)
blk[, color := ifelse(refChr == "", opt$other, pal[refChr])]
if (opt$mono) blk[, color := opt$other]
blk[highlight == 1, color := opt$highlight]
bg <- blk[highlight == 0, blkID]
blk[, refChr := NULL]          # riparian_engine would append it to blkID, breaking theseBlocksFirst

# ---- draw: riparian_engine sizes the braids from the open device, so open it first ----
pdf(paste0(opt$out, ".pdf"), width = opt$width, height = opt$height)
rp <- riparian_engine(
  blk = blk, bed = NULL,
  refGenome = ref, genomeIDs = gids,
  theseBlocksFirst = bg,
  useOrder = FALSE, reorderBySynteny = FALSE,
  minChrLen2plot = opt$min_chr_len,
  braidAlpha = opt$braid_alpha,
  inversionColor = opt$inversion_colour,
  chrFill = "white", chrBorderCol = "grey30", chrBorderLwd = 0.2,
  chrLabFontSize = opt$font,
  chrLabFun = function(x) sub("^(chr|chromosome|scaffold|super)_?", "", x, ignore.case = TRUE),
  xlabel = "Chromosomes scaled by physical position (Cactus / halSynteny blocks)")
print(rp$ggplotObj)
invisible(dev.off())
ggsave(paste0(opt$out, ".png"), rp$ggplotObj, width = opt$width, height = opt$height, dpi = opt$dpi)
fwrite(rp$sourceData$chromosomes, paste0(opt$out, "_chromosomes.tsv"), sep = "\t")
message(sprintf("[info] wrote %s.pdf, %s.png, %s_chromosomes.tsv", opt$out, opt$out, opt$out))
