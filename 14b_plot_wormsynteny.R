#!/usr/bin/env Rscript
# Step 14: WormSynteny-style riparian plot of the HelSus SLC4A8-AQP5 interval, drawn the way
# the WormSynteny app draws it (gggenomes: geom_seq + geom_gene filled by orthogroup +
# geom_link), from the tables of 14a_wormsynteny_tables.py.
#
#   * rows    = the longest stitched fragment of each genome (WormSynteny's main fragment),
#               in species-tree order; with 14a --extras, the other fragments are added
#               as "(filtered)" rows at the bottom, as in WormSynteny.
#   * genes   = each genome's own TOGA2 annotation, filled by "orthogroup": the step-5 name
#               for exendin-family genes (fixed project colours), else the HelSus gene
#               symbol TOGA2 projected from (WormSynteny's orthogroup palette).
#               Hollow = TOGA2 status not intact.
#   * links   = second-wave halLiftover alignments between neighbouring rows.
#   * sync()  flips each genome to match its links, as WormSynteny's plots are oriented.
#
# Usage
#   Rscript 14b_plot_wormsynteny.R --in wormsynteny --out wormsynteny_region \
#       [--width 10] [--height 7] [--labels named|all|none] [--legend named|all]

suppressPackageStartupMessages({
  library(gggenomes)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

a <- commandArgs(trailingOnly = TRUE)
opt <- list(`in` = NULL, out = NULL, width = 10, height = 7, labels = "named", legend = "named",
            dpi = 600, link_colour = "grey70")
i <- 1
while (i <= length(a)) {
  k <- gsub("-", "_", sub("^--", "", a[i]))
  if (!k %in% names(opt)) stop(sprintf("unknown option %s", a[i]))
  opt[[k]] <- if (is.numeric(opt[[k]])) as.numeric(a[i + 1]) else a[i + 1]
  i <- i + 2
}
if (is.null(opt$`in`) || is.null(opt$out)) stop("--in <14a prefix> and --out <prefix> are required")

rd <- function(x) read_tsv(paste0(opt$`in`, "_", x, ".tsv"), show_col_types = FALSE,
                           col_types = cols(.default = col_character(), start = col_double(),
                                            end = col_double(), length = col_double(),
                                            start2 = col_double(), end2 = col_double(),
                                            named = col_integer()))
seqs <- rd("seqs")
genes <- rd("genes")
links <- rd("links")
message(sprintf("[info] %d rows, %d genes, %d links", nrow(seqs), nrow(genes), nrow(links)))

# ---- orthogroup colours: project colours for the named exendin array, WormSynteny's palette for the rest ----
ws_palette <- c("#92C5DEFF", "#EEB479FF", "#AFE1AF", "#E88471FF", "#E9E29CFF", "#CF597EFF", "#E78AC3FF",
                "#FFD92FFF", "#009392FF", "#FF9E9D", "#FF3D7FFF", "#EBA07EFF", "#A8554EFF", "#4393C3FF",
                "#39B185FF", "#9CCB86FF", "#A6D854FF", "#803233FF", "#ED3F39FF", "#E46844FF",
                "#DC876CFF", "#FFCC65FF", "#FF9932FF")           # assign_colors() in WormSynteny
named <- genes %>% filter(named == 1) %>% distinct(og, colour)
others <- setdiff(unique(genes$og), named$og)
others <- others[order(-table(genes$og)[others], others)]      # most widespread first
oc <- setNames(rep_len(ws_palette, length(others)), others)
if (length(others) > length(ws_palette))
  message(sprintf("[warn] %d other orthogroups for %d palette colours: colours repeat", length(others), length(ws_palette)))
cols <- c(setNames(named$colour, named$og), oc)
genes <- genes %>% mutate(intact = is.na(status) | status %in% c("I", "FI", "PI", ""),
                          tag = case_when(opt$labels == "all" ~ og,
                                          opt$labels == "named" & named == 1 ~ label,
                                          TRUE ~ ""))
breaks <- if (opt$legend == "all") names(cols) else named$og

# ---- plot ----
p <- gggenomes(genes = genes, seqs = seqs, links = links) |>
  sync() +
  geom_link(fill = opt$link_colour, colour = NA, alpha = 0.6) +
  geom_seq() +
  geom_bin_label(fontface = "italic") +
  geom_gene(aes(fill = og), data = genes(intact), stroke = 0.3) +
  geom_gene(aes(colour = og), data = genes(!intact), fill = "white", stroke = 0.6, show.legend = FALSE) +
  geom_gene_tag(aes(label = tag), data = genes(tag != ""), nudge_y = 0.1, check_overlap = TRUE, size = 2.5) +
  scale_fill_manual("Orthogroup", values = cols, breaks = breaks) +
  scale_colour_manual(values = cols, guide = "none") +
  labs(caption = "Links: halLiftover between neighbouring rows. Hollow genes: TOGA2 not intact.") +
  theme(legend.position = "right")

ggsave(paste0(opt$out, ".pdf"), p, width = opt$width, height = opt$height)
ggsave(paste0(opt$out, ".png"), p, width = opt$width, height = opt$height, dpi = opt$dpi)
message(sprintf("[info] wrote %s.pdf and %s.png", opt$out, opt$out))
