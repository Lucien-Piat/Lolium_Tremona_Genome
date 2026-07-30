#!/usr/bin/env Rscript
# Goal 4 : trace du LD decay. Courbes binnees pour l'affichage,
# Usage: Rscript goal4_ld_plot.R <ld_decay_dir> <pop_palette_tsv>

suppressMessages(library(tidyverse))

args   <- commandArgs(trailingOnly = TRUE)
lddir  <- normalizePath(args[1])
pal.fn <- normalizePath(args[2])
setwd(lddir)

pal <- read.table(pal.fn, comment.char = "", col.names = c("col", "pop"))
col.pop <- setNames(pal$col, pal$pop)

read_ld <- function(p) {
  f <- paste0(p, ".stat.gz")
  d <- read.table(gzfile(f), header = TRUE, comment.char = "",
                  stringsAsFactors = FALSE)
  data.frame(pop = p, dist = d[[1]], r2 = d[[2]])
}

pops <- names(col.pop)
ld <- do.call(rbind, lapply(pops, read_ld))

halfdist <- ld |>
  group_by(pop) |>
  summarise(
    r2_start = r2[which.min(dist)],
    r2_plateau = mean(r2[dist > 250000], na.rm = TRUE),
    r2_half = (r2_start + r2_plateau) / 2,
    half_dist_kb = {
      idx <- which(r2 <= (r2_start + r2_plateau) / 2)[1]
      if (is.na(idx)) NA else dist[idx] / 1000
    },
    .groups = "drop") |>
  arrange(desc(half_dist_kb))
write.csv(halfdist, "ld_halfdist.csv", row.names = FALSE)
print(halfdist)

p3 <- ggplot(halfdist, aes(reorder(pop, half_dist_kb), half_dist_kb, fill = pop)) +
  geom_col() +
  scale_fill_manual(values = col.pop, guide = "none") +
  xlab("Population") + ylab("Demi-distance de decay du LD (kb)") +
  coord_flip() +
  theme_test()
ggsave("ld_halfdist.pdf", p3, width = 7, height = 5)

BIN <- 500
ld_binned <- ld |>
  mutate(dist_bin = floor(dist / BIN) * BIN + BIN / 2) |>
  group_by(pop, dist_bin) |>
  summarise(r2 = mean(r2, na.rm = TRUE), .groups = "drop") |>
  rename(dist = dist_bin)

p1 <- ggplot(ld_binned, aes(dist / 1000, r2, colour = pop)) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = col.pop) +
  xlab("Distance (kb)") + ylab(expression(LD~(r^2))) +
  coord_cartesian(xlim = c(0, 300)) +
  theme_test()
ggsave("ld_decay_all.pdf", p1, width = 9, height = 6)

