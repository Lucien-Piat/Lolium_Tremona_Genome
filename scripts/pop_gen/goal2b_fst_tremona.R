#!/usr/bin/env Rscript
# Goal 2b: FST de Tremona (groupe) contre chaque accession CH individuelle,
# Usage: goal2b_fst_tremona.R <vcf> <pop> <outdir> <chr_prefix>

suppressMessages({
  library(SNPRelate)
  library(tidyverse)
})

args    <- commandArgs(trailingOnly = TRUE)
vcf.fn  <- args[1]
pop.fn  <- args[2]
outdir  <- args[3]
chr.pre <- if (length(args) >= 4) args[4] else ""

setwd(outdir)

pub_theme <- theme_test(base_size = 22) +
  theme(axis.title = element_text(size = 26),
        axis.text  = element_text(size = 20, colour = "black"))

snpgdsVCF2GDS(vcf.fn, "fst_trem.gds", method = "biallelic.only",
              ignore.chr.prefix = chr.pre)
genofile <- snpgdsOpen("fst_trem.gds")
samp.id  <- read.gdsn(index.gdsn(genofile, "sample.id"))

pop <- read.table(pop.fn, col.names = c("sample.id", "population"))
pop.vec <- setNames(pop$population[match(samp.id, pop$sample.id)], samp.id)

trem_ids <- names(pop.vec)[pop.vec == "TREM"]
ch_ids   <- names(pop.vec)[pop.vec == "CH"]
cat("Tremona:", length(trem_ids), "ind ; CH:", length(ch_ids), "accessions\n")

fst_pair <- function(ids_a, ids_b, label_a, label_b) {
  keep <- c(ids_a, ids_b)
  grp  <- factor(c(rep(label_a, length(ids_a)), rep(label_b, length(ids_b))))
  f <- snpgdsFst(genofile, sample.id = keep, population = grp,
                 method = "W&C84", autosome.only = FALSE, verbose = FALSE)
  f$Fst
}

trem_vs_ch <- sapply(ch_ids, function(ch)
  fst_pair(trem_ids, ch, "TREM", "CH_ind"))
df_trem <- data.frame(comparison = "Tremona vs CH accession",
                      pair = paste0("TREM-", ch_ids),
                      fst = trem_vs_ch)

ch_pairs <- combn(ch_ids, 2, simplify = FALSE)
ch_vs_ch <- sapply(ch_pairs, function(p) fst_pair(p[1], p[2], "a", "b"))
df_ch <- data.frame(comparison = "CH vs CH accession",
                    pair = sapply(ch_pairs, function(p) paste(p, collapse="-")),
                    fst = ch_vs_ch)

closefn.gds(genofile)

all_df <- rbind(df_trem, df_ch)
write.csv(all_df, "fst_tremona_vs_ch.csv", row.names = FALSE)

summ <- all_df |> group_by(comparison) |>
  summarise(median_fst = median(fst, na.rm=TRUE),
            mean_fst   = mean(fst, na.rm=TRUE),
            min_fst    = min(fst, na.rm=TRUE),
            max_fst    = max(fst, na.rm=TRUE),
            n = n(), .groups="drop")
print(summ)
write.csv(summ, "fst_tremona_vs_ch_summary.csv", row.names = FALSE)

p <- ggplot(all_df, aes(comparison, fst, fill = comparison)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.5) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
  scale_fill_manual(values = c("Tremona vs CH accession" = "#e9162d",
                               "CH vs CH accession"      = "#8f2be7"),
                    guide = "none") +
  xlab(NULL) + ylab(expression(F[ST])) +
  pub_theme +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("fst_tremona_vs_ch.pdf", p, width = 9, height = 7)

cat("Done. -> fst_tremona_vs_ch.csv, _summary.csv, .pdf\n")