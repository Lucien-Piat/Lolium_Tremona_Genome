#!/usr/bin/env Rscript
# Verification de la normalite pour les tests t apparies de 05_diversity_plot.r.
# Usage: Rscript 05b_normality_qq.R <diversity_dir> <pop_palette_tsv>

suppressMessages(library(tidyverse))

args   <- commandArgs(trailingOnly = TRUE)
divdir <- normalizePath(args[1])
pal.fn <- normalizePath(args[2])
setwd(divdir)
REF <- "TREM"

pal        <- read.table(pal.fn, comment.char = "", col.names = c("col", "pop"))
pop.levels <- pal$pop

# agregation identique a celle de 05_diversity_plot.r
chrom_means <- function(df, value_col, weighted) {
  df |>
    filter(no_sites >= 5000, is.finite(.data[[value_col]])) |>
    group_by(pop, chromosome) |>
    summarise(v = if (weighted) sum(.data[[value_col]] * no_sites) / sum(no_sites)
                  else mean(.data[[value_col]]),
              .groups = "drop") |>
    pivot_wider(names_from = pop, values_from = v)
}

pi_tab <- read.table("pixy_pi.txt", header = TRUE, stringsAsFactors = FALSE)
td_tab <- read.table("pixy_tajima_d.txt", header = TRUE, stringsAsFactors = FALSE) |>
  filter(no_sites >= 5000, is.finite(tajima_d), abs(tajima_d) < 5)

agg_pi <- chrom_means(pi_tab, "avg_pi", weighted = TRUE)
agg_td <- chrom_means(td_tab, "tajima_d", weighted = FALSE)
others <- setdiff(names(agg_pi), c("chromosome", REF))

SET_PI <- "pi difference"
SET_D  <- "Tajima's D"

dat <- bind_rows(
  map_dfr(others, ~tibble(set = SET_PI, pop = .x, value = agg_pi[[REF]] - agg_pi[[.x]])),
  map_dfr(pop.levels, ~tibble(set = SET_D, pop = .x, value = agg_td[[.x]]))
) |>
  mutate(set = factor(set, levels = c(SET_PI, SET_D)),
         pop = factor(pop, levels = pop.levels))

# ---- Shapiro-Wilk, p brut ------------------------------------------------
sw <- dat |>
  group_by(set, pop) |>
  summarise(n = n(), W = shapiro.test(value)$statistic,
            p = shapiro.test(value)$p.value, .groups = "drop")

# une ligne par population, une colonne par test
wide <- sw |> select(set, pop, p) |>
  pivot_wider(names_from = set, values_from = p) |>
  rename(P_pi_difference = !!SET_PI, P_D_chrom_means = !!SET_D) |>
  arrange(pop)

cat("\n== Shapiro-Wilk (p brut, sans correction) : une ligne par population ==\n")
print(as.data.frame(wide |> mutate(across(where(is.numeric), ~round(.x, 3)))),
      row.names = FALSE)
cat(sprintf("\n%d tests ; p minimal = %.3f ; nombre de p < 0.05 : %d\n",
            nrow(sw), min(sw$p), sum(sw$p < 0.05)))

write.table(wide, "shapiro_tests.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE, na = "")

# ---- QQ-plots : 2 lignes (les 2 tests) x 8 colonnes (les populations) -----
lab <- sw |> mutate(txt = sprintf("P = %.3f", p),
                    col = ifelse(p < 0.05, "#9E1B1B", "grey25"),
                    face = ifelse(p < 0.05, "bold", "plain"))

p <- ggplot(dat, aes(sample = value)) +
  stat_qq_line(colour = "grey55", linewidth = 0.5) +
  stat_qq(size = 2, colour = "black") +
  facet_grid(set ~ pop, scales = "free_y", switch = "y") +
  geom_text(data = lab, inherit.aes = FALSE,
            aes(x = -Inf, y = Inf, label = txt, fontface = face, colour = col),
            hjust = -0.12, vjust = 1.5, size = 3.6) +
  scale_colour_identity() +
  scale_x_continuous(breaks = c(-1, 0, 1)) +
  labs(x = "Theoretical quantiles", y = NULL) +
  theme_bw(base_size = 13) +
  theme(strip.background = element_blank(),
        strip.placement  = "outside",
        panel.grid.minor = element_blank())

ggsave("qq_normality.pdf", p, width = 15, height = 4.8)

version groupée
res <- dat |> group_by(set, pop) |> mutate(r = value - mean(value)) |> ungroup()

sw_pool <- res |> group_by(set) |>
  summarise(n = n(), W = shapiro.test(r)$statistic,
            p = shapiro.test(r)$p.value, .groups = "drop")

cat("\n== Shapiro-Wilk sur les residus centres, populations regroupees ==\n")
print(as.data.frame(sw_pool |> mutate(W = round(W, 3), p = round(p, 4))),
      row.names = FALSE)
write.csv(sw_pool, "shapiro_pooled.csv", row.names = FALSE)

labp <- sw_pool |> mutate(txt  = sprintf("P = %.3f", p),
                          col  = ifelse(p < 0.05, "#9E1B1B", "grey25"),
                          face = ifelse(p < 0.05, "bold", "plain"))

pp <- ggplot(res, aes(sample = r)) +
  stat_qq_line(colour = "grey55", linewidth = 0.5) +
  stat_qq(size = 1.8, colour = "black") +
  facet_wrap(~ set, scales = "free", nrow = 1) +
  geom_text(data = labp, inherit.aes = FALSE,
            aes(x = -Inf, y = Inf, label = txt, fontface = face, colour = col),
            hjust = -0.15, vjust = 1.6, size = 4.2) +
  scale_colour_identity() +
  labs(x = "Theoretical quantiles", y = "Centred residuals") +
  theme_bw(base_size = 13) +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank())

ggsave("qq_normality_pooled.pdf", pp, width = 12, height = 3.4)
cat("Ecrit : qq_normality.pdf, qq_normality_pooled.pdf,",
    "shapiro_tests.tsv, shapiro_pooled.csv\n")
