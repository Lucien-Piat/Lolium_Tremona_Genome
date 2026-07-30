#!/usr/bin/env Rscript
# Goal 5: plot Pi, Watterson's theta, Tajima's D per population.
# Usage: Rscript goal5_plot.R <diversity_dir> <pop_palette_tsv>

suppressMessages(library(tidyverse))

args   <- commandArgs(trailingOnly = TRUE)
divdir <- normalizePath(args[1])
pal.fn <- normalizePath(args[2])
setwd(divdir)

pal <- read.table(pal.fn, comment.char = "", col.names = c("col", "pop"))
col.pop <- setNames(pal$col, pal$pop)
pop.levels <- pal$pop

pub_theme <- theme_test(base_size = 22) +
  theme(
    axis.title   = element_text(size = 26),
    axis.text    = element_text(size = 22, colour = "black"),
    plot.title   = element_text(size = 28, face = "bold"),
    axis.ticks   = element_line(linewidth = 0.8)
  )
plot_stat <- function(df, value_col, xlab_txt, file, vline = NULL, val_lim = NULL) {
  df <- df |>
    filter(no_sites >= 5000) |>
    mutate(pop = factor(pop, levels = rev(pop.levels)))
  p <- ggplot(df, aes(pop, .data[[value_col]], fill = pop)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, linewidth = 0.7) +
    scale_fill_manual(values = col.pop, guide = "none") +
    coord_flip(ylim = val_lim) + 
    xlab("Population") + ylab(xlab_txt) +
    pub_theme
  if (!is.null(vline)) p <- p + geom_hline(yintercept = vline,
                                           linetype = "dashed", linewidth = 0.8)
  ggsave(file, p, width = 9, height = 6.5)
}

# Pi (max 0.035)
pi_tab <- read.table("pixy_pi.txt", header = TRUE, stringsAsFactors = FALSE)
plot_stat(pi_tab, "avg_pi",
          expression(Nucleotide~diversity~(pi)), "pi_by_pop.pdf",
          val_lim = c(0, 0.035))

# theta (max 0.045)
wt_tab <- read.table("pixy_watterson_theta.txt", header = TRUE, stringsAsFactors = FALSE)
plot_stat(wt_tab, "avg_watterson_theta",
          expression(Watterson*"'"*s~theta[W]), "watterson_by_pop.pdf",
          val_lim = c(0, 0.045))

# D
td_tab <- read.table("pixy_tajima_d.txt", header = TRUE, stringsAsFactors = FALSE) |>
  filter(is.finite(tajima_d), tajima_d > -5, tajima_d < 5)
plot_stat(td_tab, "tajima_d", expression(Tajima*"'"*s~italic(D)),
          "tajimaD_by_pop.pdf", vline = 0)

summ_pi <- pi_tab |> filter(no_sites >= 5000) |> group_by(pop) |>
  summarise(pi = sum(avg_pi*no_sites, na.rm=TRUE)/sum(no_sites, na.rm=TRUE), .groups="drop")
summ_wt <- wt_tab |> filter(no_sites >= 5000) |> group_by(pop) |>
  summarise(theta_w = sum(avg_watterson_theta*no_sites, na.rm=TRUE)/sum(no_sites, na.rm=TRUE), .groups="drop")
summ_td <- td_tab |> group_by(pop) |>
  summarise(tajimaD = median(tajima_d, na.rm=TRUE), .groups="drop")

summary <- summ_pi |> left_join(summ_wt, by="pop") |> left_join(summ_td, by="pop")
summary$pop <- factor(summary$pop, levels = pop.levels)
summary <- summary[order(summary$pop), ]
write.csv(summary, "diversity_summary.csv", row.names = FALSE)
print(summary)
