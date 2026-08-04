#!/usr/bin/env Rscript
# Diversity: pi and Tajima's D per population, with chromosome-level tests.
# Usage: Rscript 05_diversity_plot.R <diversity_dir> <pop_palette_tsv>

suppressMessages(library(tidyverse))

args   <- commandArgs(trailingOnly = TRUE)
divdir <- normalizePath(args[1])
pal.fn <- normalizePath(args[2])
setwd(divdir)

pal        <- read.table(pal.fn, comment.char = "", col.names = c("col", "pop"))
col.pop    <- setNames(pal$col, pal$pop)
pop.levels <- pal$pop
REF        <- "TREM"
TOST_MARGIN <- 0.20   # equivalence margin for pi, as a fraction of the grand mean

pub_theme <- theme_test(base_size = 22) +
  theme(axis.title    = element_text(size = 26),
        axis.text     = element_text(size = 22, colour = "black"),
        plot.subtitle = element_text(size = 20),
        axis.ticks    = element_line(linewidth = 0.8))

fmt <- function(p) {
  if (p >= 0.05) "P > 0.05"
  else if (p >= 0.001) sprintf("P = %.3f", p) else "P < 0.001"
}

# Under TOST a small P establishes equivalence
equiv_lab <- function(p) if (p < 0.05) "equiv." else "not equiv."

# Neighbouring 20 kb windows are linked and are not independent replicates.
# Windows are therefore aggregated into one value per chromosome, and the seven
# chromosomes are used as the paired units, chromosome means average thousands of windows
# and are approximately normal, so a t-test is used.
chrom_means <- function(df, value_col, weighted) {
  df |>
    filter(no_sites >= 5000, is.finite(.data[[value_col]])) |>
    group_by(pop, chromosome) |>
    summarise(v = if (weighted) sum(.data[[value_col]] * no_sites) / sum(no_sites)
                  else mean(.data[[value_col]]),
              .groups = "drop") |>
    pivot_wider(names_from = pop, values_from = v)
}

# ---- effect sizes: per-population mean and 95% CI over the 7 chromosomes ----
pop_summary <- function(agg, label) {
  pops <- setdiff(names(agg), "chromosome")
  out <- map_dfr(pops, function(p) {
    tt <- t.test(agg[[p]])
    tibble(pop = p, mean = tt$estimate,
           ci_lo = tt$conf.int[1], ci_hi = tt$conf.int[2],
           p_vs_0 = tt$p.value)
  })
  cat(sprintf("\n== %s : per-population mean, 95%% CI, and test vs 0 ==\n", label))
  print(as.data.frame(out |> mutate(across(where(is.numeric), ~signif(.x, 3)))),
        row.names = FALSE)
  out
}

# ---- pi : equivalence (TOST), two one-sided tests against +/- margin ----
pi_tab <- read.table("pixy_pi.txt", header = TRUE, stringsAsFactors = FALSE)
agg_pi <- chrom_means(pi_tab, "avg_pi", weighted = TRUE)
sum_pi <- pop_summary(agg_pi, "Nucleotide diversity (pi)")

delta  <- TOST_MARGIN * mean(unlist(agg_pi[, setdiff(names(agg_pi), "chromosome")]))
others <- setdiff(names(agg_pi), c("chromosome", REF))

tost_p <- sapply(others, function(p) {
  d <- agg_pi[[REF]] - agg_pi[[p]]
  max(t.test(d, mu = -delta, alternative = "greater")$p.value,
      t.test(d, mu =  delta, alternative = "less")$p.value)
})
tost_padj <- p.adjust(tost_p, "holm")

cat(sprintf("\n== pi : equivalence of %s vs each population (TOST, margin = %.2g) ==\n",
            REF, delta))
print(data.frame(pop = others, p_holm = signif(unname(tost_padj), 3),
                 equivalent = unname(tost_padj) < 0.05), row.names = FALSE)

# ord decroissant => crochet "not equiv." d'abord, puis "equiv."
pw_pi <- data.frame(pop = others,
                    label = unname(sapply(tost_padj, equiv_lab)),
                    sig   = unname(tost_padj) < 0.05,
                    ord   = unname(tost_padj))

# Tajima's D : is it closer to zero in TREM ? 
td_tab <- read.table("pixy_tajima_d.txt", header = TRUE, stringsAsFactors = FALSE)
n_all  <- sum(td_tab$no_sites >= 5000)
td_tab <- td_tab |> filter(no_sites >= 5000, is.finite(tajima_d), abs(tajima_d) < 5)
cat(sprintf("\nTajima's D: %d of %d windows discarded as non-finite or |D| >= 5\n",
            n_all - nrow(td_tab), n_all))

agg_td <- chrom_means(td_tab, "tajima_d", weighted = FALSE)
sum_td <- pop_summary(agg_td, "Tajima's D")

# one-sided paired tests on the distance of each chromosome mean from zero
abs_td  <- agg_td |> mutate(across(-chromosome, abs))
td_p    <- sapply(others, function(p)
             t.test(abs_td[[REF]], abs_td[[p]], paired = TRUE,
                    alternative = "less")$p.value)
td_padj <- p.adjust(td_p, "holm")

cat(sprintf("\n== Tajima's D : |D| in %s vs each population (paired, one-sided) ==\n", REF))
print(data.frame(pop = others, p_holm = signif(unname(td_padj), 3),
                 closer_to_0 = unname(td_padj) < 0.05), row.names = FALSE)

pw_td <- data.frame(pop = others,
                    label = unname(sapply(td_padj, fmt)),
                    sig   = unname(td_padj) < 0.05)

plot_stat <- function(df, value_col, xlab_txt, file, vline = NULL,
                      val_lim = NULL, annot = NULL, pw = NULL, ref_lab = NULL,
                      width = 9.5) {
  df <- df |> filter(no_sites >= 5000) |>
    mutate(pop = factor(pop, levels = rev(pop.levels)))

  lim <- val_lim; ann <- NULL; refdf <- NULL
  if (!is.null(pw) || !is.null(ref_lab)) {
    w <- df |> group_by(pop) |>
      summarise(hi = quantile(.data[[value_col]], .75, na.rm = TRUE) +
                     1.5 * IQR(.data[[value_col]], na.rm = TRUE), .groups = "drop")
    lo   <- if (!is.null(val_lim)) val_lim[1] else min(df[[value_col]], na.rm = TRUE)
    dhi  <- max(w$hi); span <- dhi - lo
    step <- 0.052 * span; tick <- 0.026 * span
  }
  if (!is.null(pw)) {
    ann  <- pw |> mutate(pop  = factor(pop, levels = levels(df$pop)),
                         ypos = as.integer(pop),
                         col  = ifelse(sig, "black", "grey45"))
    grp  <- ann |> group_by(label, col) |>
      summarise(ymax = max(ypos), ord = max(ord), .groups = "drop") |>
      arrange(desc(ord)) |>
      mutate(x0 = dhi + 0.035 * span + (row_number() - 1) * step)
    ann  <- left_join(ann, select(grp, label, x0), by = "label")
    lim  <- c(lo, max(grp$x0) + 0.045 * span)
  }
  if (!is.null(ref_lab)) {
    refdf <- data.frame(ypos = which(levels(df$pop) == REF),
                        x = dhi + 0.04 * span, label = ref_lab)
    lim <- c(lo, dhi + 0.30 * span)
  }

  p <- ggplot(df, aes(pop, .data[[value_col]], fill = pop)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, linewidth = 0.7) +
    scale_fill_manual(values = col.pop, guide = "none") +
    coord_flip(ylim = lim) + xlab("Population") + ylab(xlab_txt) + pub_theme
  if (!is.null(ann) || !is.null(refdf)) {
    br <- pretty(c(lim[1], max(w$hi)), n = 4)
    p <- p + scale_y_continuous(breaks = br[br <= max(w$hi)])
  }
  if (!is.null(vline)) p <- p + geom_hline(yintercept = vline,
                                           linetype = "dashed", linewidth = 0.8)
  if (!is.null(annot)) p <- p + labs(subtitle = annot)
  if (!is.null(ann)) {
    p <- p +
      geom_segment(data = grp, inherit.aes = FALSE,
                   aes(x = 1, xend = ymax, y = x0, yend = x0, colour = col),
                   linewidth = 0.5) +
      geom_segment(data = grp, inherit.aes = FALSE,
                   aes(x = 1, xend = 1, y = x0 - tick, yend = x0, colour = col),
                   linewidth = 0.5) +
      geom_segment(data = ann, inherit.aes = FALSE,
                   aes(x = ypos, xend = ypos, y = x0 - tick, yend = x0, colour = col),
                   linewidth = 0.5) +
      geom_text(data = grp, inherit.aes = FALSE,
                aes(x = (1 + ymax) / 2, y = x0, label = label, colour = col),
                angle = 90, vjust = -0.5, hjust = 0.5, size = 4.6) +
      scale_colour_identity()
  }
  if (!is.null(refdf)) {
    p <- p + geom_text(data = refdf, inherit.aes = FALSE,
                       aes(x = ypos, y = x, label = label),
                       hjust = 0, size = 5.4, colour = "black")
  }
  ggsave(file, p, width = width, height = 6.5)
}

plot_stat(pi_tab, "avg_pi", expression(Nucleotide~diversity~(pi)),
          "pi_by_pop.pdf", val_lim = c(0, 0.035), pw = pw_pi,
          annot = sprintf("Equivalence of Tremona (TOST, margin = %.2g)", delta))

p_trem_0 <- sum_td$p_vs_0[sum_td$pop == REF]
plot_stat(td_tab, "tajima_d", expression(Tajima*"'"*s~italic(D)),
          "tajimaD_by_pop.pdf", vline = 0,
          ref_lab = fmt(p_trem_0),
          annot = "Tremona: Tajima's D vs. 0")

summary_tab <- sum_pi |> rename_with(~paste0("pi_", .x), -pop) |>
  left_join(sum_td |> rename_with(~paste0("D_", .x), -pop), by = "pop") |>
  mutate(pop = factor(pop, levels = pop.levels)) |> arrange(pop)
write.csv(summary_tab, "diversity_summary.csv", row.names = FALSE)

#  test 1 : equivalence de pi (TOST) 
tost_tab <- tibble(pop = others,
                   mean_pi = sapply(others, function(p) mean(agg_pi[[p]])),
                   mean_pi_TREM = mean(agg_pi[[REF]]),
                   diff = mean_pi_TREM - mean_pi,
                   margin = delta,
                   p_raw = unname(tost_p),
                   p_holm = unname(tost_padj),
                   equivalent = unname(tost_padj) < 0.05,
                   label = unname(sapply(tost_padj, equiv_lab))) |>
  mutate(pop = factor(pop, levels = pop.levels)) |> arrange(pop)
write.table(tost_tab, "pi_tost_tests.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)

#  test 2 : Tajima's D (moyenne vs 0, et |D| vs Tremona) 
td_tab_out <- sum_td |>
  rename(mean_D = mean, ci_lo = ci_lo, ci_hi = ci_hi, p_vs_0 = p_vs_0) |>
  left_join(tibble(pop = others,
                   absD_p_raw = unname(td_p),
                   absD_p_holm = unname(td_padj),
                   closer_to_0_than_TREM = unname(td_padj) < 0.05),
            by = "pop") |>
  mutate(pop = factor(pop, levels = pop.levels)) |> arrange(pop)
write.table(td_tab_out, "tajimaD_tests.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)

cat("\nEcrit : pi_tost_tests.tsv, tajimaD_tests.tsv, diversity_summary.csv\n")
print(as.data.frame(summary_tab))