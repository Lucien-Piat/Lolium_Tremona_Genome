#!/usr/bin/env Rscript
# Manhattan plots: Pi (Tremona + all pops) and FST, fully rasterised,
# auto-zoom on FST outlier regions with rolling mean.
# Usage: Rscript goal8_manhattan.R <div_dir> <fst_dir> <out_dir>

suppressMessages({
  library(tidyverse); library(zoo); library(ggrastr); library(patchwork)
})

args   <- commandArgs(trailingOnly = TRUE)
divdir <- normalizePath(args[1])
fstdir <- normalizePath(args[2])
outdir <- normalizePath(args[3])

MIN_SNPS   <- 50
ROLL_N     <- 15        # 10 x 20 kb = 200 kb
ROLL_ZOOM  <- 8         # rolling window inside zoom panels (100 kb)
PI_QUANT   <- 0.01
RAST_DPI   <- 300
ZOOM_FLANK <- 7000000
MERGE_GAP  <- 2000000
PT_SIZE    <- 0.3       # solid small points
PT_OUT     <- 0.9       # outlier points
LINE_W     <- 0.4       # rolling mean line width

chr_cols <- c(chr1="#1f77b4", chr2="#ff7f0e", chr3="#2ca02c",
              chr4="#9467bd", chr5="#8c564b", chr6="#e377c2", chr7="#17becf")

minimal_theme <- theme_classic(base_size = 22) +
  theme(axis.title = element_text(size = 28),
        axis.text  = element_text(size = 20, colour = "black"),
        strip.text = element_text(size = 22),
        axis.line  = element_line(colour = "black", linewidth = 0.6),
        panel.border = element_blank(), panel.grid = element_blank(),
        legend.position = "none")

add_cumpos <- function(df) {
  chr_len <- df |> group_by(chromosome) |>
    summarise(len = max(window_pos_2, na.rm=TRUE), .groups="drop") |>
    arrange(chromosome) |> mutate(offset = lag(cumsum(len), default = 0))
  df |> left_join(chr_len, by="chromosome") |>
    mutate(cumpos = window_pos_1 + offset) |> arrange(chromosome, window_pos_1)
}
add_roll <- function(df, col, n = ROLL_N) {
  df |> group_by(chromosome) |>
    mutate(roll = rollmean(.data[[col]], k=n, fill=NA, align="center")) |> ungroup()
}

pi_tab <- read.table(file.path(divdir,"pixy_pi.txt"), header=TRUE) |>
  filter(no_sites >= 5000)
pi_trem <- pi_tab |> filter(pop=="TREM") |> add_cumpos() |> add_roll("avg_pi")

lo <- quantile(pi_trem$avg_pi, PI_QUANT,   na.rm=TRUE)
hi <- quantile(pi_trem$avg_pi, 1-PI_QUANT, na.rm=TRUE)
pi_trem <- pi_trem |> mutate(status = case_when(
  avg_pi <= lo ~ "low", avg_pi >= hi ~ "high", TRUE ~ "normal"))

chr_mid <- pi_trem |> group_by(chromosome) |>
  summarise(mid = mean(range(cumpos)), .groups="drop")

p_pi <- ggplot(pi_trem, aes(cumpos, avg_pi)) +
  rasterise(geom_point(aes(colour = chromosome), size = PT_SIZE), dpi = RAST_DPI) +
  rasterise(geom_line(aes(y = roll), colour = "black", linewidth = LINE_W),
            dpi = RAST_DPI) +
  scale_colour_manual(values = chr_cols) +
  scale_x_continuous(breaks = chr_mid$mid, labels = chr_mid$chromosome,
                     expand = c(0.01, 0)) +
  xlab("Chromosome") + ylab(expression(pi~"(Tremona, 20 kb)")) +
  minimal_theme
ggsave(file.path(outdir,"manhattan_pi_tremona.pdf"), p_pi, width = 22, height = 5)

write.csv(pi_trem |> filter(status != "normal") |>
            select(chromosome, window_pos_1, window_pos_2, avg_pi, no_sites, status),
          file.path(outdir,"pi_outlier_windows.csv"), row.names=FALSE)

pi_all <- pi_tab |> add_cumpos() |> add_roll("avg_pi")

p_pi_all <- ggplot(pi_all, aes(cumpos, avg_pi)) +
  rasterise(geom_point(aes(colour = chromosome), size = PT_SIZE), dpi = RAST_DPI) +
  rasterise(geom_line(aes(y = roll), colour = "black", linewidth = LINE_W),
            dpi = RAST_DPI) +
  scale_colour_manual(values = chr_cols) +
  scale_x_continuous(breaks = chr_mid$mid, labels = chr_mid$chromosome,
                     expand = c(0.01, 0)) +
  facet_wrap(~ pop, ncol = 1, scales = "free_y") +
  xlab("Chromosome") + ylab(expression(pi~"(20 kb)")) +
  minimal_theme
ggsave(file.path(outdir,"manhattan_pi_all_pops.pdf"), p_pi_all,
       width = 22, height = 20)

fst_tab <- read.table(file.path(fstdir,"pixy_fst.txt"), header=TRUE)
fst_col <- "avg_wc_fst"; snp_col <- "no_snps"

fst_tc <- fst_tab |>
  filter((pop1=="TREM" & pop2=="CH") | (pop1=="CH" & pop2=="TREM")) |>
  filter(.data[[snp_col]] >= MIN_SNPS, .data[[fst_col]] >= 0) |>
  add_cumpos() |> add_roll(fst_col)

thr <- quantile(fst_tc[[fst_col]], 0.999, na.rm=TRUE)
cat(sprintf("Fenetres : %d ; seuil FST (99.9e pct) : %.4f\n", nrow(fst_tc), thr))
fst_out <- fst_tc |> filter(.data[[fst_col]] >= thr)

FST_YLIM <- range(fst_tc[[fst_col]], na.rm = TRUE)
cat(sprintf("Echelle y FST commune : %.3f - %.3f\n", FST_YLIM[1], FST_YLIM[2]))

p_fst <- ggplot(fst_tc, aes(cumpos, .data[[fst_col]])) +
  rasterise(geom_point(aes(colour = chromosome), size = PT_SIZE), dpi = RAST_DPI) +
  rasterise(geom_point(data = fst_out, colour = "#cc0000", size = PT_OUT),
            dpi = RAST_DPI) +
  rasterise(geom_line(aes(y = roll), colour = "black", linewidth = LINE_W),
            dpi = RAST_DPI) +
  geom_hline(yintercept = thr, linetype = "dashed",
             colour = "#cc0000", linewidth = 0.6) +
  scale_colour_manual(values = chr_cols) +
  scale_x_continuous(breaks = chr_mid$mid, labels = chr_mid$chromosome,
                     expand = c(0.01, 0)) +
  coord_cartesian(ylim = FST_YLIM) +
  xlab("Chromosome") + ylab(expression(F[ST]~"(Tremona vs CH)")) +
  minimal_theme
ggsave(file.path(outdir,"manhattan_fst_trem_ch.pdf"), p_fst, width = 22, height = 5)

write.csv(fst_out |> select(chromosome, window_pos_1, window_pos_2,
                            all_of(c(fst_col, snp_col))) |>
            arrange(desc(.data[[fst_col]])),
          file.path(outdir,"fst_outlier_windows.csv"), row.names=FALSE)

regions <- fst_out |>
  arrange(chromosome, window_pos_1) |>
  group_by(chromosome) |>
  mutate(prev_end = lag(window_pos_2),
         gap = if_else(is.na(prev_end), Inf, as.numeric(window_pos_1 - prev_end)),
         region_id = cumsum(gap > MERGE_GAP)) |>
  group_by(chromosome, region_id) |>
  summarise(start = min(window_pos_1), end = max(window_pos_2),
            n_win = n(), max_fst = max(.data[[fst_col]]), .groups="drop")
cat(sprintf("Regions outliers : %d\n", nrow(regions)))
print(regions)
write.csv(regions, file.path(outdir,"fst_outlier_regions.csv"), row.names=FALSE)

make_zoom <- function(i) {
  r  <- regions[i, ]
  x1 <- max(0, r$start - ZOOM_FLANK); x2 <- r$end + ZOOM_FLANK
  d  <- fst_tc |> filter(chromosome == r$chromosome,
                         window_pos_1 >= x1, window_pos_2 <= x2) |>
    arrange(window_pos_1) |>
    mutate(roll_z = rollmean(.data[[fst_col]], k = ROLL_ZOOM,
                             fill = NA, align = "center"))
  ggplot(d, aes(window_pos_1/1e6, .data[[fst_col]])) +
    rasterise(geom_point(size = 1.2, colour = chr_cols[[r$chromosome]]),
              dpi = RAST_DPI) +
    rasterise(geom_point(data = filter(d, .data[[fst_col]] >= thr),
                         colour = "#cc0000", size = 2.2), dpi = RAST_DPI) +
    rasterise(geom_line(aes(y = roll_z), colour = "black", linewidth = LINE_W),
              dpi = RAST_DPI) +
    geom_hline(yintercept = thr, linetype = "dashed",
               colour = "#cc0000", linewidth = 0.5) +
    coord_cartesian(ylim = FST_YLIM) +
    labs(title = sprintf("%s : %.2f-%.2f Mb (n=%d, max=%.2f)",
                         r$chromosome, r$start/1e6, r$end/1e6, r$n_win, r$max_fst),
         x = "Position (Mb)", y = expression(F[ST])) +
    theme_classic(base_size = 16) +
    theme(plot.title  = element_text(size = 15, face = "bold"),
          axis.title  = element_text(size = 17),
          axis.text   = element_text(size = 14, colour = "black"),
          panel.grid  = element_blank())
}

zooms  <- lapply(seq_len(nrow(regions)), make_zoom)
ncol_z <- 2                       # 2 colonnes => panneaux plus larges
nrow_z <- ceiling(length(zooms) / ncol_z)
panel  <- wrap_plots(zooms, ncol = ncol_z)
ggsave(file.path(outdir,"fst_outlier_zooms.pdf"), panel,
       width = 10 * ncol_z, height = 4.5 * nrow_z, limitsize = FALSE)

cat("Done.\n")