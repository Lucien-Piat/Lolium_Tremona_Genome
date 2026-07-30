#!/usr/bin/env Rscript
# Goal 6: Isolation by distance among Swiss 
# Usage: Rscript 06_ibd.R <dxy_dir> <coords_csv>

suppressMessages({
  library(tidyverse)
  library(geosphere)
  library(vegan)
})

args     <- commandArgs(trailingOnly = TRUE)
dxydir   <- normalizePath(args[1])
coord.fn <- normalizePath(args[2])
setwd(dxydir)

# mean Dxy per pair
dxy <- read.table("ch_dxy.txt", header = TRUE, stringsAsFactors = FALSE)
pair_dxy <- dxy |>
  filter(no_sites >= 5000) |>
  group_by(pop1, pop2) |>
  summarise(dxy = sum(avg_dxy * no_sites, na.rm = TRUE) / sum(no_sites, na.rm = TRUE),
            .groups = "drop")

coords   <- read.csv(coord.fn, stringsAsFactors = FALSE)
dxy_accs <- unique(c(dxy$pop1, dxy$pop2))
coords   <- coords[coords$accession %in% dxy_accs, ]
accs     <- coords$accession
n        <- length(accs)

gen_mat <- matrix(0, n, n, dimnames = list(accs, accs))
geo_mat <- matrix(0, n, n, dimnames = list(accs, accs))

for (i in 1:(n-1)) for (j in (i+1):n) {
  a <- accs[i]; b <- accs[j]
  d <- pair_dxy$dxy[(pair_dxy$pop1==a & pair_dxy$pop2==b) |
                    (pair_dxy$pop1==b & pair_dxy$pop2==a)]
  gen_mat[i,j] <- gen_mat[j,i] <- ifelse(length(d)==1, d, NA)
  km <- distHaversine(c(coords$lon[i], coords$lat[i]),
                      c(coords$lon[j], coords$lat[j])) / 1000
  geo_mat[i,j] <- geo_mat[j,i] <- km
}

gen_d <- as.dist(gen_mat)
geo_d <- as.dist(geo_mat)
mant  <- mantel(gen_d, geo_d, permutations = 9999, na.rm = TRUE)
print(mant)

plot_df <- data.frame(geo_km = as.vector(geo_d), dxy = as.vector(gen_d))

pub_theme <- theme_test(base_size = 22) +
  theme(
    axis.title = element_text(size = 26),
    axis.text  = element_text(size = 22, colour = "black"),
    plot.title = element_text(size = 24, face = "bold")
  )

p <- ggplot(plot_df, aes(geo_km, dxy)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE, colour = "#8f2be7", linewidth = 1.2) +
  xlab("Geographic distance (km)") +
  ylab(expression(D[XY]~"(divergence per site)")) +
  labs(title = sprintf("Mantel r = %.3f, p = %.4f",
                       mant$statistic, mant$signif)) +
  pub_theme
ggsave("ibd_ch.pdf", p, width = 9, height = 7.5)

write.csv(plot_df, "ibd_pairs.csv", row.names = FALSE)