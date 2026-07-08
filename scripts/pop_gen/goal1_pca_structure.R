#!/usr/bin/env Rscript

suppressMessages({
  library(SNPRelate)
  library(LEA)
  library(tidyverse)
})

args    <- commandArgs(trailingOnly = TRUE)
vcf.fn  <- args[1]
pop.fn  <- args[2]
outdir  <- args[3]
chr.pre <- if (length(args) >= 4) args[4] else ""
ploidy  <- if (length(args) >= 5) as.integer(args[5]) else 2L

setwd(outdir)

## Parametres structure
k.range <- 1:12
n.rep   <- 10
n.cpu   <- 4

# Conversion VCF vers GDS
snpgdsVCF2GDS(vcf.fn, "lmf.gds", method = "biallelic.only",
              ignore.chr.prefix = chr.pre)

genofile <- snpgdsOpen("lmf.gds")
samp.id  <- read.gdsn(index.gdsn(genofile, "sample.id"))

# PCA 
my_pca <- snpgdsPCA(genofile, autosome.only = FALSE, num.thread = n.cpu)
pc.pct <- round(my_pca$varprop * 100, 2)

snpgdsGDS2PED(genofile, "lmf.pruned")
closefn.gds(genofile)

pop <- read.table(pop.fn, col.names = c("sample.id", "population"))

pca_df <- data.frame(sample.id = my_pca$sample.id,
                     PC1 = my_pca$eigenvect[, 1],
                     PC2 = my_pca$eigenvect[, 2],
                     PC3 = my_pca$eigenvect[, 3],
                     PC4 = my_pca$eigenvect[, 4],
                     PC5 = my_pca$eigenvect[, 5],
                     PC6 = my_pca$eigenvect[, 6]) |>
  left_join(pop, by = "sample.id")

write.csv(pca_df, "pca_table.csv", row.names = FALSE)

scree_df <- data.frame(PC = factor(seq_len(10)),
                       var = pc.pct[seq_len(10)])
p_scree <- ggplot(scree_df, aes(PC, var)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = paste0(var, "%")), vjust = -0.4, size = 3) +
  xlab("Composante principale") +
  ylab("Variance expliquee (%)") +
  theme_test()
ggsave("scree_plot.pdf", p_scree, width = 8, height = 5)

pops    <- sort(unique(pca_df$population))
n.pop   <- length(pops)
my.cols <- c("#96127d", "#5DC863", "#FDE725", "#3B528B", "#21908c",
             "#e07b39", "#540101", "#ff0000")[seq_len(n.pop)]
names(my.cols) <- pops

plot_pca <- function(df, pcx, pcy, ix, iy, file) {
  p <- ggplot(df, aes(.data[[pcx]], .data[[pcy]], colour = population)) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_colour_manual(values = my.cols) +
    xlab(paste0(pcx, " : ", pc.pct[ix], "%")) +
    ylab(paste0(pcy, " : ", pc.pct[iy], "%")) +
    theme_test()
  ggsave(file, p, width = 8, height = 6)
}

plot_pca(pca_df, "PC1", "PC2", 1, 2, "PC1-PC2_population.pdf")
plot_pca(pca_df, "PC2", "PC3", 2, 3, "PC2-PC3_population.pdf")
plot_pca(pca_df, "PC3", "PC4", 3, 4, "PC3-PC4_population.pdf")
plot_pca(pca_df, "PC4", "PC5", 4, 5, "PC4-PC5_population.pdf")
plot_pca(pca_df, "PC5", "PC6", 5, 6, "PC5-PC6_population.pdf")

#  Structure
ped2geno("lmf.pruned.ped")

project <- snmf("lmf.pruned.geno",
                K = k.range, entropy = TRUE, repetitions = n.rep,
                project = "new", CPU = n.cpu, ploidy = ploidy)

# Selection de K 
pdf("crossentropy.pdf")
plot(project, col = "steelblue", pch = 19, cex = 1.2)
dev.off()

K.best <- 6
best   <- which.min(cross.entropy(project, K = K.best))
Q.mat  <- as.matrix(Q(project, K = K.best, run = best))
rownames(Q.mat) <- samp.id
write.table(Q.mat, paste0("Q_K", K.best, ".tbl"), quote = FALSE)

pop.order  <- pop$population[match(samp.id, pop$sample.id)]
pop.levels <- c("GULF",  "L46", "L60", "PR", "SLB", "L31", "CH", "TREM")
ord     <- order(factor(pop.order, levels = pop.levels))

Q.ord   <- Q.mat[ord, ]
lab.ord <- samp.id[ord]
pop.ord <- pop.order[ord]

my.colors <- c("#96127d", "#5DC863", "#FDE725", "#3B528B", "#21908c",
               "#e07b39", "#440154", "#31688e", "#c94c4c", "#8fbc8f",
               "#4682b4", "#d1a3ff")

pdf(paste0("structure_K", K.best, ".pdf"), width = 16, height = 5)
bp <- barplot(t(Q.ord), col = my.colors[1:K.best],
              border = NA, space = 0,
              xlab = "Individus", ylab = "Proportions d'ascendance",
              names.arg = lab.ord, las = 2, cex.names = 0.4)

sep <- cumsum(table(factor(pop.ord, levels = pop.levels)))
abline(v = sep[-length(sep)], col = "white", lwd = 1.5)

mid <- (c(0, sep[-length(sep)]) + sep) / 2
axis(3, at = mid, labels = pop.levels, tick = FALSE, line = -0.5, cex.axis = 0.8)
dev.off()

cat("Termine. Sorties dans", outdir, "\n")