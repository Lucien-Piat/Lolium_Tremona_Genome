#!/usr/bin/env Rscript
# Usage: 01_pca_structure.R <vcf> <pop> <outdir> <chr_prefix> <ploidy> <pal_pop> <pal_anc> <K> <recompute>

suppressMessages({
  library(SNPRelate)
  library(LEA)
  library(tidyverse)
})

args      <- commandArgs(trailingOnly = TRUE)
vcf.fn    <- args[1]
pop.fn    <- args[2]
outdir    <- args[3]
chr.pre   <- if (length(args) >= 4) args[4] else ""
ploidy    <- if (length(args) >= 5) as.integer(args[5]) else 2L
pal.pop.fn<- args[6]
pal.anc.fn<- args[7]
K.best    <- if (length(args) >= 8) as.integer(args[8]) else 6L
recompute <- if (length(args) >= 9) as.logical(args[9]) else TRUE

setwd(outdir)

k.range <- 3:12
n.rep   <- 10
n.cpu   <- 4

#  colonne 1 = couleur, colonne 2 = cle 
pal.pop <- read.table(pal.pop.fn, comment.char = "", col.names = c("col", "pop"))
pal.anc <- read.table(pal.anc.fn, comment.char = "", col.names = c("col", "k"))
col.pop <- setNames(pal.pop$col, pal.pop$pop)
col.anc <- pal.anc$col

# VCF vers GDS 
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

pub_theme <- theme_test(base_size = 22) +
  theme(
    axis.title   = element_text(size = 26),
    axis.text    = element_text(size = 20, colour = "black"),
    legend.title = element_text(size = 22),
    legend.text  = element_text(size = 20)
  )

# Scree plot
scree_df <- data.frame(PC = factor(seq_len(10)), var = pc.pct[seq_len(10)])
ggsave("scree_plot.pdf",
       ggplot(scree_df, aes(PC, var)) +
         geom_col(fill = "steelblue") +
         geom_text(aes(label = paste0(var, "%")), vjust = -0.4, size = 6) +
         xlab("Principal component") + ylab("Variance explained (%)") +
         pub_theme,
       width = 9, height = 6)

# PCA plots
plot_pca <- function(df, pcx, pcy, ix, iy, file) {
  p <- ggplot(df, aes(.data[[pcx]], .data[[pcy]], colour = population)) +
    geom_point(size = 3.5, alpha = 0.85) +
    scale_colour_manual(values = col.pop) +
    xlab(paste0(pcx, " (", pc.pct[ix], "%)")) +
    ylab(paste0(pcy, " (", pc.pct[iy], "%)")) +
    labs(colour = "Population") +
    pub_theme
  ggsave(file, p, width = 9, height = 7)
}
plot_pca(pca_df, "PC1", "PC2", 1, 2, "PC1-PC2_population.pdf")
plot_pca(pca_df, "PC2", "PC3", 2, 3, "PC2-PC3_population.pdf")
plot_pca(pca_df, "PC3", "PC4", 3, 4, "PC3-PC4_population.pdf")
plot_pca(pca_df, "PC4", "PC5", 4, 5, "PC4-PC5_population.pdf")
plot_pca(pca_df, "PC5", "PC6", 5, 6, "PC5-PC6_population.pdf")

#  Structure
ped2geno("lmf.pruned.ped")
if (recompute || !file.exists("lmf.pruned.snmfProject")) {
  project <- snmf("lmf.pruned.geno", K = k.range, entropy = TRUE,
                  repetitions = n.rep, project = "new",
                  CPU = n.cpu, ploidy = ploidy)
} else {
  cat("Rechargement du projet sNMF existant (pas de recalcul)\n")
  project <- load.snmfProject("lmf.pruned.snmfProject")
}

pdf("crossentropy.pdf"); plot(project, col = "steelblue", pch = 19, cex = 1.2); dev.off()

best  <- which.min(cross.entropy(project, K = K.best))
Q.mat <- as.matrix(Q(project, K = K.best, run = best))
rownames(Q.mat) <- samp.id
write.table(Q.mat, paste0("Q_K", K.best, ".tbl"), quote = FALSE)

pop.order  <- pop$population[match(samp.id, pop$sample.id)]
pop.levels <- names(col.pop)
ord <- order(factor(pop.order, levels = pop.levels))
Q.ord <- Q.mat[ord, ]; lab.ord <- samp.id[ord]; pop.ord <- pop.order[ord]

pdf(paste0("structure_K", K.best, ".pdf"), width = 18, height = 7)
par(mar = c(9, 6, 3, 2), cex.lab = 2.2, cex.axis = 1.6)

barplot(t(Q.ord), col = col.anc[1:K.best], border = NA, space = 0,
        names.arg = lab.ord, las = 2, cex.names = 0.9,
        xlab = "", ylab = "Ancestry proportion")
sep <- cumsum(table(factor(pop.ord, levels = pop.levels)))
abline(v = sep[-length(sep)], col = "white", lwd = 2)
mid <- (c(0, sep[-length(sep)]) + sep) / 2
axis(3, at = mid, labels = pop.levels, tick = FALSE, line = -0.5, cex.axis = 1.8)
mtext("Individuals", side = 1, line = 7, cex = 2.2)
dev.off()
