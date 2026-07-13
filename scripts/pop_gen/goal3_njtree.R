#!/usr/bin/env Rscript
# Usage: goal3_njtree.R <vcf> <pop> <outdir> <chr_prefix> <pal_pop>

suppressMessages({
  library(SNPRelate)
  library(ape)
})

args    <- commandArgs(trailingOnly = TRUE)
vcf.fn  <- args[1]
pop.fn  <- args[2]
outdir  <- args[3]
chr.pre <- if (length(args) >= 4) args[4] else ""
pal.fn  <- args[5]

setwd(outdir)

pal <- read.table(pal.fn, comment.char = "", col.names = c("col", "pop"))
col.pop <- setNames(pal$col, pal$pop)

snpgdsVCF2GDS(vcf.fn, "nj.gds", method = "biallelic.only",
              ignore.chr.prefix = chr.pre)
genofile <- snpgdsOpen("nj.gds")
samp.id  <- read.gdsn(index.gdsn(genofile, "sample.id"))
ibs <- snpgdsIBS(genofile, autosome.only = FALSE, num.thread = 4)
closefn.gds(genofile)

dist.mat <- as.dist(1 - ibs$ibs)
attr(dist.mat, "Labels") <- samp.id
tree <- nj(dist.mat)
tree <- ladderize(tree)
write.tree(tree, "njtree.nwk")

# Population et couleur de chaque feuille
pop <- read.table(pop.fn, col.names = c("sample.id", "population"))
tip.pop <- pop$population[match(tree$tip.label, pop$sample.id)]
tip.col <- col.pop[tip.pop]

edge.col <- rep("grey30", nrow(tree$edge))
tip.idx  <- match(seq_along(tree$tip.label), tree$edge[, 2])
edge.col[tip.idx] <- tip.col

pdf("njtree_phylo.pdf", width = 10, height = 18)
plot(tree, type = "phylogram",
     tip.color = tip.col, edge.color = edge.col,
     align.tip.label = TRUE, 
     cex = 0.45, no.margin = FALSE, label.offset = 0.002)

tiplabels(pch = 15, col = tip.col, cex = 0.9, offset = 0.001)

dev.off()

pdf("njtree_fan.pdf", width = 11, height = 11)
plot(tree, type = "fan", tip.color = tip.col, edge.color = edge.col,
     cex = 0.5, no.margin = TRUE, label.offset = 0.004)
tiplabels(pch = 15, col = tip.col, cex = 1.1, offset = 0.002)
legend("topright", legend = names(col.pop), col = col.pop,
       pch = 15, cex = 0.9, bty = "n")
dev.off()
