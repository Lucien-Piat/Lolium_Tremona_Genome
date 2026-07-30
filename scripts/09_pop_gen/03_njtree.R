#!/usr/bin/env Rscript
# Usage: 03_njtree.R <vcf> <pop> <outdir> <chr_prefix> <pal_pop>

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

# Population and color mapping
pop <- read.table(pop.fn, col.names = c("sample.id", "population"))
tip.pop <- pop$population[match(tree$tip.label, pop$sample.id)]
tip.col <- col.pop[tip.pop]

edge.col <- rep("grey30", nrow(tree$edge))
tip.idx  <- match(seq_along(tree$tip.label), tree$edge[, 2])
edge.col[tip.idx] <- tip.col

# Phylogram 
pdf("njtree_phylo.pdf", width = 10, height = 18)
plot(tree, type = "phylogram",
     tip.color = tip.col, edge.color = edge.col,
     align.tip.label = TRUE, 
     cex = 0.45, no.margin = FALSE, label.offset = 0.002)
tiplabels(pch = 15, col = tip.col, cex = 0.9, offset = 0.001)
dev.off()

# Fan 
pdf("njtree_fan.pdf", width = 16, height = 16)
plot(tree, type = "fan", 
     tip.color = tip.col, 
     edge.color = edge.col,
     edge.width = 1.5, 
     cex = 2.0, 
     no.margin = TRUE, 
     label.offset = 0.01)
tiplabels(pch = 15, col = tip.col, cex = 1.5, offset = 0.004)
legend("topright", legend = names(col.pop), col = col.pop,
       pch = 15, pt.cex = 1.5, cex = 1.2, bty = "n", inset = 0.02)

dev.off()