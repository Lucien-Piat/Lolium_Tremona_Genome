#!/usr/bin/env Rscript
# Goal 2 : FIS par population (hierfstat). Palette externe.
# Usage: goal2_fis.R <vcf> <pop> <outdir> <chr_prefix> <pal_pop>

suppressMessages({
  library(SNPRelate)
  library(hierfstat)
  library(tidyverse)
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

snpgdsVCF2GDS(vcf.fn, "fis.gds", method = "biallelic.only",
              ignore.chr.prefix = chr.pre)
genofile <- snpgdsOpen("fis.gds")
samp.id  <- read.gdsn(index.gdsn(genofile, "sample.id"))
geno <- snpgdsGetGeno(genofile, with.id = FALSE)   # individus x SNP, 0/1/2
closefn.gds(genofile)

pop <- read.table(pop.fn, col.names = c("sample.id", "population"))
pop.vec <- pop$population[match(samp.id, pop$sample.id)]

# Recodage format hierfstat (11/12/22)
recode <- geno
recode[geno == 0] <- 11
recode[geno == 1] <- 12
recode[geno == 2] <- 22
recode[is.na(geno)] <- NA
dat <- data.frame(pop = as.integer(factor(pop.vec)), recode)

bs <- basic.stats(dat)
lvl <- levels(factor(pop.vec))
res <- data.frame(population = lvl,
                  Ho  = round(colMeans(bs$Ho,  na.rm = TRUE), 4),
                  He  = round(colMeans(bs$Hs,  na.rm = TRUE), 4),
                  FIS = round(colMeans(bs$Fis, na.rm = TRUE), 4))
res <- res[order(res$population), ]
write.csv(res, "fis_by_population.csv", row.names = FALSE)
print(res)

p <- ggplot(res, aes(reorder(population, FIS), FIS, fill = population)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = col.pop, guide = "none") +
  xlab("Population") + ylab("FIS (1 - Ho/He)") +
  theme_test() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("fis_by_population.pdf", p, width = 8, height = 5)

cat("Termine. FIS -> fis_by_population.csv / .pdf\n")