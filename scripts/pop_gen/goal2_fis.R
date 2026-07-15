#!/usr/bin/env Rscript
# Goal 2: per-population FIS (hierfstat) + pairwise FST (SNPRelate, W&C84).
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
col.pop    <- setNames(pal$col, pal$pop)
pop.levels <- pal$pop

pub_theme <- theme_test(base_size = 22) +
  theme(
    axis.title  = element_text(size = 26),
    axis.text   = element_text(size = 20, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 22),
    legend.text  = element_text(size = 18)
  )

snpgdsVCF2GDS(vcf.fn, "fis.gds", method = "biallelic.only",
              ignore.chr.prefix = chr.pre)
genofile <- snpgdsOpen("fis.gds")
samp.id  <- read.gdsn(index.gdsn(genofile, "sample.id"))
geno     <- snpgdsGetGeno(genofile, with.id = FALSE)   # individuals x SNP, 0/1/2

pop <- read.table(pop.fn, col.names = c("sample.id", "population"))
pop.vec <- pop$population[match(samp.id, pop$sample.id)]

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
res$population <- factor(res$population, levels = pop.levels)
res <- res[order(res$population), ]
write.csv(res, "fis_by_population.csv", row.names = FALSE)
print(res)

p_fis <- ggplot(res, aes(population, FIS, fill = population)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
  scale_fill_manual(values = col.pop, guide = "none") +
  xlab("Population") + ylab(expression(F[IS]~"(1 - "*H[o]*"/"*H[e]*")")) +
  pub_theme
ggsave("fis_by_population.pdf", p_fis, width = 9, height = 6.5)

pops    <- levels(factor(pop.vec))
n.pop   <- length(pops)
fst_mat <- matrix(NA, n.pop, n.pop, dimnames = list(pops, pops))

for (i in 1:(n.pop-1)) for (j in (i+1):n.pop) {
  keep <- samp.id[pop.vec %in% c(pops[i], pops[j])]
  grp  <- factor(pop.vec[pop.vec %in% c(pops[i], pops[j])])
  fst <- snpgdsFst(genofile,
                   sample.id = keep,
                   population = grp,
                   method = "W&C84",
                   autosome.only = FALSE,
                   verbose = FALSE)
  fst_mat[i, j] <- fst_mat[j, i] <- fst$Fst
}

closefn.gds(genofile)

write.csv(round(fst_mat, 4), "fst_pairwise.csv")
print(round(fst_mat, 4))

fst_long <- as.data.frame(as.table(fst_mat)) |>
  rename(pop1 = Var1, pop2 = Var2, fst = Freq) |>
  mutate(pop1 = factor(pop1, levels = pop.levels),
         pop2 = factor(pop2, levels = rev(pop.levels)))

p_fst <- ggplot(fst_long, aes(pop1, pop2, fill = fst)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(is.na(fst), "", sprintf("%.3f", fst))), size = 6) +
  scale_fill_gradient(low = "#f7f7f7", high = "#b2182b", na.value = "grey90",
                      name = expression(F[ST])) +
  xlab(NULL) + ylab(NULL) +
  coord_fixed() +
  pub_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("fst_pairwise.pdf", p_fst, width = 10, height = 8.5)

cat("Done. -> fis_by_population.csv/.pdf, fst_pairwise.csv/.pdf\n")