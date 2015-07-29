## ----load_data, message=FALSE--------------------------------------------
library(ggplot2)
library(GenomicFiles)
library(BSgenome.Hsapiens.UCSC.hg19)
library(fmcorrelationbreastcaparp1)
library(fmdatabreastcaparp1)
data("histone_marks") # mcf-7 histone modifications
data("parp1_ln4_unique") # mcf-7 parp1 reads

## ----tile----------------------------------------------------------------
genome_tiles <- tileGenome(seqinfo(Hsapiens), tilewidth = 500, cut.last.tile.in.chrom = TRUE)

## ----parp1_counts--------------------------------------------------------
min_sum <- 25
parp1_cov <- coverage(parp1_ln4_unique, weight = "n_count")
genome_tiles <- binned_function(genome_tiles, parp1_cov, "sum", "parp1")
mcols(genome_tiles)$parp1[mcols(genome_tiles)$parp1 < min_sum] <- 0
mcols(genome_tiles)$parp1[mcols(genome_tiles)$parp1 >= min_sum] <- 1

## ----nucleosome_counts---------------------------------------------------
mcols(genome_tiles)$histmod <- 0
histone_tiles <- genome_tiles
for (i_name in names(histone_marks)){
  histone_cov <- coverage(histone_marks[[i_name]], weight = "mcols.signal")
  histone_tiles <- binned_function(histone_tiles, histone_cov, "sum", i_name)
  mcols(genome_tiles)$histmod[mcols(histone_tiles)[, i_name] > 0] <- 1 
}

## ----sample--------------------------------------------------------------
n_sample <- 1000
out_counts <- lapply(seq(1, 1000), function(x){
  use_loc <- sample(length(genome_tiles), n_sample)
  tmp <- colSums(as.matrix(mcols(genome_tiles[use_loc])))
  data.frame(counts = tmp, which = names(tmp))
})
out_counts <- do.call(rbind, out_counts)

## ----plotit--------------------------------------------------------------
ggplot(out_counts, aes(x = counts, fill = which)) + geom_histogram(alpha = 0.8, binwidth = 1, position = "identity")

## ----count_specificity---------------------------------------------------
presence <- mcols(genome_tiles)
parp_only <- sum(as.logical(presence$parp1) & !as.logical(presence$histmod))
histmod_only <- sum(!as.logical(presence$parp1) & as.logical(presence$histmod))
both <- sum(as.logical(presence$parp1) & as.logical(presence$histmod))

## ----load_nuc_signal-----------------------------------------------------
data(nuc_signal)
genome_tiles <- subsetByOverlaps(genome_tiles, nuc_signal, type = "equal")

mcols(genome_tiles)[, "gm12878"] <- 0
mcols(genome_tiles)[mcols(nuc_signal)$nuc > 100, "gm12878"] <- 1

## ----overlaps------------------------------------------------------------
library(limma)
presence2 <- as.matrix(mcols(genome_tiles))
nuc_overlap_counts <- vennCounts(presence2)
attr(nuc_overlap_counts, "class") <- NULL

## ------------------------------------------------------------------------
knitr::kable(nuc_overlap_counts)

