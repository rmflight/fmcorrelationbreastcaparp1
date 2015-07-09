library(GenomicFiles)
library(BSgenome.Hsapiens.UCSC.hg19)

genome_tiles <- tileGenome(seqinfo(Hsapiens), tilewidth = 500, cut.last.tile.in.chrom = TRUE)

genome_tiles <- genome_tiles[seqnames(genome_tiles) %in% seqnames(seqinfo(BigWigFile("../wgEncodeSydhNsomeGm12878Sig.bigWig")))]
seqlevels(genome_tiles) <- paste0("chr", c(seq(1, 22), "X", "M"))

MAP <- function(RANGE, FILE, ...) {
  tmp <- import.bw(FILE, selection=RANGE, as="GRanges")
  unlist(sum(runValue(coverage(tmp, weight = "score")[RANGE])))
}

bw_counts <- reduceByRange(genome_tiles, files = "../wgEncodeSydhNsomeGm12878Sig.bigWig", MAP)
bw_counts2 <- unlist(bw_counts)
mcols(genome_tiles)[, "gm12878_nuc"] <- bw_counts2
save(genome_tile, file = "gm12878_nuc.RData")
