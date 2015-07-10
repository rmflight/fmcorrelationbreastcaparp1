library(GenomicFiles)
library(BSgenome.Hsapiens.UCSC.hg19)
library(fmcorrelationbreastcaparp1)

genome_tiles <- tileGenome(seqinfo(Hsapiens), tilewidth = 500, cut.last.tile.in.chrom = TRUE)

genome_tiles <- genome_tiles[seqnames(genome_tiles) %in% seqnames(seqinfo(BigWigFile("../wgEncodeSydhNsomeGm12878Sig.bigWig")))]
seqlevels(genome_tiles) <- paste0("chr", c(seq(1, 22), "X", "M"))

genome_tiles <- genome_tiles[sample(length(genome_tiles), 10000)]
split_tiles <- split(genome_tiles, seqnames(genome_tiles))

MAP <- function(RANGE, FILE, ...) {
  if (length(RANGE) > 0){
    tmp <- import.bw(FILE, selection=RANGE, as="GRanges")
    tmp_cov <- coverage(tmp, weight = "score")
    tmp_range <- binned_function(RANGE, tmp_cov, "sum", "nuc")
  } else {
    tmp_range <- GRanges()
  }
  tmp_range
}

bw_counts <- reduceByRange(split_tiles, files = "../wgEncodeSydhNsomeGm12878Sig.bigWig", MAP)
bw_counts

nuc_counts <- GRanges()

for (icount in seq(1, length(bw_counts))){
  if (length(bw_counts[[icount]][[1]]) > 0){
    nuc_counts <- append(nuc_counts, bw_counts[[icount]][[1]])
  }
}

save(nuc_counts, file = "gm12878_nuc2.RData")

