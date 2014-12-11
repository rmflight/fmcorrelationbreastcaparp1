## ----setupKnitr, echo=FALSE, results='hide'------------------------------
knitr::opts_chunk$set(dev='png')

## ----load_packages-------------------------------------------------------
library(GenomicRanges)
library(ggplot2)
library(fmcorrelationbreastcaparp1)
library(fmdatabreastcaparp1)
library(BSgenome.Hsapiens.UCSC.hg19)

## ----load_tss_parp1------------------------------------------------------
data(parp1_ln4_unique)
data(parp1_ln5_unique)
data(histone_marks)
data(tss_windows)

## ----mcf7_mdamb231_counts------------------------------------------------
mcf7_cov <- coverage(parp1_ln4_unique, weight = "n_count")
mdamb231_cov <- coverage(parp1_ln5_unique, weight = "n_count")

tss_windows <- binned_function(tss_windows, mcf7_cov, "sum", "parp1_mcf7")
tss_windows <- binned_function(tss_windows, mdamb231_cov, "sum", "parp1_mdamb231")

## ----average_histone_peaks-----------------------------------------------
for (i_name in names(histone_marks)){
  histone_cov <- coverage(histone_marks[[i_name]], weight = "mcols.signal")
  tss_windows <- binned_function(tss_windows, histone_cov, "mean_nozero", i_name)
}

## ----set_non_zero--------------------------------------------------------
non_zero <- "both"

## ----graphit-------------------------------------------------------------
h3k4me3k_v_mcf7 <- subsample_nonzeros(mcols(tss_windows), c("H3k4me3_r1", "parp1_mcf7"), non_zero = non_zero, n_points = 10000)
ggplot(h3k4me3k_v_mcf7, aes(x = H3k4me3_r1, y = parp1_mcf7)) + geom_point() + scale_y_log10() + scale_x_log10()
cor(log10(h3k4me3k_v_mcf7[,1]+1), log10(h3k4me3k_v_mcf7[,2]+1))

h3k27ac_v_mcf7 <- subsample_nonzeros(mcols(tss_windows), c("H3k27ac", "parp1_mcf7"), non_zero = non_zero, n_points = 10000)
ggplot(h3k27ac_v_mcf7, aes(x = H3k27ac, y = parp1_mcf7)) + geom_point() + scale_y_log10() + scale_x_log10()
cor(log10(h3k27ac_v_mcf7[,1]+1), log10(h3k27ac_v_mcf7[,2]+1))

## ----grab_all------------------------------------------------------------
all_comb <- expand.grid(names(histone_marks), c("parp1_mcf7", "parp1_mdamb231"), stringsAsFactors = FALSE)
out_cor <- lapply(seq(1, nrow(all_comb)), function(i_row){
  #print(i_row)
  correlate_non_zero(mcols(tss_windows), as.character(all_comb[i_row,]), log_transform = TRUE, non_zero = non_zero, test = TRUE)
})
all_comb_names <- paste(all_comb[,1], all_comb[,2], sep = "_v_")
out_cor <- do.call(rbind, out_cor)
rownames(out_cor) <- all_comb_names

## ----tss_correlations, results='asis', echo=FALSE------------------------
knitr::kable(out_cor)

## ----graphs, fig.keep='all'----------------------------------------------
out_graphs <- lapply(seq(1, nrow(all_comb)), function(i_row){
  use_vars <- as.character(all_comb[i_row,])
  subpoints <- subsample_nonzeros(mcols(tss_windows), use_vars, non_zero = non_zero, n_points = 10000)
  ggplot(subpoints, aes_string(x = use_vars[1], y = use_vars[2])) + geom_point() + scale_y_log10() + scale_x_log10()
})
out_graphs

## ----genome_tiles--------------------------------------------------------
genome_tiles <- tileGenome(seqinfo(Hsapiens), tilewidth = 2000, cut.last.tile.in.chrom = TRUE)

genome_tiles <- binned_function(genome_tiles, mcf7_cov, "sum", "parp1_mcf7")
genome_tiles <- binned_function(genome_tiles, mdamb231_cov, "sum", "parp1_mdamb231")

## ----genome_histone_marks------------------------------------------------
for (i_name in names(histone_marks)){
  histone_cov <- coverage(histone_marks[[i_name]], weight = "mcols.signal")
  genome_tiles <- binned_function(genome_tiles, histone_cov, "mean_nozero", i_name)
}

## ----genome_graph_test---------------------------------------------------
genome_h3k4me3k_v_mcf7 <- subsample_nonzeros(mcols(genome_tiles), c("H3k4me3_r1", "parp1_mcf7"), non_zero = non_zero, n_points = 10000)
ggplot(genome_h3k4me3k_v_mcf7, aes(x = H3k4me3_r1, y = parp1_mcf7)) + scale_x_log10() + scale_y_log10() + geom_point()
cor(log(genome_h3k4me3k_v_mcf7[,1]+1), log(genome_h3k4me3k_v_mcf7[,2]+1))

## ----genome_cor----------------------------------------------------------
genome_cor <- lapply(seq(1, nrow(all_comb)), function(i_row){
  #print(i_row)
  correlate_non_zero(mcols(genome_tiles), as.character(all_comb[i_row,]), log_transform = TRUE, non_zero = non_zero, test = TRUE)
})
all_comb_names <- paste(all_comb[,1], all_comb[,2], sep = "_v_")
genome_cor <- do.call(rbind, genome_cor)
rownames(genome_cor) <- all_comb_names

## ----genome_cor_display, results='asis', echo=FALSE----------------------
knitr::kable(genome_cor)

## ----save_tables---------------------------------------------------------
saveloc <- "../inst/correlation_tables"
write.table(out_cor, file = file.path(saveloc, "histone_marks_tss.txt"), sep = "\t")
write.table(genome_cor, file = file.path(saveloc, "histone_marks_genome.txt"), sep = "\t")

## ----session_info, echo=FALSE--------------------------------------------
Sys.time()
sessionInfo()

