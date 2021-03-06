## ----setupKnitr, echo=FALSE, results='hide'------------------------------
knitr::opts_chunk$set(dev='png')

## ----load_packages-------------------------------------------------------
library(GenomicRanges)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)
library(fmcorrelationbreastcaparp1)
library(fmdatabreastcaparp1)

## ----load_tss_parp1------------------------------------------------------
data(tss_windows)
data(parp1_ln4_unique)
data(parp1_ln5_unique)

## ----parp1_counts--------------------------------------------------------
mcf7_cov <- coverage(parp1_ln4_unique, weight = "n_count")
mdamb231_cov <- coverage(parp1_ln5_unique, weight = "n_count")

tss_windows <- binned_function(tss_windows, mcf7_cov, "sum", "parp1_mcf7")
tss_windows <- binned_function(tss_windows, mdamb231_cov, "sum", "parp1_mdamb231")

## ----read_ctcf-----------------------------------------------------------
data(ctcf_rep1)
data(ctcf_rep2)
ctcf_r1_cov <- coverage(ctcf_rep1, weight = "mcols.signal")
ctcf_r2_cov <- coverage(ctcf_rep2, weight = "mcols.signal")
tss_windows <- binned_function(tss_windows, ctcf_r1_cov, "mean_nozero", "ctcf_r1")
tss_windows <- binned_function(tss_windows, ctcf_r2_cov, "mean_nozero", "ctcf_r2")

## ----set_non_zero--------------------------------------------------------
non_zero <- "both"

## ----graphit-------------------------------------------------------------
r1_v_mcf7 <- subsample_nonzeros(mcols(tss_windows), c("ctcf_r1", "parp1_mcf7"), non_zero = non_zero, n_points = 10000)
ggplot(r1_v_mcf7, aes(x = ctcf_r1, y = parp1_mcf7)) + geom_point() + scale_y_log10() + scale_x_log10()
cor(log10(r1_v_mcf7[,1]+1), log10(r1_v_mcf7[,2]+1))

r2_v_mcf7 <- subsample_nonzeros(mcols(tss_windows), c("ctcf_r2", "parp1_mcf7"), non_zero = non_zero, n_points = 10000)
ggplot(r2_v_mcf7, aes(x = ctcf_r2, y = parp1_mcf7)) + geom_point() + scale_y_log10() + scale_x_log10()
cor(log10(r2_v_mcf7[,1]+1), log10(r2_v_mcf7[,2]+1))

## ----grab_all------------------------------------------------------------
all_comb <- expand.grid(c("ctcf_r1", "ctcf_r2"), c("parp1_mcf7", "parp1_mdamb231"), stringsAsFactors = FALSE)
out_cor <- lapply(seq(1, nrow(all_comb)), function(i_row){
  #print(i_row)
  correlate_non_zero(mcols(tss_windows), as.character(all_comb[i_row,]), log_transform = TRUE, non_zero = non_zero, test = TRUE)
})
all_comb_names <- paste(all_comb[,1], all_comb[,2], sep = "_v_")
out_cor <- do.call(rbind, out_cor)
rownames(out_cor) <- all_comb_names

## ----tss_correlations, results='asis', echo=FALSE------------------------
knitr::kable(out_cor)

## ----graphs--------------------------------------------------------------
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
genome_tiles <- binned_function(genome_tiles, ctcf_r1_cov, "mean_nozero", "ctcf_r1")
genome_tiles <- binned_function(genome_tiles, ctcf_r2_cov, "mean_nozero", "ctcf_r2")

## ----genome_graph_test---------------------------------------------------
genome_r1_v_mcf7 <- subsample_nonzeros(mcols(genome_tiles), c("ctcf_r1", "parp1_mcf7"), non_zero = non_zero, n_points = 10000)
ggplot(genome_r1_v_mcf7, aes(x = ctcf_r1, y = parp1_mcf7)) + scale_x_log10() + scale_y_log10() + geom_point()
cor(log(genome_r1_v_mcf7[,1]+1), log(genome_r1_v_mcf7[,2]+1))

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
write.table(out_cor, file = file.path(saveloc, "ctcf_tss.txt"), sep = "\t")
write.table(genome_cor, file = file.path(saveloc, "ctcf_genome.txt"), sep = "\t")

## ----session_info, echo=FALSE--------------------------------------------
Sys.time()
sessionInfo()

