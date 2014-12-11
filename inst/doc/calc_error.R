## ----setupKnitr, echo=FALSE, results='hide'------------------------------
knitr::opts_chunk$set(dev='png')

## ----setup---------------------------------------------------------------
graph_dir <- "/mlab/data/rmflight/Documents/projects/work/fondufe-mittendorf_lab/parp1/graphs"
library(GenomicRanges)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)
library(pracma)
library(fmcorrelationbreastcaparp1)
library(fmdatabreastcaparp1)

## ----plotX11, eval=FALSE-------------------------------------------------
#  X11()

## ----load_parp1----------------------------------------------------------
data(parp1_ln4_unique)
data(parp1_ln5_unique)
genome_tiles <- tileGenome(seqinfo(Hsapiens), tilewidth = 2000, cut.last.tile.in.chrom = TRUE)

## ----parp1_abundance-----------------------------------------------------
mcols(parp1_ln4_unique)$norm <- mcols(parp1_ln4_unique)$n_count / sum(mcols(parp1_ln4_unique)$n_count)
mcols(parp1_ln5_unique)$norm <- mcols(parp1_ln5_unique)$n_count / sum(mcols(parp1_ln5_unique)$n_count)
ln4_cov <- coverage(parp1_ln4_unique, weight = "n_count")
ln5_cov <- coverage(parp1_ln5_unique, weight = "n_count")
genome_tiles <- binned_function(genome_tiles, ln4_cov, "sum", "parp1_mcf7")
genome_tiles <- binned_function(genome_tiles, ln5_cov, "sum", "parp1_mdamb231")

## ----sample--------------------------------------------------------------
non_zero <- "both"
r1_v_r2 <- subsample_nonzeros(mcols(genome_tiles), c("parp1_mcf7", "parp1_mdamb231"), non_zero = non_zero, n_points = 10000)

ggplot(r1_v_r2, aes(x = parp1_mcf7, y = parp1_mdamb231)) + geom_point() + scale_y_log10() + scale_x_log10()

## ----odr-----------------------------------------------------------------
both_non <- find_non_zeros(mcols(genome_tiles), c("parp1_mcf7", "parp1_mdamb231"), log_transform = TRUE, non_zero = non_zero)
parp1_data <- mcols(genome_tiles[both_non])
parp1_data[,1] <- log10(parp1_data[,1] + 1)
parp1_data[,2] <- log10(parp1_data[,2] + 1)
parp1_odr <- odregress(parp1_data[,1], parp1_data[,2])

## ----plot_everything-----------------------------------------------------
plot(parp1_data[,1], parp1_data[,2], asp = 1)
lines(parp1_data[,1], parp1_odr$fitted, col = "red")

## ----save_parp1_everything, eval=FALSE-----------------------------------
#  png(file = file.path(graph_dir, "parp1_reps.png")); plot(parp1_data[,1], parp1_data[,2], asp = 1); lines(parp1_data[,1], parp1_odr$fitted, col = "red"); dev.off();

## ----parp1_centile_correlation-------------------------------------------
parp1_cor_quantiles <- run_cum_quantiles(parp1_data)

odregress_fun <- function(x, y){
  tmp <- pracma::odregress(x, y)
  tmp$coef[1]
}

parp1_odr_quantiles <- run_cum_quantiles(parp1_data, similarity = odregress_fun)

## ----plotthem, fig.keep='all'--------------------------------------------
p <- ggplot(parp1_cor_quantiles, aes(x = cut, y = cor, color = type)) + geom_point()
print(p)
p + facet_grid(type ~ ., scales = "free")

## ----load_methyl---------------------------------------------------------
data(methyl_rep1)
data(methyl_rep2)

mcols(methyl_rep1)$methyl_read <- mcols(methyl_rep1)$mcols.readCount * mcols(methyl_rep1)$mcols.percentMeth / 100
mcols(methyl_rep2)$methyl_read <- mcols(methyl_rep2)$mcols.readCount * mcols(methyl_rep2)$mcols.percentMeth / 100

methyl_r1_cov <- coverage(methyl_rep1, weight = "methyl_read")
methyl_r2_cov <- coverage(methyl_rep2, weight = "methyl_read")
genome_tiles <- binned_function(genome_tiles, methyl_r1_cov, "sum", "methyl_r1")
genome_tiles <- binned_function(genome_tiles, methyl_r2_cov, "sum", "methyl_r2")

## ----methyl_transform----------------------------------------------------
methyl_nonzero <- find_non_zeros(mcols(genome_tiles), c("methyl_r1", "methyl_r2"), log_transform = TRUE, non_zero = non_zero)

methyl_data <- mcols(genome_tiles[methyl_nonzero])[,c("methyl_r1", "methyl_r2")]
methyl_data[,1] <- log10(methyl_data[,1] + 1)
methyl_data[,2] <- log10(methyl_data[,2] + 1)

## ----methyl_plot---------------------------------------------------------
plot(methyl_data[,1], methyl_data[,2], asp = 1)

## ----methyl_ord----------------------------------------------------------
methyl_ord <- odregress(methyl_data[,1], methyl_data[,2])
plot(methyl_data[,1], methyl_data[,2], asp = 1)
lines(methyl_data[,1], methyl_ord$fitted, col = "red")

## ----save_methyl_plot, eval = FALSE--------------------------------------
#  png(file = file.path(graph_dir, "methyl_reps.png")); plot(methyl_data[,1], methyl_data[,2], asp = 1); lines(methyl_data[,1], methyl_ord$fitted, col = "red"); dev.off();

## ----methyl_centiles-----------------------------------------------------
methyl_cor_quantiles <- run_cum_quantiles(methyl_data)
methyl_odr_quantiles <- run_cum_quantiles(methyl_data, similarity = odregress_fun)

## ----methyl_plot2--------------------------------------------------------
m <- ggplot(methyl_cor_quantiles, aes(x = cut, y = cor, color = type)) + geom_point()
print(m)
m + facet_grid(type ~ ., scales = "free")

## ----load_ctcf_data------------------------------------------------------
data(ctcf_rep1)
data(ctcf_rep2)

ctcf_r1_cov <- coverage(ctcf_rep1, weight = "mcols.signal")
ctcf_r2_cov <- coverage(ctcf_rep2, weight = "mcols.signal")

genome_tiles <- binned_function(genome_tiles, ctcf_r1_cov, "mean_nozero", "ctcf_r1")
genome_tiles <- binned_function(genome_tiles, ctcf_r2_cov, "mean_nozero", "ctcf_r2")

## ----ctcf_nozero---------------------------------------------------------
ctcf_nonzero <- find_non_zeros(mcols(genome_tiles), c("ctcf_r1", "ctcf_r2"), log_transform = TRUE, non_zero = non_zero)

ctcf_data <- mcols(genome_tiles[ctcf_nonzero])[, c("ctcf_r1", "ctcf_r2")]
ctcf_data[,1] <- log10(ctcf_data[,1] + 1)
ctcf_data[,2] <- log10(ctcf_data[,2] + 1)

## ----ctcf_plot-----------------------------------------------------------
plot(ctcf_data[,1], ctcf_data[,2], asp = 1)

## ----ctcf_cum------------------------------------------------------------
ctcf_cor_quantiles <- run_cum_quantiles(ctcf_data)
ctcf_odr_quantiles <- run_cum_quantiles(ctcf_data, similarity = odregress_fun)

## ----ctcf_plots----------------------------------------------------------
cplot <- ggplot(ctcf_cor_quantiles, aes(x = cut, y = cor, color = type)) + geom_point()
print(cplot)
cplot + facet_grid(type ~ ., scales = "free")

## ----session_info, echo=FALSE--------------------------------------------
Sys.time()
sessionInfo()

