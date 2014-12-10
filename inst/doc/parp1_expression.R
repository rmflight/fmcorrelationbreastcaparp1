## ----setupKnitr, echo=FALSE, results='hide'------------------------------
knitr::opts_chunk$set(dev='png')

## ----use_dir-------------------------------------------------------------
library(GenomicRanges)
library(ggplot2)
library(hgu133plus2.db)
library(fmanalysisbreastcaparp1)

## ----load_data-----------------------------------------------------------
data(tss_windows)
data(parp1_ln4_unique)
data(parp1_ln5_unique)
data(expr_data)

## ----ln4_ln5_counts------------------------------------------------------
mcf7_cov <- coverage(parp1_ln4_unique, weight = "n_count")
mdamb231_cov <- coverage(parp1_ln5_unique, weight = "n_count")

tss_windows <- binned_function(tss_windows, mcf7_cov, "sum", "mcf7_read")
tss_windows <- binned_function(tss_windows, mdamb231_cov, "sum", "mdamb231_read")

## ----translate_affy------------------------------------------------------
affy_2_ensembl <- select(hgu133plus2.db, keys = rownames(expr_data), columns = "ENSEMBLTRANS")
expr_data$ENSEMBL <- ""
expr_data <- expr_data[affy_2_ensembl[,1],]
expr_data$ENSEMBL <- affy_2_ensembl[,2]

## ----average_by_transcript-----------------------------------------------
mean_transcript <- tapply(expr_data$VALUE, expr_data$ENSEMBL, mean)

## ----add_to_tss----------------------------------------------------------
keep_transcript <- intersect(names(mean_transcript), names(tss_windows))
tss_windows <- tss_windows[keep_transcript]
mean_transcript <- mean_transcript[keep_transcript]
mcols(tss_windows)$mcf7_expr <- mean_transcript

## ----set_non_zero--------------------------------------------------------
non_zero <- "both"

## ----graphit-------------------------------------------------------------
expr_v_mcf7 <- subsample_nonzeros(mcols(tss_windows), c("mcf7_expr", "mcf7_read"), non_zero = non_zero, n_points = 10000)
ggplot(expr_v_mcf7, aes(x = mcf7_expr, y = mcf7_read)) + geom_point() + scale_y_log10() + scale_x_log10()
cor(log10(expr_v_mcf7[,1]+1), log10(expr_v_mcf7[,2]+1))

## ------------------------------------------------------------------------
out_cor <- correlate_non_zero(mcols(tss_windows), c("mcf7_expr", "mcf7_read"), log_transform = TRUE, non_zero = non_zero, test = TRUE)
out_cor

## ----save_tables---------------------------------------------------------
saveloc <- "../inst/correlation_tables"
cat(out_cor, file = file.path(saveloc, "expression_tss.txt"))

## ----session_info, echo=FALSE--------------------------------------------
Sys.time()
sessionInfo()

