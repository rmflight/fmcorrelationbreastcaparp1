#' calculate Rle mean without zeros
#' 
#' Given an \code{Rle} with possible zeros, first remove the zeros, and then calculate the weighted mean
#' 
#' @param in_rle the \code{Rle} to calculate the non-zero mean for
#' @return a numeric value
#' @export
rle_mean_nozero <- function(in_rle){
  n_values <- length(in_rle@values)
  zero_vals <- in_rle@values == 0
  
  if ((sum(zero_vals) == n_values) || (in_rle@lengths == 0)){
    mean_val <- 0
  } else {
    mean_val <- weighted.mean(in_rle@values[!zero_vals], in_rle@lengths[!zero_vals])
  }
  
  return(mean_val)
}

#' calculate binned values
#' 
#' given a set of genomic bins and an RLE-list object, return value
#' 
#' Adapted from documentation of \code{tileGenome}
#' 
#' @param bins a \code{GRanges} object representing genomic bins
#' @param numvar a named \code{RleList} object representing numerical variable along the genome
#' @param binfun character string of the function to apply
#' @param mcolname the name of the metadata column containing the binned value of the object
#' @export
#' @return bins with the \code{mcolname} containing the binned value
binned_function <- function(bins, numvar, binfun = "mean", mcolname = "avg"){
  stopifnot(is(bins, "GRanges"))
  stopifnot(is(numvar, "RleList"))

  bin_names <- levels(factor(as.character(seqnames(bins))))
  numvar_names <- names(numvar)

  calc_names <- intersect(bin_names, numvar_names)

  stopifnot(length(calc_names) > 0)

  zero_names <- setdiff(bin_names, numvar_names)
  has_bin_names <- FALSE
  if (length(names(bins)) != 0){
	  has_bin_names <- TRUE
  }

  
  bins_per_seqname <- split(ranges(bins), as.character(seqnames(bins)))
  
  calc_list <- mclapply(calc_names,
                         function(seqname){
                           views <- Views(numvar[[seqname]],
                                         bins_per_seqname[[seqname]])
                           switch(binfun,
                                  mean_nozero = viewApply(views, rle_mean_nozero),
                                  mean = viewMeans(views),
                                  sum = viewSums(views),
                                  min = viewMins(views),
                                  max = viewMaxs(views))
                         })

  zero_list <- lapply(zero_names,
		      function(seqname){
			      zero_bin <- bins_per_seqname[[seqname]]
			      n_zero_bin <- length(zero_bin)
			      out_bin <- vector("numeric", n_zero_bin)
			      if (has_bin_names){
				      names(out_bin) <- names(zero_bin)
			      }
		      	      out_bin
		      })
  names(calc_list) <- calc_names
  names(zero_list) <- zero_names
  bin_list <- c(calc_list, zero_list)
  bin_list <- bin_list[bin_names]
  new_mcol <- unsplit(bin_list, as.factor(as.character(seqnames(bins))))
  mcols(bins)[[mcolname]] <- new_mcol
  return(bins)
}

#' convert to counted GRanges
#' 
#' given a set of reads in a delimited file, find the unique locations, count the number of reads at each location,
#' apply a hard threshold if required, and save the data to an RData file.
#' 
#' @param reads_file the file with the reads
#' @param delim the delimiter to use
#' @param read_start the column containing the read starts
#' @param width how wide are the reads
#' @param scramble_strand should the strand be scrambled (default is \code{TRUE})
#' @param max_count what is the maximum value that the counts should be
#' @param prepend text to prepend the out_file name with
#' @import GenomicRanges
#' @import magrittr
#' @export
#' @return the RData filename
convert_to_counted <- function(reads_file, delim = ",", read_start = "startx", 
                               width = 1, strand = "strand", scramble_strand = TRUE, max_count = 6, prepend = ""){
  
  tmp_reads <- read.table(reads_file, sep = delim, header = TRUE, stringsAsFactors = FALSE)
  
  use_strand <- tmp_reads[, strand]
  
  if (scramble_strand){
    use_strand <- "*"
  }
  
  use_chr <- get_chr(reads_file, "_")
  
  tmp_locs <- GRanges(seqnames = use_chr,
                      ranges = IRanges(start = tmp_reads[, read_start], width = width),
                      strand = use_strand)
  
  unique_locs <- unique(tmp_locs)
  unique_overlap <- countOverlaps(unique_locs, tmp_locs)
  unique_overlap[(unique_overlap > max_count)] <- max_count
  mcols(unique_locs)[, "count"] <- unique_overlap
  
  data_path <- dirname(reads_file)
  out_file <- paste(prepend, use_chr, sep = "_") %>% paste(., ".RData", sep = "")
  save(unique_locs, file = file.path(data_path, out_file))
  return(out_file)
}

#' subsample non-zeros
#' 
#' Takes a \code{DataFrame} instance and for the two columns indicated, returns a subset of points that are non-zero in one or both of 
#' the variables. 
#' 
#' @param data a \code{DataFrame}
#' @param data_columns which columns to use
#' @param log_transform whether or not to log-transform the data
#' @param non_zero "either" one of the columns, or "both"
#' @param n_points how many points to sample
#' @export
#' @return data.frame
subsample_nonzeros <- function(data, data_columns, log_transform = TRUE, non_zero = "either", n_points = 1000){
  nz_index <- find_non_zeros(data, data_columns, log_transform, non_zero)
  
  sample_index <- sample(nz_index, size = n_points, replace = FALSE)
  return(as.data.frame(data[sample_index, data_columns]))
}

#' weighted average by width
#' 
#' given a \code{GRanges} object and an \code{mcols} column, a weighted average of the \code{mcols} column will
#' be calculated where the weight will be the fractional number of bases indicated by the width of each entry in 
#' the \code{GRanges} object.
#' 
#' @param granges the \code{GRanges} object
#' @param use_column which column from \code{mcols} has the signal to average
#' @return numeric
#' @export
average_value <- function(granges, use_column = "mcols.signal"){
  if (length(granges) == 0){
    return(0)
  }
  range_width <- width(granges)
  total_width <- sum(range_width)
  range_weights <- range_width / total_width
  
  range_value <- mcols(granges)[,use_column]
  return(weighted.mean(range_value, range_weights))
}

#' get overlap counts
#' 
#' given a set of files (assumed to be chromosomes), and a set of \code{GRanges}, count the number of overlaps with
#' the \code{GRanges}.
#' 
#' @param file_list the list of files to read in corresponding to ranges
#' @param granges the genomic ranges to count the overlaps with
#' @param offset a value by which to \code{shift} the reads by, positively for "+" and negatively for "-"
#' 
#' @return vector of counts
#' @export
#' @import GenomicRanges
get_overlap_counts <- function(file_list, granges, offset = 0){
  out_overlap <- mclapply(file_list, function(x){
    load(x)
    if (!(exists("unique_locs"))){
      stop("unique_locs does not exist!", call. = FALSE)
    }
    
    ul_granges_overlap <- findOverlaps(unique_locs, granges)
    ul_granges_list <- split(queryHits(ul_granges_overlap), subjectHits(ul_granges_overlap))
    
    unlist(lapply(ul_granges_list, function(in_overlap){
      sum(mcols(unique_locs[in_overlap])[, "count"])
    }))
  })
  
  
  out_tss <- unlist(out_overlap)
  names(out_tss) <- names(granges)[as.numeric(names(out_tss))]
  return(out_tss)
}

#' get chromosome
#' 
#' from a filename, get the chromosome part of the filename
#' 
#' @param filename the filename to parse
#' @param split_chr the character separating the chr* from the rest
#' @return string of the chromosome
#' @export
get_chr <- function(filename, split_chr){
  file_split <- strsplit(filename, split_chr)[[1]]
  n_split <- length(file_split)
  chr_part <- file_split[n_split]
  substr(chr_part, 1, nchar(chr_part)-4)
}

#' correlate non-zeros
#' 
#' from a \code{DataFrame} object, generate correlation of the non-zero entries
#' 
#' @param data the data we are working with
#' @param data_columns which columns to use
#' @param log_transform do a log transformation on the data before calculating the correlation
#' @param non_zero "either" or "both" values need to be non-zero
#' @param test also return a p-value of the correlation?
#' @return numeric vector
#' @export
correlate_non_zero <- function(data, data_columns, log_transform = TRUE, non_zero = "either", test = TRUE){
  
  nz_index <- find_non_zeros(data, data_columns, log_transform, non_zero)
    
  x <- data[nz_index, data_columns[1]]
  y <- data[nz_index, data_columns[2]]
  
  if (log_transform){
    x <- log(x + 1)
    y <- log(y + 1)
  }
  
  c_value <- cor(x, y)
  
  
  if (!test){
    return(corr_value = c_value)
  } else {
    t_value <- cor.test(x, y)$p.value
    return(c(corr_value = c_value, p_value = t_value))
  }
}

#' find non-zeros
#' 
#' given a \code{data.frame}, find the non-zero entries in each case
#' 
#' @param data the data we are working with
#' @param data_columns which columns to use
#' @param log_transform do a log transformation on the data before identifying zeros
#' @param non_zero "either" or "both", how much data *must* be non-zero
#' 
#' @export
#' @return indices into the original data that are not zero
find_non_zeros <- function(data, data_columns, log_transform = TRUE, non_zero = "either"){
  
  # this list holds the non-zero, and other stuff
  non_zero_entries <- lapply(data_columns, function(x){
    tmp_data <- data[, x]
    if (log_transform){
      tmp_data <- log(tmp_data + 1)
    }
    (tmp_data != 0) & (!(is.infinite(tmp_data))) & (!(is.na(tmp_data))) & (!(is.nan(tmp_data))) 
  })
  
  # name them so we can access them
  names(non_zero_entries) <- data_columns
  keep_data <- non_zero_entries[data_columns]
  
  non_zero_logical <- switch(non_zero,
                             either = do.call("|", non_zero_entries),
                             both = do.call("&", non_zero_entries))
  
  nz_index <- which(non_zero_logical)
  names(nz_index) <- NULL
  
  return(nz_index)
}

#' run bootstrapped correlation
#' 
#' given a set of points for correlations, generate bootstrapped samples and calculate a 95%CI or SD
#' 
#' @param data the data we are working with
#' @param data_columns which columns to use
#' @param log_transform do a log transformation on the data before calculating the correlation
#' @param n_boot how many bootstrap samples to generate
#' @return correlation value and standard deviation
#' @export
bootstrap_correlation <- function(data, data_columns, log_transform = TRUE, non_zero = "either", n_boot = 1000){
  nz_index <- find_non_zeros(data, data_columns, log_transform, non_zero)
  
  x <- data[nz_index, data_columns[1]]
  y <- data[nz_index, data_columns[2]]
  
  if (log_transform){
    x <- log(x + 1)
    y <- log(y + 1)
  }
  
  c_value <- cor(x, y)
  
  n_sample <- length(x)
  
  bootstrap_cor <- lapply(seq(1, n_boot), function(in_rep){
    use_boot <- sample(n_sample, n_sample, replace = TRUE)
    cor(x[use_boot], y[use_boot])
  })
  bootstrap_cor <- unlist(bootstrap_cor)
  sd_cor <- sd(bootstrap_cor)
  ci_val <- sd_cor * 1.96
  ci_range <- c(c_value - ci_val, c_value + ci_val)
  return(c(cor = c_value, ci = ci_range))
}

#' boostrap orthogonal regression
#' 
#' for a set of data, do a bootstrapped orthogonal regression
#'
#' @param x_data the X-data
#' @param y_data the Y-data
#' @param n_boot how many bootstrap samples to do
#' @export 
#' @return vector of the coefficients
#' @importFrom pracma odregress
bootstrap_odregress <- function(x_data, y_data, n_boot = 500){
  n_point <- length(x_data)
  
  boot_data <- lapply(seq(1, n_boot), function(x){
		      use_sample <- sample(n_point, n_point, replace = TRUE)
		      odres <- odregress(x_data[use_sample], y_data[use_sample])
		      odres$coeff
  })
  boot_data <- do.call(rbind, boot_data)
  return(boot_data)
}

#' quantile index
#' 
#' Given a set of values, generate the indices of the data for supplied quantiles
#' 
#' @param x the data to generate quantiles for
#' @param n_quantile how many quantiles to split into
#' @export
#' @importFrom stats quantile
#' @return list with indices for each quantile
generate_quantile_indices <- function(x, n_quantile = 101){
  gen_quantile <- unique(quantile(x, seq(0, 1, length.out = n_quantile)))
  x_split <- cut(x, gen_quantile, right = FALSE)
  x_index <- split(seq(1, length(x)), x_split)
  return(x_index)
}

#' cumulative indices
#' 
#' Given a list of indices assumed to in ranked order, generate a set of cumulative
#' indices wherein the indices are essentially collapsed together.
#' 
#' @param indices_list the list of indices we want to work with
#' @param direction which way to put the indices together (low_to_high or high_to_low)
#' @examples
#' indices_list <- list("1" = c(1,2,3), "2" = c(4,5,6), "3" = c(7,8,9))
#' cum_indices(indices_list, "low_to_high")
#' cum_indices(indices_list, "high_to_low")
#' @export
#' @return new indices_list with cumulative indices
cum_indices <- function(indices_list, direction = "low_to_high"){
  n_indices <- length(indices_list)
  
  seq_list <- switch(direction,
                     low_to_high = seq(1, n_indices, 1),
                     high_to_low = seq(n_indices, 1, -1))
  
  new_indices <- lapply(seq(1, n_indices), function(index){
    grab_index <- seq_list[1:index]
    unlist(indices_list[grab_index], use.names = FALSE)
  })
  names(new_indices) <- names(indices_list)[seq_list]
  return(new_indices)
}

#' list correlation
#' 
#' Given a data.frame of data, and a list of indices, calculate the correlation
#' between the columns of the data.frame using the indices, iterating over the list
#' of indices.
#' 
#' @param x the data.frame of data (should be two columns only)
#' @param indices_list the list of indices to be iterated over
#' @param similarity the function used to calculate the similarity
#' @export
#' @return numerical vector of correlations
#' @importFrom stats cor
list_correlation <- function(x, indices_list, similarity = stats::cor){
  out_cor <- lapply(indices_list, function(in_ind){
    similarity(x[in_ind, 1], x[in_ind, 2])
  })
  out_cor <- do.call(c, out_cor)
  return(out_cor)
}

#' cut names 2 value
#' 
#' Given a set of \code{cuts} as character names of a vector, return either the \code{left},
#' \code{right} or \code{mid} as a numerical value that can be used for plotting.
#' 
#' @param cuts the set of character \code{cuts} to process
#' @param type which type of value to return
#' 
#' @examples
#' 
#' # set up some fake data
#' cuts <- c("[0.301,0.4771)", "[0.4771,0.6276)")
#' 
#' cut_2_value(cuts, "left")
#' cut_2_value(cuts, "right")
#' cut_2_value(cuts, "mid")
#' 
#' @return numerical vector
#' @export
cut_2_value <- function(cuts, type = "left"){
  return_function <- switch(type,
                            left = function(x){as.numeric(x[1])},
                            right = function(x){as.numeric(x[2])},
                            mid = function(x){mean(as.numeric(x))})
  
  cut_sub <- gsub("\\[", "", cuts)
  cut_sub <- gsub("\\)", "", cut_sub)
  
  cut_split <- strsplit(cut_sub, ",")
  
  out_val <- sapply(cut_split, return_function)
  return(out_val)
}

#' run adding quantiles
#' 
#' Given a data set, generate the quantiles, the \code{cum_indices}, run the correlation,
#' and return a data.frame for plotting and diagnostics.
#' 
#' @param x a two column data.frame
#' @param n_quantile how many quantiles
#' @param similarity_function what similarity function to use
#' @param cut_loc which way to generate cut values on the data ("left", "right", "mid")
#' 
#' @export
#' @return data.frame with the values, the corresponding cut value, and which cumulative indices generated it
#' 
#' @importFrom stats cor
run_cum_quantiles <- function(x, n_quantile = 101, similarity = stats::cor, cut_loc = "left"){
  x_mean <- rowMeans(as.matrix(x))
  
  x_q <- generate_quantile_indices(x_mean, n_quantile)
  x_low2hi <- cum_indices(x_q, "low_to_high")
  x_hi2low <- cum_indices(x_q, "high_to_low")
  
  cor_low2hi <- list_correlation(x, x_low2hi, similarity)
  cor_hi2low <- list_correlation(x, x_hi2low, similarity)
  
  cut_low2hi <- cut_2_value(names(cor_low2hi), type = cut_loc)
  cut_hi2low <- cut_2_value(names(cor_hi2low), type = cut_loc)
  
  n_low2hi <- length(cor_low2hi)
  n_hi2low <- length(cor_hi2low)
  
  data.frame(cut = c(cut_low2hi, cut_hi2low),
             cor = c(cor_low2hi, cor_hi2low),
             type = c(rep("low2hi", n_low2hi), rep("hi2low", n_hi2low)))
}