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
}

#' subsample non-zeros
#' 
#' Takes a \code{DataFrame} instance and for the two columns indicated, returns a subset of points that are non-zero in one or both of 
#' the variables. 
#' 
#' @param data a \code{DataFrame}
#' @param data_columns which columns to use
#' @param non_zero either one of the columns, or "both"
#' @param n_points how many points to sample
#' @export
#' @return data.frame
subsample_nonzeros <- function(data, data_columns, non_zero = "both", n_points = 1000){
  non_zero_entries <- lapply(data_columns, function(x){
    data[, x] != 0
  })
  names(non_zero_entries) <- data_columns
  
  if (non_zero == "both"){
    use_data <- data_columns
  } else {
    use_data <- non_zero
  }
  
  keep_data <- non_zero_entries[use_data]
  
  non_zero_entries <- do.call("&", non_zero_entries)
  nz_index <- which(non_zero_entries)
  
  if (length(nz_index) < n_points){
    n_points <- length(nz_index)
  }
  
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
    tmp <- read.table(x, header = TRUE, sep = ",")
    tmpRange <- GRanges(seqnames = get_chr(x),
                        ranges = IRanges(start = tmp$startx, width = 1),
                        strand = tmp$strand)
    
    if (offset != 0){
      pos_loc <- strand(tmpRange) == "+"
      tmpRange[pos_loc] <- shift(tmpRange[pos_loc], offset)
      tmpRange[!pos_loc] <- shift(tmpRange[!pos_loc], -1 * offset)
    }
    countOverlaps(tss_windows, tmpRange, ignore.strand = TRUE)
  })
  
  out_tss <- lapply(out_overlap, function(in_count){
    in_count[in_count != 0]
  })
  out_tss <- unlist(out_tss)
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
#' @param test also return a p-value of the correlation?
#' @return numeric vector
#' @export
correlate_non_zero <- function(data, data_columns, log_transform = TRUE, test = TRUE){
  non_zero_entries <- lapply(data_columns, function(x){
    data[, x] != 0
  })
  
  names(non_zero_entries) <- data_columns
  keep_data <- non_zero_entries[data_columns]
  
  non_zero_entries <- do.call("&", non_zero_entries)
  nz_index <- which(non_zero_entries)
  
  x <- data[nz_index, data_columns[1]]
  y <- data[nz_index, data_columns[2]]
  
  if (log_transform){
    x <- log(x)
    y <- log(y)
  }
  
  c_value <- cor(x, y)
  
  
  if (!test){
    return(corr_value = c_value)
  } else {
    t_value <- cor.test(x, y)$p.value
    return(c(corr_value = c_value, p_value = t_value))
  }
}

