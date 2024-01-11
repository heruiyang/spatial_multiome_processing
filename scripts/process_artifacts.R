library(seqinr)
library(dplyr)
library(data.table)

# Parse arguments
args <- commandArgs(trailingOnly=TRUE)
genome_fasta <- args[1]
peaks_file <- args[2]
test_type <- args[3]
window <- args[4]
p_len <- args[5]
output_dir <- args[6]

genome <- read.fasta(genome_fasta)
all_peaks <- read.table(peaks_file)
peaks <- all_peaks[,c(1:3)]
colnames(peaks) <- c('chr','start','end')

# Given parameters, return logical vector of peaks to filter out
filter_peaks <- function(peaks, p_len, window, genome){
  pa_seq <- paste0(rep('a',p_len), collapse='')
  pt_seq <- paste0(rep('t',p_len), collapse='')


  artifact_peaks <- rep(FALSE, dim(peaks)[1])
  for (i in 1:dim(peaks)[1]){
    peak <- peaks[i,]
    upstream_region <- paste0(genome[[peak$chr]][(peak$start - window):(peak$start)], collapse='')
    downstream_region <- paste0(genome[[peak$chr]][(peak$end):(peak$end + window)], collapse='')


    if (grepl(pa_seq, upstream_region, fixed=T) || grepl(pt_seq, upstream_region, fixed=T) || grepl(pa_seq, downstream_region, fixed=T) || grepl(pa_seq, downstream_region, fixed=T)){
      artifact_peaks[i] <- TRUE
    }
  }

  return(artifact_peaks)
}


if (!all(unique(peaks$chr) %in% names(genome))) {
  stop('Chromosome names in genome fasta file do not match peak names. Make sure that the same genome sequence is used for both alignment and artifact filtering.')
}

if (test_type == 'debug'){
  p_len <- seq(5, 20, 1)
  window <- seq(20,170,50)

  mat <- matrix(0,nrow=length(p_len), ncol=length(window))
  colnames(mat) <- window
  rownames(mat) <- p_len

  for (i in 1:length(p_len)){
    for (j in 1:length(window)){
      mat[i, j] <- sum(filter_peaks(peaks=peaks, p_len=p_len[i], window=window[j], genome=genome))
    }
  }

  filename <- paste0(output_dir, '/artifact_debug_plot.pdf',collapse='')
  plot_df <- as.data.frame(mat)
  plot_df[['plen']] <- as.numeric(rownames(plot_df))
  plot_df <- as.data.frame(melt(as.data.table(plot_df), id.vars=c('plen'), value.name=c('n_filtered_out')))
  pdf(filename, width=7, height=5)
  ggplot(plot_df) +
    geom_line(aes(x=plen,y=n_filtered_out,group=variable,color=variable)) +
    labs(color = 'Window size') +
    xlab('Poly A/T length threshold') +
    ylab('No. of peaks filtered out') +
    theme_classic()
  dev.off()

} else {
  barcodes <- readLines(bc_file)
  counts <- Matrix::readMM(mtx_file)
  peaks_to_remove <- filter_peaks(peaks=peaks, p_len=p_len, window=window, genome=genome)

  print(paste0('Filtered out ',sum(peaks_to_remove),' peaks using parameters p_len=',p_len,' and window=',window,'.', collapse=''))

  # Filter matrix and save to output

  peaks_filtered <- all_peaks[!peaks_to_remove,]

  # peaks.bed
  write.table(peaks_filtered, paste0(output_dir, '/filtered_peaks.narrowPeak', collapse=''), quote=F, sep='\t', row.names=F, col.names=F)
}
