library(jsonlite)
library(Matrix)
library(dplyr)

# Parse arguments
args <- commandArgs(trailingOnly=TRUE)
fiducial_json <- args[1]
bc_file <- args[2]
peaks_file <- args[3]
mtx_file <- args[4]
spot_coords <- args[5]
output_dir <- args[6]

json_data <- fromJSON(txt=fiducial_json)
barcodes <- readLines(bc_file)
peaks <- read.table(peaks_file)
counts <- Matrix::readMM(mtx_file)
visium_coords <- read.table(spot_coords, sep='\t', row.names=1)

# Select spots in tissue
colnames(visium_coords) <- c('col','row')
visium_coords[['barcode']] <- rownames(visium_coords)

merged_res <- dplyr::left_join(visium_coords, json_data$oligo, by=c('col','row'))
barcodes_keep <- merged_res[!is.na(merged_res$tissue),'barcode']

print(paste0('Kept ',length(barcodes_keep),' barcodes in tissue.', collapse=''))

# Create filtered matrix
colnames(counts) <- barcodes
  
counts_out <- counts[,barcodes_keep]
bc_out <- barcodes_keep
peaks_out <- peaks[rowSums(as.matrix(counts_out)) > 0,]

# Write output
# barcodes.tsv
write.table(bc_out, paste(output_dir, '/barcodes.tsv', sep=''), quote=F, sep='\t', row.names=F, col.names=F)
# peaks.bed
write.table(peaks_out, paste(output_dir, '/peaks.bed', sep=''), quote=F, sep='\t', row.names=F, col.names=F)
# matrix.mtx
writeMM(counts_out, paste(output_dir, '/matrix.mtx', sep=''))
