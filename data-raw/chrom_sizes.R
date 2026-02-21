# Script to download and process GRCh37 and GRCh38 chromosome lengths

# UCSC chromosome sizes
# hg38 (GRCh38)
# hg19 (GRCh37)

fetch_ucsc_sizes <- function(url) {
  df <- read.table(url, header = FALSE, col.names = c("chrom", "size"), stringsAsFactors = FALSE)
  # Keep only main chromosomes
  main_chroms <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
  df <- df[df$chrom %in% main_chroms, ]
  
  # Strip 'chr' for the named vector to match user style
  names_clean <- sub("^chr", "", df$chrom)
  names_clean[names_clean == "M"] <- "MT"
  
  sizes <- df$size
  names(sizes) <- names_clean
  
  # Standard order
  order_main <- c(as.character(1:22), "X", "Y", "MT")
  sizes <- sizes[match(order_main, names(sizes))]
  sizes <- sizes[!is.na(sizes)]
  
  sizes
}

message("Fetching GRCh38 (hg38) sizes...")
grch38 <- fetch_ucsc_sizes("http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes")

message("Fetching GRCh37 (hg19) sizes...")
grch37 <- fetch_ucsc_sizes("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes")

# Save to data/
usethis::use_data(grch38, grch37, overwrite = TRUE)
