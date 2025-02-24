#' Compute Minor Allele Frequency (MAF)
#'
#' This function calculates the Minor Allele Frequency (MAF) for each SNP in the genotype matrix.
#' It ensures that all columns are numeric before computing allele frequencies.
#'
#' @param X A genotype matrix (taxa in rows, SNPs in columns)
#' @return A numeric vector of MAF values for each SNP
#' @export

compute_maf <- function(X) {
  # Convert data to numeric
  X_maf <- X  # Create a copy
  for (col in colnames(X_maf)) {
    X_maf[[col]] <- suppressWarnings(as.numeric(as.character(X_maf[[col]])))
  }
  
  # Convert to matrix to avoid errors
  X_maf <- as.matrix(X_maf)
  # Compute allele frequency
  allele_freq <- colMeans(X_maf, na.rm = TRUE) / 2
  # Compute Minor Allele Frequency (MAF)
  maf <- pmin(allele_freq, 1 - allele_freq)
  
  return(maf)
}
