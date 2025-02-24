#' PCA Adjustment Function
#'
#' Performs PCA on the genotype matrix and returns adjusted covariates.
#'
#' @param X Genotype matrix.
#' @param C Covariate matrix (optional).
#' @param npcs Number of principal components to include.
#' @return List with adjusted covariates.
#' @export
pca_adjustment <- function(X, C = NULL, npcs = 5) {
  npcs <- min(npcs, ncol(X))
  pcs <- prcomp(X)$x[, 1:npcs, drop = FALSE]
  
  if (!is.null(C)) {
    combined <- cbind(C, pcs)
    rank_combined <- qr(combined)$rank
    rank_C <- qr(C)$rank
    
    if (rank_combined == rank_C) {
      warning("Some PCs are linearly dependent with covariates. Excluding dependent PCs.")
      pcs <- pcs[, (rank_C + 1):npcs, drop = FALSE]
    }
  }
  return(list(adjusted_covariates = if (!is.null(C)) cbind(C, pcs) else pcs, pcs = pcs))
}
