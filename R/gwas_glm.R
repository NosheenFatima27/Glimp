#' Perform GWAS using a Generalized Linear Model (GLM)
#'
#' This function performs GWAS using a GLM, with optional PCA adjustment for population structure.
#'
#' @param y A numeric vector of phenotypes.
#' @param X A numeric matrix of genotype data (SNPs in columns).
#' @param C A numeric matrix of covariates (optional).
#' @param include_pcs Logical. If TRUE, includes principal components as covariates.
#' @param npcs Number of principal components to include (default: 3).
#' @return A numeric vector of p-values for each SNP.
#' @examples
#' # Simulated example
#' set.seed(123)
#' y <- rnorm(100)  # Simulated phenotype data
#' X <- matrix(rbinom(1000, 2, 0.5), nrow = 100, ncol = 10)  # Simulated genotype data
#' C <- matrix(rnorm(100 * 3), nrow = 100, ncol = 3)  # Simulated covariates
#' pvals <- gwas_glm(y, X, C, include_pcs = TRUE, npcs = 3)
#' head(pvals)  # Show first few p-values
#'
#' @export
gwas_glm <- function(y, X, C, include_pcs = TRUE, npcs = 3) {
  if (length(y) != nrow(X)) stop("Length of y must match the number of rows in X.")
  if (!is.null(C) && nrow(C) != length(y)) stop("Number of rows in C must match length of y.")
  
  if (include_pcs) {
    C <- pca_adjustment(X, C, npcs)$adjusted_covariates
  }
  
  m <- ncol(X)
  pvals <- numeric(m)
  
  for (i in 1:m) {
    model_data <- data.frame(y = y, marker = X[, i])
    if (!is.null(C)) model_data <- cbind(model_data, C)
    
    fit <- tryCatch(
      glm(y ~ ., data = model_data, family = gaussian()),
      error = function(e) {
        warning(paste("GLM failed for marker", i, "- setting p-value to 1"))
        return(NULL)
      }
    )
    
    if (!is.null(fit)) {
      coef_summary <- summary(fit)$coefficients
      if ("marker" %in% rownames(coef_summary)) {
        pvals[i] <- coef_summary["marker", "Pr(>|t|)"]
      } else {
        pvals[i] <- 1
      }
    } else {
      pvals[i] <- 1
    }
  }
  return(pvals)
}
