#' Generate Manhattan Plot
#'
#' This function creates a Manhattan plot for visualizing GWAS results.
#'
#' @param pvals Vector of p-values from GWAS.
#' @param myGM Data frame with SNP names, chromosome, and position.
#' @export
manhattan_plot <- function(pvals, myGM) {
  plot_data <- data.frame(
    SNP = myGM$SNP,
    Chromosome = factor(myGM$Chromosome),
    Position = myGM$Position,
    P = pvals
  )
  
  plot_data$logP <- -log10(plot_data$P)
  threshold <- -log10(0.05 / nrow(myGM))
  
  ggplot(plot_data, aes(x = Chromosome, y = logP, color = Chromosome)) +
    geom_jitter(width = 0.4, alpha = 0.7, size = 1.5) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
    labs(title = "Manhattan Plot", x = "Chromosome", y = "-log10(p-value)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
}

#' Generate QQ Plot
#'
#' This function creates a QQ plot to assess the distribution of p-values.
#'
#' @param pvals Vector of p-values from GWAS.
#' @export
qq_plot <- function(pvals) {
  obs <- -log10(sort(pvals))
  exp <- -log10(ppoints(length(pvals)))
  
  ggplot(data.frame(Expected = exp, Observed = obs), aes(x = Expected, y = Observed)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(title = "QQ Plot", x = "Expected -log10(p-value)", y = "Observed -log10(p-value)") +
    theme_minimal()
}
