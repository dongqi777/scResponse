#' Plot enrichment score curve for a single analysis unit
#'
#' This function visualizes the cumulative enrichment score curve
#' for a specified analysis unit, highlighting the dominant deviation,
#' enrichment direction, and contributing gene ranks.
#'
#' An analysis unit can represent a single cell in single-cell data,
#' a spatial spot in spatial transcriptomics, or any other expression unit
#' represented by a column in the input matrix.
#'
#' @param mat A numeric expression matrix,
#'   with genes as rows and analysis units as columns.
#' @param cell_id Character string specifying the target analysis unit ID
#'   (must match a column name in \code{mat}).
#' @param cancer Character string specifying cancer type used to
#'   load pre-trained gene differential weights.
#'   Cancer type must be one of the supported cancer codes
#'   ("NSCLC","MEL","HCC","BLCA","RCC","BC","GC","ESCA","HNSC","GBM","CRC").
#' @param topN Integer. Number of top-ranked genes used to build
#'   the enrichment curve. Default is 3000.
#' @param nGenes Integer. Total number of detected genes used for
#'   position weighting. Default is 10000.
#'
#' @returns Invisibly returns \code{NULL}. The function is called
#'   for its side effect of generating a plot.
#'
#' @export
Plot_EnrichmentScore <- function(mat, cell_id, cancer,
                                 topN = 1500, nGenes = 10000) {

  if (!cell_id %in% colnames(mat)) {
    stop("cell_id not found in matrix colnames.", call. = FALSE)
  }

  d <- .load_cancer_data_internal(cancer)
  gene_Diff <- d$gene_Diff

  genes <- intersect(rownames(mat), names(gene_Diff))
  if (length(genes) < 10) {
    stop("Too few overlapping genes between expression matrix and gene_Diff.", call. = FALSE)
  }

  mat <- mat[genes, , drop = FALSE]
  gene_Diff <- gene_Diff[genes]

  expr_vec <- mat[, cell_id, drop = TRUE]
  expr_vec <- as.numeric(expr_vec)
  names(expr_vec) <- genes

  gene_list <- names(sort(expr_vec, decreasing = TRUE))

  ## ====== deviation ======
  perm_estimates <- calculate_deviation(
    gene_Diff,
    nperm = 3000,
    gene_list = gene_list,
    topN = topN
  )

  order_idx <- order(expr_vec, decreasing = TRUE)
  topN <- min(topN, length(order_idx))
  topN_genes <- genes[order_idx[seq_len(topN)]]

  cumulative_scores <- cumsum(gene_Diff[topN_genes])

  result <- calculate_enrichment_score(
    cumulative_scores,
    random_mean = perm_estimates$random_mean,
    random_sd   = perm_estimates$random_sd,
    nGenes      = nGenes
  )

  n <- length(result$cumulative_scores)
  x <- seq_len(n)
  y <- result$cumulative_scores

  main_color <- if (result$direction == "positive") "#456882" else "#C75D2C"

  par(mar = c(5, 4, 4, 4) + 0.1)
  plot(x, y, type = "l", col = main_color, lwd = 3,
       xlab = "Gene Rank", ylab = "Cumulative Score",
       main = sprintf(
         "Enrichment Score = %.3f\nDirection: %s, Max %s Dev = %.1f @ Rank %d",
         result$normalized_score,
         result$direction,
         result$direction,
         result$max_deviation,
         result$max_position
       ))

  if (result$positive_area > 0) {
    polygon(c(x, rev(x)), c(pmax(y, 0), rep(0, n)),
            col = "#C4E1E6", border = NA)
  }
  if (result$negative_area > 0) {
    polygon(c(x, rev(x)), c(pmin(y, 0), rep(0, n)),
            col = "#DCA06D", border = NA)
  }

  abline(h = 0, col = "gray", lty = 2, lwd = 1.5)
  abline(v = result$max_position, col = "#8A784E", lty = 2, lwd = 1.5)
  points(result$max_position, result$max_deviation,
         pch = 19, col = "#8A0000", cex = 1.5)

  legend(
    if (result$direction == "positive") "topleft" else "bottomleft",
    legend = c("Cumulative Score", "Max Deviation", "Positive Area", "Negative Area"),
    col = c(main_color, "#8A0000", "#C4E1E6", "#DCA06D"),
    lty = c(1, NA, NA, NA),
    lwd = c(2, NA, NA, NA),
    pch = c(NA, 19, 15, 15),
    pt.cex = c(NA, 1.5, 2, 2),
    bg = adjustcolor("white", alpha.f = 0.5),
    border = NA
  )

  invisible(NULL)
}




