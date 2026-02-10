#' Compute response enrichment score for each analysis unit
#'
#' This function calculates a normalized enrichment score
#' (AUC-like score) for each analysis unit based on ranked gene
#' expression and pre-trained gene differential weights
#' for a given cancer type.
#'
#' @param mat A numeric expression matrix,
#'   with genes as rows and analysis units as columns.
#' @param topN Integer. Number of top-ranked genes used to compute
#'   the cumulative enrichment score. Default is 3000.
#' @param nGenes Integer. Total number of detected genes used for
#'   position weighting. Default is 10000.
#' @param cancer Character string specifying cancer type.
#'   Must be one of the supported cancer codes
#'   ("NSCLC","MEL","HCC","BLCA","RCC","BC","GC","ESCA","HNSC","GBM","CRC").
#'
#' @returns A numeric vector of normalized enrichment scores,
#'   named by analysis unit IDs.
#'
#' @export
Response_Score <- function(mat, topN = 1500, nGenes = 10000, cancer) {

  d <- .load_cancer_data_internal(cancer)
  gene_Diff <- d$gene_Diff

  genes <- intersect(rownames(mat), names(gene_Diff))
  mat <- mat[genes, , drop = FALSE]
  gene_Diff <- gene_Diff[genes]

  gene_list <- rowMeans(as.matrix(mat))
  gene_list <- sort(gene_list, decreasing = TRUE)
  gene_list <- names(gene_list)

  perm_estimates <- calculate_deviation(
    gene_Diff,
    nperm = 3000,
    gene_list = gene_list,
    topN = topN
  )

  sample_rankings <- build_rankings_simple(mat)

  AUCs <- numeric(ncol(mat))
  names(AUCs) <- colnames(mat)

  for (i in seq_len(ncol(mat))) {

    order_idx <- sample_rankings[, i, drop = TRUE]
    order_idx <- as.numeric(order_idx)
    names(order_idx) <- rownames(sample_rankings)

    order_idx <- order_idx[order_idx <= topN]
    topN_genes <- names(order_idx)[order(order_idx)]

    cumulative_scores <- cumsum(gene_Diff[topN_genes])

    result <- calculate_enrichment_score(
      cumulative_scores,
      random_mean = perm_estimates$random_mean,
      random_sd   = perm_estimates$random_sd,
      nGenes      = nGenes
    )

    AUCs[i] <- result$normalized_score
  }

  return(AUCs)
}

