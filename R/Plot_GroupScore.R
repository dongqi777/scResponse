#' Plot group-wise distribution of enrichment scores
#'
#' Visualize the distribution of enrichment (response) scores across
#' predefined cell groups using violin and box plots.
#'
#' This function operates directly on a metadata \code{data.frame},
#' making it independent of any specific single-cell framework
#' (e.g. Seurat). Each row in \code{meta} is assumed to represent
#' one cell or sample.
#'
#' @param meta A \code{data.frame} containing per-cell metadata.
#'   Rows correspond to cells or samples.
#' @param group Character scalar. Column name in \code{meta} defining
#'   the grouping variable on the x-axis.
#'   Default is \code{"CellType"}.
#' @param value Character scalar. Column name in \code{meta} containing
#'   enrichment scores to be visualized.
#'   Default is \code{"RScore"}.
#' @param fill Character scalar. Column name in \code{meta} used to
#'   determine fill colors of the violins and boxplots.
#'   Default is \code{"Response"}.
#' @param ViolinWidth Numeric. Width of the violin plots.
#'   Default is \code{0.9}.
#' @param ViolinAlpha Numeric. Transparency level of the violin plots.
#'   Default is \code{0.9}.
#' @param BoxplotWidth Numeric. Width of the boxplots.
#'   Default is \code{0.2}.
#' @param BoxAlpha Numeric. Transparency level of the boxplots.
#'   Default is \code{1}.
#' @param BoxColor Character. Color of the boxplot outlines.
#'   Default is \code{"white"}.
#' @param bgcolor Character. Background color of the plot panel.
#'   Default is \code{"#F1F0E4"}.
#'
#' @details
#' Groups on the x-axis are ordered by the median value of the
#' enrichment score within each group, in decreasing order.
#'
#' @return A \code{ggplot2} object.
#'
#' @export
Plot_GroupScore <- function(meta, group= "CellType",  value= "RScore" , fill="Response",
                            ViolinWidth=0.9,  ViolinAlpha=0.9,
                            BoxplotWidth=0.2, BoxAlpha=1 , BoxColor="white" ,bgcolor="#F1F0E4"   ) {
  # -------- 安全检查 --------
  if (!is.data.frame(meta)) {
    stop("`meta` must be a data.frame.", call. = FALSE)
  }

  for (col in c(group, value, fill)) {
    if (!col %in% colnames(meta)) {
      stop(sprintf("Column '%s' not found in meta.", col), call. = FALSE)
    }
  }

  colors=c("#DE8F5F",    "#8FA31E",    "#c1d7ae" ,    "#C5B0CD" ,   "#9DD2CBFF" ,  "#789461" ,   "#e1ca96",    "#91adc2",    "#C1856D",    "#F5C9B0" ,
           "#39737C",    "#928167FF",  "#896C6C",     "#538C48FF",  "#BFCD70FF",   "#580c1f",    "#82a3a1",    "#BCA88D" ,   "#617764FF" , "#D6A99D" ,
           "#A6B28B" ,   "#7D8D86",    "#7E748DFF" ,  "#9CAFAA",    "#aa998f",     "#78B9B5",    "#B67352"  ,  "#A6B28B",    "#DEE791",    "#d1d2f9" ,
           "#a9ddd6" ,   "#9ba0bc" ,   "#ffb977",     "#acd98d",    "#c1b8c8",     "#32a251",    "#A9D1E1"  ,  "#82853b"    ,"#ccc94d",    "#3cb7cc" ,
           "#b85a0d",    "#FFE797" ,   "#98d9e4" ,    "#D1A980",    "#5E936C"  ,   "#FFDE63"   , "#93DA97"  ,  "#86bbd8"    ,"#F5CBCB"    ,"#EAC8A6" )

  sortgroup <- tapply(meta[[value]], meta[[group]], stats::median,na.rm = TRUE)
  meta[[group]] <- factor(
    meta[[group]],
    levels = names(sort(sortgroup, decreasing = TRUE))
  )

  n_fill <- length(unique(meta[[fill]]))
  if (n_fill > length(colors)) {
    stop(
      sprintf("Number of fill groups (%d) exceeds available colors (%d).",
              n_fill, length(colors)),
      call. = FALSE
    )
  }

  ggplot(
    meta,
    aes(
      x = .data[[group]],
      y = .data[[value]],
      fill = .data[[fill]]
    )
  ) +
    geom_violin(
      position = position_dodge(width = ViolinWidth),
      scale = "width",
      alpha = ViolinAlpha,
      drop = FALSE
    ) +
    geom_boxplot(
      alpha = BoxAlpha,
      outlier.size = 0,
      size = 0.4,
      position = position_dodge(0.9),
      width = BoxplotWidth,
      color = BoxColor
    )+
    scale_fill_manual(
      values = colors[1:length(unique(meta[, fill]))],
      name = fill
    )+
    theme_light()+
    theme(#axis.line = element_line(colour = "black", linewidth  = 1),
      axis.text.x = element_text(angle = 45, hjust = 1 , color="black"),   #axis.line.x.bottom = element_line(color="black"),
      axis.text.y = element_text(color="black"),    #  axis.line.y.left = element_line(color="black"),
      panel.background = element_rect(fill = bgcolor )
    )  + labs(x = group, y = value )
}
