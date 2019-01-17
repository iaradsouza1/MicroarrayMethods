#' PCA plot
#'
#' @param eset ExpressionSet object with normalized expression values. Commonly the result of RMA normalization
#' @param components A vector of length 2 containing which principal components to plot
#' @param groups Column name corresponding to the groups of samples annotation on pheno data provided to the \code{import_celfiles()} function
#' @param type Column name corresponding to the types of samples (control and treated samples) on pheno data provided to the \code{import_celfiles()} function (optional)
#' @param batch Column name corresponding to each replicates on pheno data provided to the \code{import_celfiles()} function (optional)
#' @param filename File name to save plot (pdf format)
#' @param ... Other arguments to be passed to \code{pdf()} function
#'
#' @return Plot of chosen PCs
#' @export
#'
#' @examples
plot_pca <- function(eset, components, groups, type=NULL, batch=NULL, filename=NULL, ...) {

  data_exprs <- as.data.frame(na.omit(exprs(eset)))
  tdata_exprs <- t(data_exprs)
  pca <- prcomp(tdata_exprs)
  importance <- summary(pca)$importance
  importance <- round(importance[2, c(components[1], components[2])], digits = 3) * 100

  pca_df <- as.data.frame(pca$x[, c(components[1], components[2])])
  pca_df <- cbind(pca_df, pData(eset))

  if (!is.null(filename)) {
    if(is.null(type) && is.null(batch)) {
      pdf(filename, ...)
      print(ggplot2::ggplot(pca_df, ggplot2::aes_string(x = pca_df[,1], y = pca_df[,2], col = groups)) +
        ggplot2::labs(title = "PCA after normalization", x = paste0("PC ", components[1], " (", importance[1], "%)"),
             y = paste0("PC ", components[2], " (", importance[2], "%)"), col = "Group") +
        ggplot2::geom_point(size = 2) +
        ggplot2::scale_color_brewer(palette = "Set1", type = "qual") +
        ggplot2::theme_bw())
      dev.off()
    } else if (!is.null(type) && is.null(batch)) {
      pdf(filename, ...)
      print(ggplot2::ggplot(pca_df, ggplot2::aes_string(x = pca_df[,1], y = pca_df[,2], col = groups, shape = type)) +
        ggplot2::labs(title = "PCA after normalization", x = paste0("PC ", components[1]," (", importance[1], "%)"),
             y = paste0("PC ", components[2], " (", importance[2], "%)"), col = "Group", shape = "Type") +
        ggplot2::geom_point(size = 2) +
        ggplot2::scale_color_brewer(palette = "Set1", type = "qual") +
        ggplot2::theme_bw())
      dev.off()
    } else if (is.null(type) && !is.null(batch)) {
      pdf(filename, ...)
      print(ggplot2::ggplot(pca_df, ggplot2::aes_string(x = pca_df[,1], y = pca_df[,2], col = groups)) +
        ggplot2::labs(title = "PCA after normalization", x = paste0("PC ", components[1], " (", importance[1], "%)"),
             y = paste0("PC ", components[2], " (", importance[2], "%)"), col = "Group") +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_text(ggplot2::aes(label = batch, hjust = 0, vjust = 1.5), show.legend = F) +
        ggplot2::scale_color_brewer(palette = "Set1", type = "qual") +
        ggplot2::theme_bw())
      dev.off()
    } else {
      pdf(filename, ...)
      print(ggplot2::ggplot(pca_df, ggplot2::aes_string(x = pca_df[,1], y = pca_df[,2], col = groups, shape = type)) +
        ggplot2::labs(title = "PCA after normalization", x = paste0("PC ", components[1], " (", importance[1], "%)"),
             y = paste0("PC ", components[2], " (", importance[2], "%)"), col = "Group", shape = "Type") +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_text(ggplot2::aes(label = batch, hjust = 0, vjust = 1.5), show.legend = F) +
        ggplot2::scale_shape_manual(values = c(19,17)) +
        ggplot2::scale_color_brewer(palette = "Set1", type = "qual") +
        ggplot2::theme_bw())
      dev.off()
    }
  } else {
    if(is.null(type) && is.null(batch)) {
      ggplot2::ggplot(pca_df, ggplot2::aes_string(x = pca_df[,1], y = pca_df[,2], col = groups)) +
        ggplot2::labs(title = "PCA after normalization", x = paste0("PC ", components[1], " (", importance[1], "%)"),
                      y = paste0("PC ", components[2], " (", importance[2], "%)"), col = "Group") +
        ggplot2::geom_point(size = 2) +
        ggplot2::scale_color_brewer(palette = "Set1", type = "qual") +
        ggplot2::theme_bw()
    } else if (!is.null(type) && is.null(batch)) {
      ggplot2::ggplot(pca_df, ggplot2::aes_string(x = pca_df[,1], y = pca_df[,2], col = groups, shape = type)) +
        ggplot2::labs(title = "PCA after normalization", x = paste0("PC ", components[1]," (", importance[1], "%)"),
                      y = paste0("PC ", components[2], " (", importance[2], "%)"), col = "Group", shape = "Type") +
        ggplot2::geom_point(size = 2) +
        ggplot2::scale_color_brewer(palette = "Set1", type = "qual") +
        ggplot2::theme_bw()
    } else if (is.null(type) && !is.null(batch)) {
      ggplot2::ggplot(pca_df, ggplot2::aes_string(x = pca_df[,1], y = pca_df[,2], col = groups)) +
        ggplot2::labs(title = "PCA after normalization", x = paste0("PC ", components[1], " (", importance[1], "%)"),
                      y = paste0("PC ", components[2], " (", importance[2], "%)"), col = "Group") +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_text(ggplot2::aes(label = batch, hjust = 0, vjust = 1.5), show.legend = F) +
        ggplot2::scale_color_brewer(palette = "Set1", type = "qual") +
        ggplot2::theme_bw()
    } else {
      ggplot2::ggplot(pca_df, ggplot2::aes_string(x = pca_df[,1], y = pca_df[,2], col = groups, shape = type)) +
        ggplot2::labs(title = "PCA after normalization", x = paste0("PC ", components[1], " (", importance[1], "%)"),
                      y = paste0("PC ", components[2], " (", importance[2], "%)"), col = "Group", shape = "Type") +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_text(ggplot2::aes(label = batch, hjust = 0, vjust = 1.5), show.legend = F) +
        ggplot2::scale_shape_manual(values = c(19,17)) +
        ggplot2::scale_color_brewer(palette = "Set1", type = "qual") +
        ggplot2::theme_bw()

    }
  }
}
