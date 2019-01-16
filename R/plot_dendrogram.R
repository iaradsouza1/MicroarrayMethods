#' Samples dendrogram
#'
#' @param eset ExpressionSet object with normalized expression values. Commonly the result of RMA normalization
#' @param groups Column name corresponding to the groups of samples annotation on pheno data provided to the `import_celfiles` function
#' @param method Method to calculate distance between samples. See ?dist()
#' @param filename File name to save plot (pdf format)
#' @param ... Other arguments to be passed to pdf() function
#'
#' @return
#' @export
#'
#' @examples
plot_dendrogram <- function(eset, groups, method = "euclidean", filename=NULL, ...) {
  data_exprs <- as.data.frame(na.omit(exprs(eset)))
  tdata_exprs <- t(data_exprs)
  rownames(tdata_exprs) <- pData(eset)[, groups]
  d <- dist(as.matrix(tdata_exprs), method = method)
  clusters <- hclust(d)

  if (!is.null(filename)) {
    pdf(file = filename, ...)
    print(plot(clusters))
    dev.off()
  } else {
    plot(clusters)
  }



}
