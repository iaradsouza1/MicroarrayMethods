#' Plot PCs proportions
#'
#' @param eset ExpressionSet object with normalized expression values. Commonly the result of RMA normalization
#' @param filename File name to save plot (pdf format)
#' @param ... Other arguments to be passed to pdf() function
#'
#' @return PCs proportions plot
#' @export
#'
#' @examples
pca_prop <- function(eset, filename=NULL, ...) {

  data_exprs <- as.data.frame(na.omit(Biobase::exprs(eset)))
  tdata_exprs <- t(data_exprs)
  importance <- data.frame(t(summary(prcomp(tdata_exprs))$importance), stringsAsFactors = F)
  colnames(importance) <- c("sd", "var_prop", "cum_prop")
  importance$pc <- factor(rownames(importance), levels = rownames(importance))
  importance[,2:3] <- sapply(importance[,2:3], function(x) {x * 100})

  if(!is.null(filename)) {
    pdf(file = filename, ...)
    print(ggplot2::ggplot(importance, ggplot2::aes(x = pc, y = var_prop)) +
      ggplot2::geom_bar(stat = "identity", fill = "coral") +
      ggplot2::geom_line(ggplot2::aes(y = cum_prop, group = 1), lty = 2, show.legend = T) +
      ggplot2::scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
      ggplot2::labs(title = "Principal components proportions", x = "PC", y = "Variation proportion (%)") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(vjust = 0.5, angle = 90)))
    dev.off()
  } else {
    ggplot2::ggplot(importance, ggplot2::aes(x = pc, y = var_prop)) +
      ggplot2::geom_bar(stat = "identity", fill = "coral") +
      ggplot2::geom_line(ggplot2::aes(y = cum_prop, group = 1), lty = 2, show.legend = T) +
      ggplot2::scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
      ggplot2::labs(title = "Principal components proportions", x = "PC", y = "Variation proportion (%)") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(vjust = 0.5, angle = 90))
  }

}
