#' Plot housekeeping expression levels
#'
#' @param eset An ExpressionSet object
#' @param actin Acting gene id
#' @param gapdh GAPDH gene id
#' @param id Column name for the gene identifier used in eset object
#'
#' @return
#' @export
#'
#' @examples
plot_hk <- function(eset, id, actin, gapdh, filename=NULL, ...) {

  if (!is.null(filename)) {
    pdf(filename, ...)
    exp_gapdh <- exprs(eset)[which(fData(eset)[, id] %in% gapdh)[1],]
    exp_actin <- exprs(eset)[which(fData(eset)[, id] %in% actin)[1],]

    exp_df <- data.frame(samples = names(exp_actin),
                         exp_actin, exp_gapdh, row.names = NULL, stringsAsFactors = F)

    if (grepl(".CEL.gz", exp_df$samples[1])) {
      exp_df$samples <- gsub(".CEL.gz", "", exp_df$samples)
    }

    exp_df$samples <- factor(exp_df$samples, levels = exp_df$samples)
    lim  <- range(c(range(exp_actin), range(exp_gapdh)))
    lim <- lim + c(-0.5, 0.5)

    print(ggplot2::ggplot(exp_df) +
      ggplot2::geom_line(aes(x = samples, y = exp_actin, group = 1, col = "Actin"), lwd = 1) +
      ggplot2::geom_point(aes(x = samples, y = exp_actin, group = 1), col = "blue") +
      ggplot2::geom_line(aes(x = samples, y = exp_gapdh, group = 2, col = "GAPDH"), lwd = 1) +
      ggplot2::geom_point(aes(x = samples, y = exp_gapdh, group = 2), col = "red") +
      ggplot2::scale_y_continuous(limits = lim) +
      ggplot2::scale_color_manual(values = c("Actin" = "blue", "GAPDH" = "red")) +
      ggplot2::labs(title = "Housekeeping genes expression levels",
           x = "Samples", y = "log2-normalized expression level", col = "") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5)))
    dev.off()

  } else {
    exp_gapdh <- exprs(eset)[which(fData(eset)[, id] %in% gapdh)[1],]
    exp_actin <- exprs(eset)[which(fData(eset)[, id] %in% actin)[1],]

    exp_df <- data.frame(samples = names(exp_actin),
                         exp_actin, exp_gapdh, row.names = NULL, stringsAsFactors = F)

    if (grepl(".CEL.gz", exp_df$samples[1])) {
      exp_df$samples <- gsub(".CEL.gz", "", exp_df$samples)
    }

    exp_df$samples <- factor(exp_df$samples, levels = exp_df$samples)
    lim  <- range(c(range(exp_actin), range(exp_gapdh)))
    lim <- lim + c(-0.5, 0.5)

    print(ggplot2::ggplot(exp_df) +
            ggplot2::geom_line(aes(x = samples, y = exp_actin, group = 1, col = "Actin"), lwd = 1) +
            ggplot2::geom_point(aes(x = samples, y = exp_actin, group = 1), col = "blue") +
            ggplot2::geom_line(aes(x = samples, y = exp_gapdh, group = 2, col = "GAPDH"), lwd = 1) +
            ggplot2::geom_point(aes(x = samples, y = exp_gapdh, group = 2), col = "red") +
            ggplot2::scale_y_continuous(limits = lim) +
            ggplot2::scale_color_manual(values = c("Actin" = "blue", "GAPDH" = "red")) +
            ggplot2::labs(title = "Housekeeping genes expression levels",
                          x = "Samples", y = "log2-normalized expression level", col = "") +
            ggplot2::theme_bw() +
            ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5)))
  }

}
