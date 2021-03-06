#' Differential expression analysis
#'
#' Differential expression analysis based on limma methods.
#'
#' @param eset ExpressionSet object. Usually the result of the normalization process
#' @param contrasts A vector containing contrasts in which differential expression will be done
#' @param method Character string specifying how genes and contrasts are to be combined in the multiple testing scheme.
#' See \code{?decideTests}
#' @param adjust.method Character string specifying p-value adjustment method. See \code{?decideTests}
#' @param p.value numeric value between 0 and 1 giving the required family-wise error rate or false discovery rate.
#' See \code{?decideTests}
#' @param degenes_only Logical value to specify whether return a dataframe with only differentially expressed genes
#' @param ... Other parameters to be passed to \code{decideTests}. See \code{?decideTests}
#'
#' @return A dataframe
#' @export
#'
#' @examples
diff_exp <- function(eset, contrasts, method, adjust.method, p.value = 0.05, degenes_only = FALSE, ...) {

  temp <- as.factor(Biobase::pData(eset)[, group])
  design <- model.matrix(~0 + temp)
  cols <- colnames(design)
  colnames(design) <- gsub("temp", "", cols)
  fit <- limma::lmFit(eset, design)
  contrasts <- limma::makeContrasts(contrasts = contrasts, levels = design)
  ct.fit <- limma::eBayes(limma::contrasts.fit(fit, contrasts))
  res.fit <- limma::decideTests(ct.fit, method = method, adjust.method = adjust.method, p.value = p.value, ...)
  sh.limma <- data.frame(ct.fit$genes, logFC = ct.fit$coef, p.value = ct.fit$p.value,
                         degenes = unclass(res.fit), stringsAsFactors = FALSE)

  if (degenes_only) {
    features <- rowSums(res.fit != 0) > 0
    features <- names(features)[features]
    de.limma <- sh.limma[features, ]
    de.limma <- de.limma[complete.cases(de.limma), ]
    return(de.limma)
  } else {
    return(sh.limma)
  }
}
