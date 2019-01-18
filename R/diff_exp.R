#' Title
#'
#' @param eset
#' @param contrasts
#' @param method
#' @param adjust.method
#' @param p.value
#' @param degenes_only
#' @param ...
#'
#' @return
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
  res.fit <- limma::decideTests(ct.fit, method = method, adjust.method = adjust.method, p.value = p.value)
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