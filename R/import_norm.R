#' Import and normalize
#'
#' Import CEL files (Affymetrix) and normalize with RMA method
#'
#' @param celfiles_path Path to the celfiles directory
#' @param pheno_data File name with metadata information from the experiment
#' @param sep \code{pheno_data} separated format
#' @param filenames a vector containg the names of .CEL files, if only a subset of files are needed to be imported. (optional)
#' @param ... Other parameters to be passed to \code{rma()} function. See \code{?rma}
#'
#' @return An ExpressionSet object with normalized expression values
#' @export
#'
#' @examples
import_norm <- function(celfiles_path, pheno_data, sep, filenames = NULL, ...) {

  raw <- import_celfiles(celfiles_path = celfiles_path, pheno_data = pheno_data, sep = sep, filenames = filenames)
  eset <- affy::rma(raw, ...)
  return(eset)

}
