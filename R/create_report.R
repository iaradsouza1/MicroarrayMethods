#' Create report
#'
#' Create quality assessment report of Affymetrix CEL files
#'
#' @param celfiles Path to the celfiles directory
#' @param pheno_data File name with metadata information from the experiment
#' @param sep \code{pheno_data} file separated format
#' @param filenames A vector containg the names of .CEL files, if only a subset of files are needed to be imported.
#' @param components A vector of length 2 containing which principal components to plot
#' @param groups Column name corresponding to the groups of samples annotation on pheno data provided to the \code{import_celfiles()} function
#' @param type Column name corresponding to the types of samples (control and treated samples) on pheno data provided to the \code{import_celfiles()} function (optional)
#' @param batch Column name corresponding to each replicates on pheno data provided to the \code{import_celfiles()} function (optional)
#' @param ... Other arguments to be passed to \code{knitr::render()} function. Useful arguments:
#' \code{output_file}, \code{output_dir}, etc. See \code{?render}
#'
#' @details \code{pheno_data} is a table containing basic information about the experiment, such as GSMid and sample groups.
#' It must have at least these two columns (GSMid and sample groups). For PCA plot, one could provide the type of each samples (control or treated)
#' and the batches (replicates) of each sample on different columns.
#'
#' @return Produces an html file with quality assessment plots.
#' @export
#'
#' @examples
create_report <- function(celfiles = NULL, pheno_data = NULL, sep = NULL, filenames = NULL,
                                        components = NULL, group = NULL, batch = NULL, type = NULL, ...) {

  if (is.null(celfiles) || is.null(pheno_data) || is.null(sep) || is.null(group) || is.null(components)) {
    stop("All of the following arguments must be provided:
CEL files path (celfiles), pheno data file name (pheno_data), pheno data text file separator (sep), PCA components (components),
         and column name for group annotation on pheno data (group)")
  }

  qc_path <- system.file("rmd", "qc.Rmd", package = "MicroarrayMethods")
  suppressWarnings(rmarkdown::render(qc_path,
                    params = list(
                      celfiles = celfiles,
                      pheno_data = pheno_data,
                      sep = sep,
                      filenames = filenames,
                      components = components,
                      group = group,
                      batch = batch,
                      type = type),
                    output_format = "html_document", quiet = T, ...))

}


