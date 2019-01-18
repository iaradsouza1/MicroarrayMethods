#' Import celfiles (Affymetrix)
#'
#' @param celfiles_path Path to the celfiles directory
#' @param pheno_data File name with metadata information from the experiment
#' @param sep \code{pheno_data} separated format
#' @param filenames a vector containg the names of .CEL files, if only a subset of files are needed to be imported. (optional)
#'
#' @details \code{pheno_data} is a table containing basic information about the experiment, such as GSMid and sample groups.
#' It must have at least these two columns (GSMid and sample groups). For PCA plot, one could provide the type of each samples (control or treated)
#' and the batches (replicates) of each sample.
#'
#'
#' @return An AffyBatch object with imported celfiles
#' @export
#'
#' @examples
import_celfiles <- function(celfiles_path, pheno_data, sep, filenames=NULL) {

  # Import only the specified cel files
  if (!is.null(filenames)) {
    previous_path <- getwd()
    setwd(celfiles_path)
    pheno_data <- read.table(file = paste0(celfiles_path, pheno_data), sep = sep,
                             header = T, row.names = list_celfiles(path = celfiles_path),
                             stringsAsFactors = F)
    pheno_data <- pheno_data[filenames,]
    raw <- affy::ReadAffy(filenames = filenames, phenoData = pheno_data)
    setwd(previous_path)
    return(raw)

  # Import all cel files
  } else {
    pheno_data <- read.table(file = paste0(celfiles_path, pheno_data), sep = sep,
                             header = T, row.names = list_celfiles(path = celfiles_path),
                             stringsAsFactors = F)

    raw <- affy::ReadAffy(celfile.path = celfiles_path, phenoData = pheno_data)
    return(raw)
  }

}
