#' List all attributes
#'
#' @param dataset ENSEMBL dataset name. A dataset name can be retrieved with list_datasets() function.
#' @param db One of the following ENSEMBL databases: "main" for ENSEMBL_MART_ENSEMBL database (the main database in ENSEMBL, provides information for the majority of species)
#' and "mouse" for ENSEMBL_MART_MOUSE database (which provides specific information of mouse strains)
#'
#' @return Dataframe containing all attributes from a dataset.
#' @export
#'
#' @examples
list_all_attributes <- function(dataset, db = "main") {
  if (db == "main") {
    ensembl <- biomaRt::useDataset(dataset = dataset, mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL"))
    att <- biomaRt::listAttributes(ensembl)
    return(att)
  } else if (db == "mouse") {
    ensembl <- biomaRt::useDataset(dataset = dataset, mart = biomaRt::useMart("ENSEMBL_MART_MOUSE"))
    att <- biomaRt::listAttributes(ensembl)
    return(att)
  } else {
    stop("Provide an ensembl value: \n 'main': corresponding to ENSEMBL_MART_ENSEMBL database or 'mouse': corresponding to ENSEMBL_MART_MOUSE database")
  }
}
