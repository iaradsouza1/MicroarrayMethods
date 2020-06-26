#' List ENSEMBL datasets
#'
#' @param db One of the following ENSEMBL databases: "main" for ENSEMBL_MART_ENSEMBL database (the main database in ENSEMBL, provides information for the majority of species)
#' and "mouse" for ENSEMBL_MART_MOUSE database (which provides specific information of mouse strains)
#' @return
#' @export
#'
#' @examples
list_datasets <- function(db = "main") {
  if (db == "main") {
    ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL")
    ds <- biomaRt::listDatasets(ensembl)
    return(ds)
  } else if (db == "mouse") {
    ensembl <- biomaRt::useMart("ENSEMBL_MART_MOUSE")
    ds <- biomaRt::listDatasets(ensembl)
    return(ds)
  } else {
    stop("Provide an ensembl value: \n 'main': corresponding to ENSEMBL_MART_ENSEMBL database or 'mouse': corresponding to ENSEMBL_MART_MOUSE database")
  }
}
