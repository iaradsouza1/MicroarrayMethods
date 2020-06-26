#' Get annotation
#'
#' Get annotation for a specific platform
#'
#' @param dataset ENSEMBL dataset name. A dataset name can be retrieved with list_datasets() function.
#' @param platform Platform name in ENSEMBL. Usually one retrieved by list_platforms() function. Use "none" if you aren't working with microarray. In that case the parameter "att" is mandatory.
#' @param db One of the following ENSEMBL databases: "main" for ENSEMBL_MART_ENSEMBL database (the main database in ENSEMBL, provides information for the majority of species)
#' and "mouse" for ENSEMBL_MART_MOUSE database (which provides specific information of mouse strains)
#' @param att Attributes to be retrieved. If not specified, retrives HGNC symbol and ENTREZ id.
#'
#' @details The \code{att} argument expects a vector containing the attributes according to biomart package. One can
#' retrieve all available attributes with list_all_attributes() function.
#'
#' @return Dataframe containing the attributes requested.
#' @export
#'
#' @examples
get_annotation <- function(dataset, platform="none", db = "main", att ) {
  if (missing("dataset")){
    stop("A dataset object must be provided\n")
  }
  if (platform == "none"){
    if (missing("att")){
      stop("An vector of attributes or a plataform must be provided\n")
    }
  }else {
    att = c(platform, "hgnc_symbol", "entrezgene")#replaces line att <- c(platform, att)
  }
  if (db == "main") {
    ensembl <- biomaRt::useDataset(dataset = dataset, mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL"))
    #att <- c(platform, att)
    id <- biomaRt::getBM(attributes = att, mart = ensembl)

    id[id == ""] <- NA
    id <- unique(na.omit(id))
    return(id)
  } else if (db == "mouse") {
    ensembl <- biomaRt::useDataset(dataset = dataset, mart = biomaRt::useMart("ENSEMBL_MART_MOUSE"))
    #att <- c(platform, att)
    id <- biomaRt::getBM(attributes = att, mart = ensembl)
    id[id == ""] <- NA
    id <- unique(na.omit(id))
    return(id)
  } else {
    stop("Provide an ensembl value: \n 'main': corresponding to ENSEMBL_MART_ENSEMBL database or 'mouse': corresponding to ENSEMBL_MART_MOUSE database")
  }
}
