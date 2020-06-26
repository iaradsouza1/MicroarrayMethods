#' Create TG dictionary
#'
#' Create a dictionary to transcriptogramer package
#'
#' @param ensembleGeneId vector of ensemble GeneId that have counts.
#' @param dataset ENSEMBL dataset name. A dataset name can be retrieved with list_datasets() function.
#'
#' @details
#'
#' @return
#' @export
#'
#' @examples
create_tg_dict <- function(ensembleGeneId, dataset = "hsapiens_gene_ensembl") {
  if (missing("ensembleGeneId")){
    stop("A ensemble GeneId list must be provided\n")
  }
  att <- c("ensembl_gene_id", "ensembl_peptide_id")
  platform<-"none"
  db <- "main"
  dic_tg <- data.frame(PROBE=ensembleGeneId, stringsAsFactors = F)
  bmList <- get_annotation(dataset = "hsapiens_gene_ensembl",att = c("ensembl_gene_id", "ensembl_peptide_id"))
  #Discard any empty reference of ensembl_peptide_id
  bmList<-bmList[grep("ENSP",bmList[,"ensembl_peptide_id"]),]
  #Dictionary with gene ID (PROBE) and ensemble ID (ENSP)
  dic_tg<-merge(dic_tg,bmList,
                               by.x="PROBE",by.y = "ensembl_gene_id")

  colnames(dic_tg)<-c("PROBE","ENSP")

  #the prefix 9606. will be needed later
  dic_tg$ENSP<-gsub("ENSP","9606.ENSP",
                                   dic_tg$ENSP,fixed = T)
  return(dic_tg)

}