#' Set feature data
#'
#' Set feature data into ExpressionSet object. Additionaly, it is possible to remove promiscuous probes that
#' map to different ids.
#'
#' @param eset ExpressionSet object. Usually the result of normalization process
#' @param fdata Feature data (dataframe) containing at least two columns: probe id and other id (gene symbol, for instance).
#' @param probe_col Column name for probes in \code{fdata}
#' @param rm_probes Logical. Whether toremove promiscuous probes form feature data dataframe.
#' @param ref_col Column name of feature data with the reference id to look for promiscuous probes.
#'
#' @return ExpressionSet object with filtered assayData slot.
#' @export
#'
#' @examples
set_fdata <- function(eset, fdata, probe_col, rm_probes = FALSE, ref_col = "hgnc_symbol") {

  fdata <- unique(na.omit(fdata))
  id1 <- dplyr::sym(probe_col)
  id2 <- dplyr::sym(ref_col)

  if(rm_probes) {
    if(is.null(ref_col)) {
      stop("Provide a reference column in order to remove probes that map to more than one id.")
    }

    # Remove promiscuous probes
    rem_probes <- fdata %>%
      dplyr::group_by((!!id1)) %>%
      dplyr::summarise(c = dplyr::n_distinct((!!id2))) %>%
      dplyr::ungroup() %>%
      dplyr::filter(c > 1) %>%
      dplyr::select((!!id1)) %>%
      unlist(use.names = F)


    # Subset fdata by probes that are not promiscuous
    probes_in <- fdata[,probe_col] %in% rem_probes
    fdata <- fdata[!probes_in,]

    # Filter fdata containing not-duplicated values
    dup <- sum(duplicated(fdata[,probe_col]))
    fdata <- fdata %>%
      dplyr::filter(!duplicated((!!id1)))

    # Set probes as rownames
    rownames(fdata) <- fdata[,probe_col]
    fdata <- fdata[, !(colnames(fdata) %in% probe_col)]

    # Calculate proportion of promiscuous probes removed from feature data
    prop <- round(sum(probes_in) / length(unique(rownames(fdata))) * 100, digits = 2)

    # Filter assayData by feature data and assign fData
    eset <- eset[rownames(eset) %in% rownames(fdata), ]
    eset <- eset[base::match(rownames(fdata), rownames(eset)), ]
    Biobase::fData(eset) <- fdata
    message(prop, "% of probes mapping to different ", ref_col, " were removed. Also, ", dup, " duplicated probes in feature data were removed.")
    return(eset)

  } else {

    # Filter fdata containing not-duplicated values
    dup <- sum(duplicated(fdata[,probe_col]))
    fdata <- fdata %>%
      dplyr::filter(!duplicated((!!id1)))

    # Set probes as rownames
    rownames(fdata) <- fdata[,probe_col]
    fdata <- fdata[, !(colnames(fdata) %in% probe_col)]

    # Filter assayData by feature data and assign fData
    eset <- eset[rownames(eset) %in% rownames(fdata), ]
    eset <- eset[base::match(rownames(fdata), rownames(eset)), ]
    Biobase::fData(eset) <- fdata
    message(dup, " duplicated probes in feature data were removed.")
    return(eset)
  }

}



