#' Remove probes
#'
#' Remove promiscuous probes that map to different reference ids (usually gene symbols).
#'
#' @param fdata Feature data (dataframe) with unique rows. The output of get_annotation() function
#' @param probe_col Probe column name in fdata dataframe
#' @param ref_col Column name of the reference id to be counted for each probe.
#'
#' @return ExpressionSet object with filtered assayData slot for promiscous probes
#' @export
#'
#' @examples
remove_probes <- function(fdata, probe_col, ref_col = "hgnc_symbol") {

  fdata <- unique(na.omit(fdata))
  id1 <- dplyr::sym(probe_col)
  id2 <- dplyr::sym(ref_col)

  # Remove promiscuous probes
  id1 <- dplyr::sym(probe_col)
  id2 <- dplyr::sym(ref_col)
  rem_probes <- fdata %>%
    dplyr::group_by((!!id1)) %>%
    dplyr::summarise(c = dplyr::n_distinct((!!id2))) %>%
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

  message(dup, " duplicated probes in feature data were removed.")

  return(fdata)
}



