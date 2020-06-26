#' Compare variances
#'
#' Compare the variance of each group with the total variance
#'
#' @param var_groups A dataframe with the variances of each group
#' @param var_total A numeric vector with the total variance
#'
#' @return A named vector with the proportion of probes from groups with greater variance
#' than the total variance.
#' @export
#'
#' @examples
compare_variances <- function(var_groups, var_total) {
  res <- lapply(var_groups, function(x) {
    x > var_total
  })
  sapply(res, mean)
}
