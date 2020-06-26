#' Calculate the total variance
#'
#' @param exp_matrix An normalized expression matrix.
#'
#' @return A numeric vector with the total variance for each probe.
#' @export
#'
#' @examples
calc_var_total <- function(exp_matrix) {

  # Calculate total variance
  apply(exp_matrix, 1, var)

}
