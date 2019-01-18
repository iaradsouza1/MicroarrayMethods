#' List celfiles in a given directory
#'
#' @param path Directory path to look for .CEL files
#'
#' @return A character vector correponding to cel file names in a given directory.
#' @export
#'
#' @examples
list_celfiles <- function(path) {
  files <- dir(path)
  files <- files[grep("\\.cel|\\.CEL", files)]
  return(files)
}
