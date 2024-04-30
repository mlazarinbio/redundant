

#' Output from list to dataframe
#'
#' @param data output list object from redFD
#'
#' @return A named dataframe with metrics of FD as columns. "Spnr" is the number of species that was eventually removed from the community.
#'
#' @example
#' ## load example data from package
#' data(cmnt) ## community matrix
#' data(traits) ## traits
#'
#' ## remove sequentially one and two species (using two random samples) and calculate the FD indices
#' \dontrun{ rem <- redFD(cmnt = cmnt, trait = traits, spn = c(1,2), rep = 3, calc.CWM = FALSE, m = 4) }
#' ## output list object to dataframe
#' \dontrun{ out <- out_dtfr(rem)}
#'
#' @export
#'
out_dtfr <- function(data) {

  while (length(grep("sp",names(data))) == 0) {
  data <- unlist(data, recursive = FALSE)
  }
  indx <- grep("sp",names(data))
  data <- data[indx]
  data <- unlist(data, recursive = FALSE)
  ch <- sapply(data, class)
  if ( length(grep("data.frame",ch)) != 0 | length(grep("matrix",ch)) != 0) {stop("Unable to produce dataframe")}
  x <- length(data)/length(indx)
  res <-  data.frame(do.call(mapply, c(cbind, split(data, cut(seq_along(data), length(data)/x, labels = FALSE)))))
  resnames <- sapply(seq_along(res), function(i) {n <- unlist(strsplit(names(data[i]),"[.]")) ; n[length(n)]})
  colnames(res) <- resnames
  return(res)
}


