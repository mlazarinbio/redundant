

#' Calculate functional diversity indices for a community after the removal of a specific number of species (based on FD package)
#'
#' @param cmnt matrix containing the abundances of the species in x (or presence-absence, i.e. 0 or 1). Rows are sites and species are columns. Can be missing, in which case dbFD assumes that there is only one community with equal abundances of all species. NAs will be replaced by 0. The number of species (columns) in a must match the number of species (rows) in x. In addition, the species labels in a and x must be identical and in the same order.
#' @param trait matrix or data frame of functional traits. Traits can be numeric, ordered, or factor. Binary traits should be numeric and only contain 0 and 1. character traits will be converted to factor. NAs are tolerated.x can also be a species-by-species distance matrix of class dist, in which case NAs are not allowed. When there is only one trait, x can be also be a numeric vector, an ordered factor, or a unordered factor. In all cases, species labels are required.
#' @param spn positive integer or vector of positive integers specifying the number or numbers of species to be removed from the community. If spn = 1, all species are sequentially removed. If spn >1, random samples (without replacement) of species equal to spn are removed.
#' @param rep positive integer specifying the number of repetitions to remove spn random species from the community.
#' @param ... arguments to be passed to \code{\link[FD]{dbFD}}.
#'
#' @return A named list with the indices calculated. The names of each element are the species names that were removed from the community. More species than initially specified, may be removed if they had no appearance in the community. Rows with less than 3 species are removed from the result with a warning. If a community has 0 rows, it is removed with a warning.
#'
#' @example
#' ## load example data from package
#' data(cmnt) ## community matrix
#' data(traits) ## traits
#'
#' ## remove sequentially one and two species (using two random samples) and calculate the FD indices
#' \dontrun{ rem <- redFD(cmnt = cmnt, trait = traits, spn = c(1,2), rep = 3, m = 4, calc.CWM = FALSE)}
#'
#' @export



redFD <- function(cmnt, trait, spn, rep, ...) {

  if (class(cmnt) != "matrix") stop("Community is not of class matrix")
  if (missing(spn)) { spn <- max(rowSums(cmnt)) - 3
  warning(paste("Missing spn, setting to maximum possible ->",spn))}
  if (max(spn) > max(rowSums(cmnt)) - 3) {
    stop(paste0("Species number set to total species number"))
  }
  if (missing(rep)) { rep <- 10 }
  if (length(spn) > 1) {
    FD.ALL <- lapply(1:length(spn), function(x) basic_removal(spn[x], ...))
  } else {FD.ALL <- basic_removal(spn, ...)}

  return(FD.ALL)
}




basic_removal <- function(x,...) eval.parent(substitute( {

  sp <- colnames(cmnt)
  cmntr <- nrow(cmnt)
  cmntc <- ncol(cmnt)

  if (x == 1) {

    indx <- lapply(seq_along(sp), function(i) sp[i])

  }
  else{

    indx <- lapply(seq_len(rep), function(i) sample(sp, x, replace = FALSE))

  }
  tabc <- lapply(seq_along(indx), function(i) cmnt[, !(sp %in% indx[[i]]), drop = FALSE])
  rmrow <- lapply(seq_along(tabc), function(i) which(rowSums(tabc[[i]]) >= 3))
  tabc <- lapply(seq_along(tabc), function(i) tabc[[i]][which(rowSums(tabc[[i]]) >= 3), , drop = FALSE])
  tabc <- tabc[which(sapply(seq_along(tabc), function(i) length(tabc[i])) > 0)]
  if ( length(tabc) == 0) {warning(paste0("No species for calculations at spn ", x)) ; next}

  test <- lapply(seq_along(tabc), function(i) colSums(tabc[[i]]))
  test <- unlist(lapply(seq_along(test), function(i) test[[i]] == 0))

  if (length(which(test == TRUE)) > 0) {

    warning("At least one species was absent in all cells, removing that species")
    try( tabc <- lapply(seq_along(tabc), function(i) tabc[[i]][ , which(colSums(tabc[[i]]) != 0), drop = FALSE]))
  }

  tabc <- tabc[which(sapply(seq_along(tabc), function(i) length(tabc[i])) > 0)]

  if ( length(tabc) == 0) {warning(paste0("No species for calculations at spn ", x)) ; next}

  ind <- lapply(seq_along(tabc), function(i) colnames(tabc[[i]]))
  tabt <- lapply(seq_along(ind), function(i) trait[row.names(trait) %in% ind[[i]], , drop = FALSE])

  n <- sapply(seq_along(tabc), function(i) setdiff(sp, ind[[i]]), simplify = FALSE)
  names <- sapply(seq_along(n), FUN = function(i) paste(n[[i]], collapse = "/"))
  FDtab <- lapply(seq_along(tabc), function(i) {
    print(i)
    db <- tryCatch(dbFD(x = tabt[[i]], a = tabc[[i]], ...), error = function(e) NULL)
    sdb <- c(db, id = list(rmrow[[i]]))
    })

  names(FDtab) <-  names
  return(FDtab)

}))
