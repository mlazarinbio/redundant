


#' Calculate functional diversity indices for a community after the removal of a specific number of species according to their rarity (based on FD package)
#'
#' @param cmnt matrix containing the abundances of the species in x (or presence-absence, i.e. 0 or 1). Rows are sites and species are columns. Can be missing, in which case dbFD assumes that there is only one community with equal abundances of all species. NAs will be replaced by 0. The number of species (columns) in a must match the number of species (rows) in x. In addition, the species labels in a and x must be identical and in the same order.
#' @param trait matrix or data frame of functional traits. Traits can be numeric, ordered, or factor. Binary traits should be numeric and only contain 0 and 1. character traits will be converted to factor. NAs are tolerated.x can also be a species-by-species distance matrix of class dist, in which case NAs are not allowed. When there is only one trait, x can be also be a numeric vector, an ordered factor, or a unordered factor. In all cases, species labels are required.
#' @param spn positive integer or vector of positive integers specifying the number or numbers of species to be removed from the community. If spn = 1, all species are sequentially removed. If spn >1, random samples (without replacement) of species equal to spn are removed.
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
#' \dontrun{ rem <- redFD(cmnt = cmnt[,1:10], trait = traits[1:10,], spn = c(1,2), rep = 3) }


redFD_rare <- function(cmnt, trait, spn, cmn_to_rare = FALSE, ...) {

  if (is.matrix(cmnt)) stop("Community is not of class matrix")
  if (missing(spn)) { spn <- max(rowSums(cmnt)) - 3
  warning(paste("Missing spn, setting to maximum possible ->",spn))}
  if (max(spn) > max(rowSums(cmnt)) - 3) {
    stop(paste0("Species number set to total species number"))
  }

  if (length(spn) > 1) {
    FD.ALL <- lapply(1:length(spn), function(x) rare_removal(spn[x], ...))
  } else {FD.ALL <- rare_removal(spn, ...)}

  return(FD.ALL)
}

rare_removal <-  function(x) eval.parent(substitute( {

  sp <- colnames(cmnt)
  Nsp <- apply(cmnt, 2, sum)
  rank <- data.frame(sp = sp, Nr = Nsp, stringsAsFactors = FALSE)
  if (cmn_to_rare == FALSE) { rem <- rank$Nr }
  else { rem <- -rank$Nr }
  rank <- rank[order(rem), ]

  if (x == 1){

    indx <- lapply(seq_along(sp), function(i) rank[i, 1])
    tabc <- lapply(seq_along(indx), function(i) cmnt[, !(colnames(cmnt) %in% indx[[i]])])
    tabc <- lapply(seq_along(tabc), function(i) tabc[[i]][which(rowSums(tabc[[i]]) >= 3), , drop = FALSE])
    tabc <- tabc[which(sapply(tabc, function(i) length(tabc[i])) > 0)]

    if ( length(tabc) == 0) {warning(paste0("No species for calculations at spn ", x)) ; next}

    test <- lapply(seq_along(tabc), function(i) colSums(tabc[[i]][1:nrow(tabc[[i]]), 1:ncol(tabc[[i]]), drop = FALSE]))
    test <- unlist(lapply(seq_along(test), function(i) test[[i]] == 0))

    if (length(test == TRUE) > 0){

      warning("At least one species had zero total abundance, removing that species")
      try( tabc <- lapply(seq_along(tabc), function(i) tabc[[i]][ , which(colSums(tabc[[i]]) != 0), drop = FALSE]))
    }

    tabc <- tabc[which(sapply(tabc, function(i) length(tabc[i])) > 0)]

    if ( length(tabc) == 0) {warning(paste0("No species for calculations at spn ", x)) ; next}

    ind <- lapply(seq_along(tabc), function(i) colnames(tabc[[i]]))
    tabt <- lapply(seq_along(ind), function(i) trait[row.names(trait) %in% ind[[i]], ])

    n <- sapply(seq_along(tabc), function(i) setdiff(sp, colnames(tabc[[i]])), simplify = FALSE)
    names <- sapply(seq_along(n), FUN =function(i) paste(n[[i]], collapse = "/"))
    FD.temp <- lapply(seq_along(tabc), function(i) {
      db <- dbFD(tabt[[i]], tabc[[i]], calc.CWM = FALSE)})
  }
  else{

    #indx <- lapply(seq_len(combi), function(i) rank[1:x])
    indx <- rank[1:x, 1]
    tabc <- cmnt[, !(colnames(cmnt) %in% indx)]
    tabc <- tabc[which(rowSums(tabc) >= 3), , drop = FALSE]


    if (nrow(tabc) == 0) {warning(paste0("No species for calculations at spn ", x)) ; next}
    test <- colSums(tabc)
    test <- test[test == 0]

    if (length(test == TRUE) > 0){

      warning("At least one species had zero total abundance, removing that species")
      try( tabc <- tabc[ , which(colSums(tabc) != 0)])
    }

    tabc <- tabc[which(rowSums(tabc) >= 3), , drop = FALSE]

    if ( nrow(tabc) == 0) {warning(paste0("No species for calculations at spn ", x)) ; next}

    ind <- colnames(tabc)
    tabt <- trait[row.names(trait) %in% ind, ]

    n <- setdiff(sp, colnames(tabc))
    names <- paste(n, collapse = "/")
    db <- dbFD(tabt, tabc, calc.CWM = FALSE)
    FD.temp <- list(db)
  }


  names(FD.temp) <-  names
  return(FD.temp)

}))


