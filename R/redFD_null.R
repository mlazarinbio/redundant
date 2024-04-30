

#' Create simulated communities for null models, remove species and calculate functional diversity indices (based on FD package)
#'
#' @param cmnt matrix containing the abundances of the species in x (or presence-absence, i.e. 0 or 1). Rows are sites and species are columns. Can be missing, in which case dbFD assumes that there is only one community with equal abundances of all species. NAs will be replaced by 0. The number of species (columns) in a must match the number of species (rows) in x. In addition, the species labels in a and x must be identical and in the same order.
#' @param trait matrix or data frame of functional traits. Traits can be numeric, ordered, or factor. Binary traits should be numeric and only contain 0 and 1. character traits will be converted to factor. NAs are tolerated.x can also be a species-by-species distance matrix of class dist, in which case NAs are not allowed. When there is only one trait, x can be also be a numeric vector, an ordered factor, or a unordered factor. In all cases, species labels are required.
#' @param spn positive integer or vector of positive integers specifying the number or numbers of species to be removed from the community. If spn = 1, all species are sequentially removed. If spn >1, random samples (without replacement) of species equal to spn are removed.
#' @param rep positive integer specifying the number of repetitions to remove spn random species from the community.
#' @param permat_args list, arguments to be passed to \code{\link[vegan]{permatfull}}.
#' @param ... arguments to be passed to \code{\link[FD]{dbFD}}.
#'
#' @return A named list with the indices calculated and the simulated community used. The names of each element are the species names that were removed from the community. More species than initially specified, may be removed if they had no appearance in the community. Rows with less than 3 species are removed from the result with a warning. If a community has 0 rows, it is removed with a warning.
#'
#' @example
#' ## load example data from package
#' data(cmnt) ## community matrix
#' data(traits) ## traits
#' \dontrun{ rem1 <- redFD_null(cmnt = cmnt, trait = traits, spn = 2, rep = 3, permat_args = list(times = 4, fixedmar = "rows"), calc.CWM = FALSE, m = 4)}
#'
#' @export
redFD_null <- function(cmnt, trait, spn, rep, permat_args = list(times = 99), ...) {

  temp <- do.call(vegan::permatfull, c(list(m = cmnt), permat_args))
  sim <- temp$perm
  FD.sim <- lapply(seq_along(sim), function(i) {
    temp <- redFD(cmnt = sim[[i]] , trait = trait, spn = spn, rep = rep, ...)
    list(metrics = temp, cmntsim = sim[i])})
  names(FD.sim) <- paste0("simulation", seq_len(permat_args$times))
  return(FD.sim)
}


