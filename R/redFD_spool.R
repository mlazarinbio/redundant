

#' Calculate functional diversity indices for a community after the removal of a specific number of species taking into account the regional species pools (based on FD package)
#'
#' @param cmnt matrix containing the abundances of the species in x (or presence-absence, i.e. 0 or 1). Rows are sites and species are columns. Can be missing, in which case dbFD assumes that there is only one community with equal abundances of all species. NAs will be replaced by 0. The number of species (columns) in a must match the number of species (rows) in x. In addition, the species labels in a and x must be identical and in the same order.
#' @param trait matrix or data frame of functional traits. Traits can be numeric, ordered, or factor. Binary traits should be numeric and only contain 0 and 1. character traits will be converted to factor. NAs are tolerated.x can also be a species-by-species distance matrix of class dist, in which case NAs are not allowed. When there is only one trait, x can be also be a numeric vector, an ordered factor, or a unordered factor. In all cases, species labels are required.
#' @param spool matrix that associates each site with a regional species pool. Sites must be ordered as in cmnt matrix. Columns must be named "site" and "region" for the sites and pool columns respectively.
#' @param spn positive integer or vector of positive integers specifying the number or numbers of species to be removed from the community. If spn = 1, all species are sequentially removed. If spn >1, random samples (without replacement) of species equal to spn are removed.
#' @param rep positive integer specifying the number of repetitions to remove spn random species from the community.
#' @param times positive integer specifying number of simulations to perform.
#' @param ... arguments to be passed to \code{\link[FD]{dbFD}}.
#'
#' @return A named list with the indices calculated. The names of each element are the species names that were removed from the community. More species than initially specified, may be removed if they had no appearance in the community. Rows with less than 3 species are removed from the result with a warning. If a community has 0 rows, it is removed with a warning.
#'
#' @example
#' ## load example data from package
#' data(cmnt) ## community matrix
#' data(traits) ## traits
#' data(spool) ## matrix assigning each site to a region
#'
#' \dontrun{ rem <- redFD(cmnt = cmnt, trait = traits, spool = spool, spn = c(1,2), rep = 3, times = 10,  m = 4, calc.CWM = FALSE)}
#'
#' @export

redFD_spool <- function(cmnt, trait, spool, spn, rep, times, ...) {

  #if (identical(cmnt[,0], spool[,"site"]) == FALSE) {stop("Sites in spool not identical to sites in cmnt")}

  sim <- lapply(seq_len(times), function(x) sim_regional())
  FD.sim <- lapply(seq_along(sim), function(i) {
   tmpcmnt <- sim[[i]]
   tmpcmnt <- as.matrix(tmpcmnt[, -ncol(tmpcmnt)])
   tempFD <- redFD(cmnt = tmpcmnt , trait = trait, spn = spn, rep = rep, ...)
   list(metrics = tempFD, cmntsim = sim[i])
  })
  names(FD.sim) <- paste0("simulation", seq_len(times))
  return(FD.sim)

}


sim_regional <- function() eval.parent(substitute({

  cmntc <- ncol(cmnt)
  cmntr <- nrow(cmnt)
  sp <- colnames(cmnt)
  tmp <- data.frame(cmnt, spool[,"region"], stringsAsFactors = FALSE)
  colnames(tmp) <- c(sp, "region")
  region <- unique(tmp[,"region"])
  tmp_region <- aggregate(. ~ region, tmp, sum)
  region_list <- apply(tmp_region[,-1], 1, function(x) colnames(tmp_region[,-1])[which(x > 0)])
  names(region_list) <- tmp_region[[1]]
  tmpmat <- sapply(seq_len(cmntr), function(i)  sample(region_list[[tmp[i,cmntc + 1]]], sum(tmp[i, 1:cmntc]), replace = FALSE))
  sm <- lapply(seq_len(cmntr), function(i) as.numeric(sp %in% tmpmat[[i]]))
  simmat <- matrix(unlist(sm), nrow = cmntr, ncol = cmntc, dimnames = dimnames(cmnt), byrow = TRUE)
  return(simmat)

}))
