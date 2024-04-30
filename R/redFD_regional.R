
#' Create simulated communities for null models based on regional species pools, remove species and calculate functional diversity indices (based on FD package)
#'
#' @param cmnt matrix containing the abundances of the species in x (or presence-absence, i.e. 0 or 1). Rows are sites and species are columns. Can be missing, in which case dbFD assumes that there is only one community with equal abundances of all species. NAs will be replaced by 0. The number of species (columns) in a must match the number of species (rows) in x. In addition, the species labels in a and x must be identical and in the same order.
#' @param trait matrix or data frame of functional traits. Traits can be numeric, ordered, or factor. Binary traits should be numeric and only contain 0 and 1. character traits will be converted to factor. NAs are tolerated.x can also be a species-by-species distance matrix of class dist, in which case NAs are not allowed. When there is only one trait, x can be also be a numeric vector, an ordered factor, or a unordered factor. In all cases, species labels are required.
#' @param region dataframe of sites and regions. First column must be the sites, as numeric, and second the regions. Sites must be in the same order as cmnt.
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
 redFD_regional <- function(cmnt, trait, region, spn, rep, permat_args = list(times = 99), ...) {

   if (class(region[[1]]) != "numeric") stop("Sites column not numeric")
  temp <- do.call(vegan::permatfull, c(list(m = cmnt), permat_args))
  sim <- temp$perm
  FD.sim <- lapply(seq_along(sim), function(i) {
    temp <- redFD(cmnt = sim[[i]] , trait = trait, spn = spn, rep = rep, ...)
    list(metrics = temp, cmntsim = sim[i])})
  names(FD.sim) <- paste0("simulation", seq_len(permat_args$times))
  return(FD.sim)
 }

 sim_region <- function() eval.parent(substitute({
   require(plyr)

   sp <- colnames(compo)
   Ncom <- nrow(compo)
   Nsp <- apply(compo, 1, sum)
   areas <- sort(unique(region[[2]]))

   tmp <- as.data.frame(compo)
   tmp$names <- paste0("grid", seq_len(Ncom))
   data <- join(tmp, region, type = "inner", by = "names")
   data <- data[,-(ncol(data)-1)]
   sp_pres <- which(data[,1:(ncol(data)-1)] == 1, arr.ind = TRUE, useNames = TRUE)
   sp_pres <- sp_pres[order(sp_pres[,1]),]

   sp_pres <- aggregate(col ~ ., sp_pres, little_f)
   sp_pres <- as.data.frame(sp_pres, stringsAsFactors = FALSE)
   colnames(sp_pres) <- c("names", "Sps")
   tmp <- merge(sp_pres, data[,ncol(data)], by.x = "names", by.y = "row.names")
   colnames(tmp) <- c("names", "Sps", "bio_region")
   tmp <- aggregate(Sps ~ bio_region, tmp, toString )
   Spools <- sapply(seq_along(areas), function(i) unique(unlist(strsplit(tmp[i,2], split = ", "))), simplify = TRUE)
   names(Spools) <- areas

   smpmat <- sapply(seq_along(Nsp), function(i) sample(Spools[[which(names(Spools) == data[i,ncol(data)])]], Nsp[i], replace = FALSE))

   temp <- lapply(seq_len(Ncom), function(i) as.numeric(sp %in% smpmat[[i]]))
   simmat <- matrix(unlist(temp), nrow = Ncom, ncol = ncol(compo), dimnames = dimnames(compo), byrow = TRUE)
   return(simmat)
 }))
