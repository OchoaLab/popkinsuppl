# shared with Hudson Fst
clean_labs <- function(labs) {
    # compute some of the WC parameters related to the populations overall (SNP-independent)
    # array of unique populations
    labs_k <- sort(unique(labs))
    # number of populations
    r <- length(labs_k)
    # some mappings
    k2is <- vector('list',r)
    k2n <- vector('integer',r)
    for (k in 1:r) {
        is <- which(labs == labs_k[k])
        k2is[[k]] <- is # indexes with this population
        k2n[k] <- length(is) # size of population
    }
    # these are the things we need to return
    # note that r is the length of both things, so we don't return it separately
    list(
        k2is = k2is,
        k2n = k2n
    )
}
