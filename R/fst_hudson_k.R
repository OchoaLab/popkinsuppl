#' The Generalized "Hudson" FST estimator
#'
#' This function implements the "Hudson" FST estimator (Bhatia, Patterson, Sankararaman, and Price 2013) generalized to K subpopulations (Ochoa and Storey 2016).
#' Handles very large data, passed as BEDMatrix or as a regular R matrix.
#' Handles missing values correctly.
#'
#' @param X The genotype matrix (BEDMatrix, regular R matrix, or function, same as `popkin`).
#' @param labs A vector of subpopulation assignments for every individual.
#' @param indexes_keep An optional vector of loci to keep (as booleans or indexes, used to subset an R matrix).
#' @param loci_on_cols Determines the orientation of the genotype matrix (by default, `FALSE`, loci are along the rows).
#' Set automatically to `TRUE` if `X` is a BEDMatrix object.
#'
#' @return A list with the following named elements, in this order:
#' 
#' - fst: The genome-wide Fst estimate (scalar).
#' - fst_loci: A vecctor of per-locus Fst estimates.
#' - data: a 2-by-m matrix of statistics used in estimating Fst.
#'   Useful to obtain a bootstrap distribution for the genome-wide Fst.
#'
#' @examples
#' # dimensions of simulated data
#' n_ind <- 100
#' m_loci <- 1000
#' k_subpops <- 10
#' n_data <- n_ind * m_loci
#' 
#' # missingness rate
#' miss <- 0.1
#' 
#' # simulate ancestral allele frequencies
#' # uniform (0,1)
#' # it'll be ok if some of these are zero
#' p_anc <- runif(m_loci)
#'
#' # simulate some binomial data
#' X <- rbinom(n_data, 2, p_anc)
#' 
#' # sprinkle random missingness
#' X[ sample(X, n_data * miss) ] <- NA
#'
#' # turn into a matrix
#' X <- matrix(X, nrow = m_loci, ncol = n_ind)
#'
#' # create subpopulation labels
#' # k_subpops groups of equal size
#' labs <- ceiling( (1 : n_ind) / k_subpops )
#'
#' # estimate FST using the "Hudson" formula
#' fst_hudson_k_obj <- fst_hudson_k(X, labs)
#' 
#' # the genome-wide FST estimate
#' fst_hudson_k_obj$fst
#' 
#' # vector of per-locus FST estimates
#' fst_hudson_k_obj$fst_loci
#' 
#' @seealso
#' The popkin package.
#'
#' @export
fst_hudson_k <- function(X, labs, indexes_keep = NULL, loci_on_cols = FALSE) {
    
    if (missing(labs))
        stop('Subpopulation labels `labs` are required!')

    # this is always given through labs
    n <- length(labs)

    # shared with WC estimator...
    obj <- preprocess_genotype_obj(X, n, loci_on_cols, indexes_keep = indexes_keep)
    isFn <- obj$isFn
    m <- obj$m
    loci_on_cols <- obj$loci_on_cols

    # global computations shared by all chunks
    
    # general label processing to simplfy things
    out <- clean_labs(labs)
    k2is <- out$k2is
    k2n <- 2 * out$k2n - 1 # NOTE: in Hudson estimator, n_k is number of alleles, or here twice the number of individuals in the populations
    r <- length( k2n )

    mc <- get_mem_lim_m_WC(m, n, r) # for now use the same formula as for WC (most major matrices are the same size)

    # initialize these things, unfortunately have to grow them for each chunk...
    FstTs <- c()
    FstBs <- c()
    
    # navigate chunks
    mci <- 1 # start of first chunk (needed for matrix inputs only; as opposed to function inputs)
    while(TRUE) { # start an infinite loop, break inside as needed
        if (isFn) {
            # get next "mc" SNPs
            Xi <- X( mc, indexes_keep )

            # stop when SNPs run out (only happens for functions X, not matrices)
            if (is.null(Xi))
                break
        } else {
            # here m is known...

            # this means all SNPs have been covered!
            if (mci > m)
                break
            
            is <- mci : min(mci + mc - 1, m) # range of SNPs to extract in this chunk

            # is this more efficient than setting indexes_keep<-1:n and subsetting?
            if (is.null(indexes_keep)) {
                if (loci_on_cols) {
                    Xi <- t(X[, is, drop = FALSE]) # transpose for our usual setup
                } else  {
                    Xi <- X[is, , drop = FALSE]
                }
            } else {
                if (loci_on_cols) {
                    Xi <- t(X[indexes_keep, is, drop = FALSE]) # transpose for our usual setup
                } else  {
                    Xi <- X[is, indexes_keep, drop = FALSE]
                }
            }
            mci <- mci + mc # update starting point for next chunk! (overshoots at the end, that's ok)
        }
        
        # this is a matrix that contains within-population allele frequencies
        k2ps <- do.call(
            cbind,
            lapply(
                k2is,
                function(is)
                    rowMeans(Xi[, is], na.rm = TRUE)
            )
        ) / 2
        
        # count number of non-NA individuals per SNP per pop k (then turn to allele counts with 2x-1)
        k2n <- 2 * do.call(
                       cbind,
                       lapply(
                           k2is,
                           function(is)
                               rowSums( !is.na(Xi[, is]), na.rm = TRUE )
                       )
                   ) - 1

        # estimate ancestral allele frequencies across the matrix
        # in WC, this is sample mean of each row of X, but in Hudson it is the sample mean of per-pop freqs
        p_anc_hat <- rowMeans( k2ps )

        # a vector of estimated variances
        # this terms appears in the top and bottom
        s2s <- rowSums( ( k2ps - p_anc_hat )^2 ) / (r-1)
        
        # compute Fst parts (tops and bottoms) for each locus using this method
        # bottom
        FstBsi <- p_anc_hat * ( 1 - p_anc_hat ) + s2s / r
        # top
        FstTsi <- s2s - rowMeans( k2ps * ( 1 - k2ps ) / k2n )

        # append to global data (for all chunks)
        FstTs <- c(FstTs, FstTsi)
        FstBs <- c(FstBs, FstBsi)
    }

    # compute global mean Fst! (the non-bootstrap final estimate across SNPs)
    fst <- sum(FstTs) / sum(FstBs)
    
    # also compute vector of per-locus Fst values
    fst_loci <- FstTs / FstBs
    
    # return an object with the mean and the raw data for bootstrapping (matrix of Fst parts: top and bottom, separating each SNP)
    return(
        list(
            fst = fst,
            fst_loci = fst_loci,
            data = cbind(FstTs, FstBs)
        )
    )
}

