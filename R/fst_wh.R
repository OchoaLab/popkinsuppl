#' The Weir-Hill FST estimator
#'
#' This function implements the full FST estimator from Weir and Hill (2002) (equation 9).
#' Handles very large data, passed as BEDMatrix or as a regular R matrix.
#' Handles missing values correctly.
#'
#' @param X The genotype matrix (BEDMatrix, regular R matrix, or function, same as `popkin`).
#' @param labs A vector of subpopulation assignments for every individual.
#' At least two subpoplations must be present.
#' @param m The number of loci, required if `X` is a function (ignored otherwise).
#' In particular, `m` is obtained from `X` when it is a BEDMatrix or a regular R matrix.
#' @param ind_keep An optional vector of individuals to keep (as booleans or indexes, used to subset an R matrix).
#' @param loci_on_cols Determines the orientation of the genotype matrix (by default, `FALSE`, loci are along the rows).
#' If `X` is a BEDMatrix object, the input value is ignored (set automatically to `TRUE` internally).
#' @param mem_factor Proportion of available memory to use loading and processing genotypes.
#' Ignored if `mem_lim` is not `NA`.
#' @param mem_lim Memory limit in GB, used to break up genotype data into chunks for very large datasets.
#' Note memory usage is somewhat underestimated and is not controlled strictly.
#' Default in Linux and Windows is `mem_factor` times the free system memory, otherwise it is 1GB (OSX and other systems).
#' @param m_chunk_max Sets the maximum number of loci to process at the time.
#' Actual number of loci loaded may be lower if memory is limiting.
#'
#' @return A list with the following named elements, in this order:
#' 
#' - fst: The genome-wide Fst estimate (scalar).
#' - fst_loci: A vector of per-locus Fst estimates.
#' - data: a 2-by-m matrix of statistics used in estimating Fst.
#'   Useful to obtain a bootstrap distribution for the genome-wide Fst.
#' - maf: a vector of marginal allele frequency estimates, one for each locus.
#'   Note that it has not been converted to *minor* allele frequencies.
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
#' # estimate FST using the Weir-Hill formula
#' fst_wc_obj <- fst_wh(X, labs)
#' 
#' # the genome-wide FST estimate
#' fst_wc_obj$fst
#' 
#' # vector of per-locus FST estimates
#' fst_wc_obj$fst_loci
#' 
#' @seealso
#' The popkin package.
#'
#' @export
fst_wh <- function(
                   X,
                   labs,
                   m = NA,
                   ind_keep = NULL,
                   loci_on_cols = FALSE,
                   mem_factor = 0.7,
                   mem_lim = NA,
                   m_chunk_max = 1000
                   ) {
    if (missing(X))
        stop('Genotype matrix `X` is required!')
    if (missing(labs))
        stop('Subpopulation labels `labs` are required!')
    
    # this is always given through labs
    n <- length(labs)

    # shared with WC and Hudson estimator...
    obj <- preprocess_genotype_obj_fst(
        X = X,
        n = n,
        m = m,
        loci_on_cols = loci_on_cols,
        ind_keep = ind_keep
    )
    isFn <- obj$isFn
    m <- obj$m
    loci_on_cols <- obj$loci_on_cols
    
    # global computations shared by all chunks
    
    # general label processing to simplfy things
    out <- clean_labs(labs)
    k2is <- out$k2is
    k2n <- out$k2n
    r <- length(k2n)
    
    if (r < 2)
        stop('There must be two or more subpopulations!')

    # now compute these WH-specific things
    # weights for subpopulations
    k2w <- k2n / sum(k2n)
    # some other sample sizes that contain weights themselves
    k2nc <- k2n * (1 - k2w)
    # the sum of these things
    nc_sum <- sum( k2nc )
    # and even more normalized things, which appear on sums often
    # sorry for the cryptic names, there's no good descriptions of these monsters
    k2ndnc <- k2n / nc_sum
    k2wc <- k2nc / nc_sum
    k2wcu <- k2wc - k2w * (2 * k2n) / (2 * k2n - 1)
    
    # initialize these vectors
    # Do before get_mem_lim_m so free memory is accounted for properly
    FstTs <- vector('numeric', m)
    FstBs <- vector('numeric', m)
    mafs <- vector('numeric', m)
    
    # estimating total memory usage in bytes
    # p_anc_hat = m*8+ao
    # k2ps = m*r*8+mo
    # Xi[, k2is[[k]] ] = m*(n/r)*8+mo # on average, let's be conservative: m*n*8+mo
    # (k2ps-Pi)^2 = m*r*8+mo # intermediate calculation, unnamed?
    # s2s = m*8+ao # final length
    # p_q_hats = m*r*8+mo # new in WH
    # FstTsi = m*8+ao
    # FstBsi = m*8+ao
    
    # the general function doesn't have an "r" input, but since "n" only enters in m*n matrices and so does "r", we can hack it all to have an effective "n" that combines both
    # Have:
    # - 1.5 m*n (input int matrix, plus 1 subpop extract)
    # - 3 m*r
    # since there's no other n's this is equivalent to
    # n = 1.5 * n + 3 * r
    
    # calculate chunk size given available memory
    data <- popkin:::solve_m_mem_lim(
                         n = 1.5 * n + 3 * r,
                         m = m,
                         mat_m_n = 1, # left simple to not confound with above `n` hack
                         vec_m = 4,
                         mem = mem_lim,
                         mem_factor = mem_factor
                     )
    m_chunk <- data$m_chunk
    # cap value to a nice performing value (very good speed, minimal memory)
    if ( m_chunk > m_chunk_max )
        m_chunk <- m_chunk_max
    
    # navigate chunks
    i_chunk <- 1 # start of first chunk (needed for matrix inputs only; as opposed to function inputs)
    while(TRUE) { # start an infinite loop, break inside as needed
        # indexes to extract loci, and also so save to FstTs and FstBs vectors
        indexes_loci_chunk <- i_chunk : min(i_chunk + m_chunk - 1, m)
        
        if (isFn) {
            # get next "m_chunk" SNPs
            Xi <- X( m_chunk, ind_keep )

            # stop when SNPs run out (only happens for functions X, not matrices)
            if (is.null(Xi))
                break
            
        } else {
            # this means all SNPs have been covered!
            if (i_chunk > m)
                break
            
            # is this more efficient than setting ind_keep<-1:n and subsetting?
            if (is.null(ind_keep)) {
                if (loci_on_cols) {
                    Xi <- t(X[, indexes_loci_chunk, drop = FALSE]) # transpose for our usual setup
                } else {
                    Xi <- X[indexes_loci_chunk, , drop = FALSE]
                }
            } else {
                if (loci_on_cols) {
                    Xi <- t(X[ind_keep, indexes_loci_chunk, drop = FALSE]) # transpose for our usual setup
                } else  {
                    Xi <- X[indexes_loci_chunk, ind_keep, drop = FALSE]
                }
            }
        }
        
        # NOTE: in the other notes below I assume there aren't rows that are entirely NA
        # rowSums of all NAs is 0 (with na.rm = TRUE)
        # rowMeans of all NA is NaN (with na.rm = TRUE)

        # estimate ancestral allele frequencies across the matrix
        # NOTE: never NA (at least one genotype not NA per row)
        p_anc_hat <- rowMeans(Xi, na.rm = TRUE) / 2
        
        # this is a matrix that contains within-population allele frequencies
        # NOTE: rowMeans of all NA is NaN, so some cells of this k2ps matrix may be NaNs
        k2ps <- do.call(
            cbind,
            lapply(
                k2is,
                function(is)
                    rowMeans(Xi[, is, drop = FALSE], na.rm = TRUE)
            )
        ) / 2
        # a vector of estimated variances
        # s2s <- drop( ( k2ps - p_anc_hat )^2 %*% k2ndnc ) # OLD, not NA-robust
        # NOTE: this does not result in NAs
        s2s <- rowSums(
            sweep(
                ( k2ps - p_anc_hat )^2,
                MARGIN = 2,
                k2ndnc,
                `*`
            ),
            na.rm = TRUE
        )
        # and the p * (1-p) estimates per subpopulation (r by m matrix)
        # some of these may be NaNs
        p_q_hats <- k2ps * ( 1 - k2ps )
        
        # the actual estimates weighted as they appear in paper
        # NOTE: this does not result in NAs
        # top
        #FstTsi <- s2s + drop( p_q_hats %*% k2wcu ) # OLD, not NA-robust
        FstTsi <- s2s + rowSums(
                            sweep(
                                p_q_hats,
                                MARGIN = 2,
                                k2wcu,
                                `*`
                            ),
                            na.rm = TRUE
                        )
        # bottom
        #FstBsi <- s2s + drop( p_q_hats %*% k2wc ) # OLD, not NA-robust
        FstBsi <- s2s + rowSums(
                            sweep(
                                p_q_hats,
                                MARGIN = 2,
                                k2wc,
                                `*`
                            ),
                            na.rm = TRUE
                        )
        
        # copy chunk data to global vectors (containing all chunks)
        FstTs[ indexes_loci_chunk ] <- FstTsi
        FstBs[ indexes_loci_chunk ] <- FstBsi
        mafs[ indexes_loci_chunk ] <- p_anc_hat
        
        # update starting point for next chunk! (overshoots at the end, that's ok)
        i_chunk <- i_chunk + m_chunk
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
            data = cbind(FstTs, FstBs),
            maf = mafs
        )
    )
}

