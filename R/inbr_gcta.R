#' 2011 GCTA inbreeding estimator III
#'
#' This function calculates the biased GCTA inbreeding estimator III described in Yang et al. (2011).
#' Though these estimates (MOR version) were the basis of the GRM diagonal according to that paper, the GCTA software history shows that this exact estimator was abandoned in version 0.93.0 (8 Jul 2011) in favor of [kinship_std()] (also MOR version), which remains in use as of writing (2022).
#'
#' @inheritParams kinship_std
#'
#' @return Inbreeding estimates
#'
#' @examples
#' # dimensions of simulated data
#' n_ind <- 100
#' m_loci <- 1000
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
#' # estimate inbreeding
#' # ... ROM version (see Ochoa and Storey (2021)).
#' inbr_gcta_rom <- inbr_gcta(X)
#' # ... MOR version (from Yang et al. (2011)).
#' inbr_gcta_mor <- inbr_gcta(X, mean_of_ratios = TRUE)
#' 
#' @seealso
#' GCTA 2011 GRM estimator ROM limit [kinship_gcta_limit()].
#' The limit of `inbr_gcta` with `mean_of_ratios = FALSE` is given by `popkin::inbr(kinship_gcta_limit(true_kinship))`.
#' 
#' Standard kinship estimator [kinship_std()] and the limit of the ROM version [kinship_std_limit()].
#' 
#' GCTA software, including history/update log.
#' https://yanglab.westlake.edu.cn/software/gcta/#Download
#'
#' @export
inbr_gcta <- function(
                      X,
                      n = NA,
                      mean_of_ratios = FALSE,
                      loci_on_cols = FALSE,
                      mem_factor = 0.7,
                      mem_lim = NA,
                      m_chunk_max = 1000 # gave good performance in tests
                      ) {
    if ( missing( X ) )
        stop( 'Genotype matrix `X` is required!' )

    # determine some behaviors depending on data type
    # first validate class and set key booleans
    isFn <- FALSE
    if ( is.function(X) ) {
        isFn <- TRUE
        if (is.na(n))
            stop('Number of individuals `n` is required when X is a function!')
    } else if ('BEDMatrix' %in% class(X)) {
        # same as general matrix but transposed
        # this is always imposed for this particular format!
        loci_on_cols <- TRUE
    } else if (!is.matrix(X)) {
        stop('X has unsupported class: ', toString( class(X) ) )
    } 

    # extract dimensions from data (not possible for function version)
    # also get individual names (IDs)
    names_X <- NULL # default
    if (isFn) {
        # have to define as NA to pass to get_mem_lim_m below
        m <- NA
    } else {
        if (loci_on_cols) {
            if (!is.na(n) && n != nrow(X)) 
                warning('User set number of samples that does not match X dimensions (will go with latter): ', n, ' != ', nrow(X))
            
            n <- nrow(X)
            m <- ncol(X)
            names_X <- rownames(X)
        } else {
            if (!is.na(n) && n != ncol(X)) 
                warning('User set number of samples that does not match X dimensions (will go with latter): ', n, ' != ', ncol(X))
            
            n <- ncol(X)
            m <- nrow(X)
            names_X <- colnames(X)
        }
    } 

    # initialize vectors of interest
    # calculate as running sum over chunks
    if ( mean_of_ratios ) {
        f_hat <- vector( 'numeric', n )
        m_c <- vector( 'numeric', n ) # complete observations
    } else {
        f_top <- vector( 'numeric', n )
        f_bot <- vector( 'numeric', n )
    }

    # estimating total memory usage in bytes
    data <- popkin:::solve_m_mem_lim(
                         n = n,
                         m = m,
                         mat_m_n = 2.5,
                         mat_n_n = 4,
                         vec_m = 3,
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
        if (isFn) {
            # get next "m_chunk" SNPs
            Xi <- X( m_chunk )

            # stop when SNPs run out (only happens for functions X, not matrices)
            if (is.null(Xi))
                break
        } else {
            # here m is known...

            # this means all SNPs have been covered!
            if (i_chunk > m)
                break

            # range of SNPs to extract in this chunk
            indexes_loci_chunk <- i_chunk : min(i_chunk + m_chunk - 1, m)
            
            if (loci_on_cols) {
                Xi <- t(X[, indexes_loci_chunk, drop = FALSE]) # transpose for our usual setup
            } else  {
                Xi <- X[indexes_loci_chunk, , drop = FALSE]
            }

            # update starting point for next chunk! (overshoots at the end, that's ok)
            i_chunk <- i_chunk + m_chunk
        }

        # calculate MAF
        p_est <- rowMeans( Xi, na.rm = TRUE ) / 2

        # identify fixed loci, exclude
        # (this is true for loci to keep)
        indexes_not_fixed <- 0 < p_est & p_est < 1
        if ( any( !indexes_not_fixed ) ) {
            p_est <- p_est[ indexes_not_fixed ]
            Xi <- Xi[ indexes_not_fixed, ]
        }

        # this is an efficient way to compute terms of interest, formulated as self-kinship
        # note this is entrywise product!
        # - to save memory, overwrite Xi (original Xi no longer needed)
        # - NA values in original Xi are NAs in this new Xi
        Xi <- ( Xi - 1 ) * ( Xi - 2 * p_est )

        # this is vector of estimates
        # note that estimates are actually 2x self-kinship at this stage
        if ( mean_of_ratios ) {
            # normalize each locus, then sum excluding missing loci
            f_hat <- f_hat + colSums( Xi / ( p_est * ( 1 - p_est ) ), na.rm = TRUE )
            # add to counts of non-missing observations
            m_c <- m_c + colSums( !is.na( Xi ) )
        } else {
            # sum excluding missing loci
            f_top <- f_top + colSums( Xi, na.rm = TRUE )
            # sum these too, but include for each individual only if its observation was not missing
            f_bot <- f_bot + colSums( ( !is.na( Xi ) ) * ( p_est * ( 1 - p_est ) ) )
        }
    }

    # when done, complete normalizations
    if ( mean_of_ratios ) {
        # get average from sum by dividing by number of complete observations
        f_hat <- f_hat / m_c / 2
    } else {
        # get final ratio by dividing both parts
        f_hat <- f_top / f_bot / 2
    }

    # transfer names from X to inbr if present
    if (!is.null(names_X))
        names( f_hat ) <- names_X
    
    # and very lastly, turn 2x-self-kinship to inbreeding by subtracting 1
    return ( f_hat - 1 )
}
