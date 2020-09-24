#' Standard kinship estimator
#'
#' This function constructs the standard kinship estimates for a given genotype matrix.
#' Handles very large data, passed as BEDMatrix or as a regular R matrix.
#' Handles missing values correctly.
#'
#' @param X The genotype matrix (BEDMatrix, regular R matrix, or function, same as `popkin`).
#' @param n The number of individuals.
#' Required if `X` is a function, ignored otherwise.
#' @param mean_of_ratios The standard kinship estimator can be computed in two broad forms.
#' If `FALSE` (default) the ratio-of-means version is computed, which behaves more favorably and has a known asymptotic bias.
#' If `TRUE`, the mean-of-ratios version is computed, which is more variable and has an uncharacterized bias, but is most common in the literature.
#' @param loci_on_cols Determines the orientation of the genotype matrix (by default, `FALSE`, loci are along the rows).
#' If `X` is a BEDMatrix object, the input value is ignored (set automatically to `TRUE` internally).
#' @param mem_factor Proportion of available memory to use loading and processing genotypes.
#' Ignored if `mem_lim` is not `NA`.
#' @param mem_lim Memory limit in GB, used to break up genotype data into chunks for very large datasets.
#' Note memory usage is somewhat underestimated and is not controlled strictly.
#' Default in Linux and Windows is `mem_factor` times the free system memory, otherwise it is 1GB (OSX and other systems).
#' @param want_M If `TRUE`, includes the matrix `M` of non-missing pair counts in the return value, which are sample sizes that can be useful in modeling the variance of estimates.
#' Default `FALSE` is to return the kinship matrix only.
#'
#' @return If `want_M` is `FALSE`, returns the estimated `n`-by-`n` kinship matrix only.
#' If `X` has names for the individuals, they will be copied to the rows and columns of this kinship matrix.
#' If `want_M` is `TRUE`, a named list is returned, containing:
#'
#' - `kinship`: the estimated `n`-by-`n` kinship matrix
#' - `M`: the `n`-by-`n` matrix of non-missing pair counts (see `want_M` option).
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
#' # estimate kinship matrices
#' # ... ratio-of-means version
#' kinship_rom <- kinship_std(X)
#' # ... mean-of-ratios version
#' kinship_mor <- kinship_std(X, mean_of_ratios = TRUE)
#' 
#' @seealso
#' The popkin package.
#'
#' @export
kinship_std <- function(
                        X,
                        n = NA,
                        mean_of_ratios = FALSE,
                        loci_on_cols = FALSE,
                        mem_factor = 0.7,
                        mem_lim = NA,
                        want_M = FALSE
                        ) {
    if ( missing( X ) )
        stop( 'Genotype matrix `X` is required!' )
    
    # SUPER low-mem version that processes input as it is read
    
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
    
    # estimating total memory usage in bytes
    # kinship = n*n*8+mo
    # M = n*n*8+mo
    # Xi = m*n*4+mo # initial int version
    # Pi = m*8+ao
    # Xi = m*n*8+mo # centered version is type double
    # is.na(Xi) = m*n*4+mo # ignore (NEWEST) or count once (NEW) or appears twice (OLD)
    # !is.na(Xi) = m*n*8+mo # gets turned into double when input to tcrossprod
    # M = n*n*8+mo # introduce again as we update it in addition
    # Xi = m*n*8+mo # ignore (NEWEST) or introduce again as it is updated (replacing NAs with zeroes; OLD)
    # kinship = n*n*8+mo # introduce again as we update it in addition
    # Pi*(1-Pi) = 2*(m*8+ao) # two temporary arrays in computing this
    
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
    
    # initialize desired matrix
    kinship <- matrix(0, nrow = n, ncol = n)
    M <- matrix(0, nrow = n, ncol = n) # normalization now varies per individual pair
    nu <- 0 # scalar normalization = mean p(1-p)
    m_nu <- 0 # to average nu (almost equal to m, but must remove completely fixed loci)
    
    # transfer names from X to kinship if present
    # this will carry over all the way to the final kinship matrix!
    # (M need not have names at all)
    if (!is.null(names_X)) {
        colnames(kinship) <- names_X
        rownames(kinship) <- names_X
    }
    
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

        # standard mean times half
        Pi <- rowMeans(Xi, na.rm = TRUE) / 2
        
        # variance estimate (length-i_chunk vector; factor of 2 or 4 added at the end)
        Vi <- Pi * ( 1 - Pi )
        
        # in the mean_of_ratios formulation it's extra critical to handle fixed loci correctly... but let's just do the same thing both ways
        indexes_not_fixed <- Vi > 0 # these are the good cases
        # filter everything if needed (will increase memory but it's easy to see things are handled correctly)
        if (any(!indexes_not_fixed)) {
            Xi <- Xi[indexes_not_fixed,]
            Pi <- Pi[indexes_not_fixed]
            Vi <- Vi[indexes_not_fixed]
        }

        # center before cross product...
        Xi <- Xi - 2 * Pi
        
        # some complicated steps are only necessary when there's missing data...
        if (anyNA(Xi)) {
            M <- M + crossprod( !is.na(Xi) ) # add non-missingness counts to pair count matrix
            Xi[is.na(Xi)] <- 0 # before applying cross product, to prevent NA errors, just set those values to zero and it works out!
        } else {
            M <- M + nrow(Xi) # this is correct denominator
        }
        
        if (mean_of_ratios) {
            # will average into A but scaling first!
            Xi <- Xi / sqrt(Vi) # this works! (scales each row as needed)
        } else {
            # work on scalar normalization
            nu <- nu + sum(Vi) # add these values, which are always defined
            # to normalize nu (every polymorphic SNP contributes a value in this case, regardless of missingness of particular individual pairs)
            m_nu <- m_nu <- length(Vi)
        }

        # cross product matrix at this SNP, add to running sum.  We'll add an extra -1 later... (this is computationally faster and maybe even more numerically stable)
        kinship <- kinship + crossprod(Xi)
        # NOTE: M and m count the same loci (including fixed loci in this case; as long as there's consistency the difference just cancels out as expected)
    }
    
    # normalize estimate!
    if (mean_of_ratios) {
        kinship <- kinship / M / 4 # just need per-pair average norm
    } else { 
        kinship <- kinship / M / ( 4 * nu / m_nu ) # turn into ratio of averages, return!
    }
    
    # figure out what to return
    if ( want_M ) {
        return(
            list(
                kinship = kinship,
                M = M
            )
        )
    } else {
        return( kinship )
    }
}

