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
#' If `TRUE` (default) the ratio-of-means version is computed, which behaves more favorably and has a known asymptotic bias.
#' If `FALSE`, the mean-of-ratios version is computed, which is more variable and has an uncharacterized bias, but is most common in the literature.
#' @param loci_on_cols Determines the orientation of the genotype matrix (by default, `FALSE`, loci are along the rows).
#' Set automatically to `TRUE` if `X` is a BEDMatrix object.
#'
#' @return The estimated kinship matrix.
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
kinship_std <- function(X, n = NA, mean_of_ratios = FALSE, loci_on_cols = FALSE) {
    # SUPER low-mem version that processes input as it is read
    # NOTE: below I kept "A" notation, but this is kinship_std
    
    # determine some behaviors depending on data type
    # first validate class and set key booleans
    isFn <- FALSE
    if (class(X) == 'function') {
        isFn <- TRUE
        if (is.na(n))
            stop('Number of individuals `n` is required when X is a function!')
    } else if (class(X) == 'BEDMatrix') {
        # same as general matrix but transposed
        # this is always imposed for this particular format!
        loci_on_cols <- TRUE
    } else if (!is.matrix(X)) {
        stop('X has unsupported class: ', class(X))
    } 

    # extract dimensions from data (not possible for function version)
    if (isFn) {
        # have to define as NA to pass to get_mem_lim_m below
        m <- NA
    } else {
        if (loci_on_cols) {
            if (!is.na(n) && n != nrow(X)) 
                warning('User set number of samples that does not match X dimensions (will go with latter): ', n, ' != ', nrow(X))
            
            n <- nrow(X)
            m <- ncol(X)
        } else {
            if (!is.na(n) && n != ncol(X)) 
                warning('User set number of samples that does not match X dimensions (will go with latter): ', n, ' != ', ncol(X))
            
            n <- ncol(X)
            m <- nrow(X)
        }
    } 
    
    mc <- get_mem_lim_m_kinship_std(m, n)
    
    # initialize desired matrix
    A <- matrix(0, nrow = n, ncol = n)
    M <- matrix(0, nrow = n, ncol = n) # normalization now varies per individual pair
    nu <- 0 # scalar normalization = mean p(1-p)
    m_nu <- 0 # to average nu (almost equal to m, but must remove completely fixed loci)
    
    # navigate chunks
    mci <- 1 # start of first chunk (needed for matrix inputs only; as opposed to function inputs)
    while(TRUE) { # start an infinite loop, break inside as needed
        if (isFn) {
            # get next "mc" SNPs
            Xi <- X( mc )

            # stop when SNPs run out (only happens for functions X, not matrices)
            if (is.null(Xi))
                break
        } else {
            # here m is known...

            # this means all SNPs have been covered!
            if (mci > m)
                break

            # range of SNPs to extract in this chunk
            is <- mci : min(mci + mc - 1, m)
            
            if (loci_on_cols) {
                Xi <- t(X[, is, drop = FALSE]) # transpose for our usual setup
            } else  {
                Xi <- X[is, , drop = FALSE]
            }

            # update starting point for next chunk! (overshoots at the end, that's ok)
            mci <- mci + mc
        }

        # standard mean times half
        Pi <- rowMeans(Xi, na.rm = TRUE) / 2
        
        # variance estimate (length-mci vector; factor of 2 or 4 added at the end)
        Vi <- Pi * ( 1 - Pi )
        
        # in the mean_of_ratios formulation it's extra critical to handle fixed loci correctly... but let's just do the same thing both ways
        is <- Vi > 0 # these are the good cases
        # filter everything if needed (will increase memory but it's easy to see things are handled correctly)
        if (any(!is)) {
            Xi <- Xi[is,]
            Pi <- Pi[is]
            Vi <- Vi[is]
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
        A <- A + crossprod(Xi)
        # NOTE: M and m count the same loci (including fixed loci in this case; as long as there's consistency the difference just cancels out as expected)
    }
    
    # normalize estimate!
    if (mean_of_ratios) {
        A <- A / M / 4 # just need per-pair average norm
    } else { 
        A <- A / M / ( 4 * nu / m_nu ) # turn into ratio of averages, return!
    }
    
    return(A)
}

