#' Matrix of between-subpopulation Hudson FST estimates
#'
#' This function applies the `fst_hudson_k` formula to every pair of subpopulations in the dataset.
#' Since it is only applied to two subpopulations at the time, the values equal the original (non-generalized) "Hudson" FST estimator of Bhatia, Patterson, Sankararaman, and Price (2013).
#'
#' @param X The genotype matrix (BEDMatrix, regular R matrix, or function, same as `popkin`).
#' @param labs A vector of subpopulation assignments for every individual.
#' @param pops An optional vector of unique subpopulation labels, in the desired order for the matrix.
#' By default the unique labels in `labs` sorted alphabetically are used.
#'
#' @return A symmetric matrix of FST estimates between every pair of subpopulations.
#' The diagonal has zero values.
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
#' # estimated pairwise FST matrix using the "Hudson" formula
#' fst_hudson_matrix <- fst_hudson_pairwise(X, labs)
#' 
#' @seealso
#' The popkin package.
#'
#' @export
fst_hudson_pairwise <- function(X, labs, pops = NULL) {
    if (missing(X))
        stop('Genotype matrix `X` is required!')
    if (missing(labs))
        stop('Subpopulation labels `labs` are required!')

    # default is to place subpops alphabetically in matrix
    if (is.null(pops))
        pops <- sort(unique(labs))
    
    n <- length(pops)

    # construct a small matrix (between pops)
    pwfst_hudson <- matrix(0, nrow = n, ncol = n)
    
    for (i in 2:n) {
        pop_i <- pops[i]

        # j < i, don't include self-comparison!
        for (j in 1:(i-1)) {
            pop_j <- pops[j]
            # these are TRUE for individuals from these two populations only
            ind_keep <- labs %in% c(pop_i, pop_j)
            
            # and this is the pairwise Fst (because there are only two kinds of labels in input)
            # NOTE: X gets filtered by ind_keep, but labs should be pre-filtered!
            fst_ij <- fst_hudson_k( X, labs[ ind_keep ], ind_keep = ind_keep )$fst
            
            # store both ways
            pwfst_hudson[i, j] <- fst_ij
            pwfst_hudson[j, i] <- fst_ij
        }
    }
    
    # add labels, can be helpful for manual inspection or some plots
    rownames( pwfst_hudson ) <- pops
    colnames( pwfst_hudson ) <- pops

    # return this!
    return( pwfst_hudson )
}

