#' Limit of standard kinship estimator
#'
#' This function calculates the biased limiting expectation of the standard kinship estimator.
#' This limit is easily calculated given the true kinship matrix and the weights used to estimate the ancestral allele frequency.
#'
#' @param kinship The true kinship matrix
#' @param weights Optional weights for individuals
#'
#' @return The biased limit of the standard kinship estimator
#'
#' @examples
#' # create a dummy kinship matrix
#' kinship <- matrix(
#'     c(
#'         0.6, 0.1, 0,
#'         0.1, 0.6, 0.1,
#'         0, 0.1, 0.6
#'         ),
#'     nrow = 3
#' )
#' # this is its biased limit
#' # (uniform weights)
#' kinship_biased_limit <- kinship_std_limit(kinship)
#'
#' @seealso
#' Standard kinship estimator [kinship_std()].
#'
#' @export
kinship_std_limit <- function(kinship, weights = NULL) {
    # make sure this is not missing
    if ( missing(kinship) )
        stop('`kinship` is required!')

    # validate input!
    popkin::validate_kinship(kinship)

    # extract number of individuals
    n_ind <- nrow(kinship)

    # formula requires a vector of 1's (even when weights are provided)
    u <- rep.int(1, n_ind)
    
    # default is uniform weights
    if ( is.null(weights) ) {
        weights <- u / n_ind
    } else {
        # validate dimensions if weights were defined
        if (length(weights) != n_ind)
            stop('Length of `weights` does not agree with kinship matrix: ', length(weights), ' != ', n_ind)
        # renormalize for good measure
        weights <- weights / sum(weights)
    }
    
    # weighted mean of each row
    mean_kinship_j <- drop( kinship %*% weights )
    
    # weighted average of entire matrix
    mean_kinship <- drop( weights %*% mean_kinship_j )
    
    # return this "centered" matrix we expect given our bias calculation
    ( kinship - tcrossprod(mean_kinship_j, u) - tcrossprod(u, mean_kinship_j) + mean_kinship ) / (1 - mean_kinship)
}
