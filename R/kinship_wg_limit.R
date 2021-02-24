#' Limit of Weir-Goudet kinship estimator
#'
#' This function calculates the biased limiting expectation of the Weir-Goudet kinship estimator.
#' This limit is easily calculated given the true kinship matrix.
#'
#' @param kinship The true kinship matrix
#'
#' @return The biased limit of the Weir-Goudet kinship estimator
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
#' kinship_wg_lim <- kinship_wg_limit(kinship)
#'
#' @export
kinship_wg_limit <- function( kinship ) {
    # make sure this is not missing
    if ( missing(kinship) )
        stop('`kinship` is required!')

    # validate input!
    popkin::validate_kinship( kinship )

    # WG bias is given by the mean of the non-diagonal elements!
    bias_wg <- mean( kinship[ lower.tri( kinship) ] )
    
    # return this rescaled kinship matrix
    ( kinship - bias_wg ) / ( 1 - bias_wg )
}
# bias is the same in terms of inbreeding!
# kinship_wg = ( kinship - bias_wg ) / ( 1 - bias_wg )
# (1 + f_wg)/2 = ( (1+f)/2 - bias_wg ) / ( 1 - bias_wg )
# 1 + f_wg = ( 1+f - 2 * bias_wg ) / ( 1 - bias_wg )
# f_wg = ( 1+f - 2 * bias_wg ) / ( 1 - bias_wg ) - 1
# f_wg = ( 1+f - 2 * bias_wg - 1 + bias_wg ) / ( 1 - bias_wg )
# f_wg = ( f - bias_wg ) / ( 1 - bias_wg )
