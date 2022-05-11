#' Limit of ROM version of the 2011 GCTA kinship (GRM) estimator
#'
#' This function calculates the biased limiting expectation of the GCTA kinship estimator described in Yang et al. (2011), which they refer to as a GRM instead, except the ratio-of-means (ROM) version (the above publication defines the mean-of-ratios (MOR) version instead; see Ochoa and Storey 2021 for the ROM version).
#' The GCTA software history shows that this exact estimator was abandoned in version 0.93.0 (8 Jul 2011) in favor of [kinship_std()], which remains in use as of writing (2022).
#' This limit is calculated from the true kinship matrix.
#' It agrees with [kinship_std_limit()] off-diagonal, but the diagonal bias is slightly different.
#'
#' @param kinship The true kinship matrix
#'
#' @return The biased limit of the GCTA kinship (GRM) estimator
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
#' kinship_gcta_lim <- kinship_gcta_limit(kinship)
#'
#' @seealso
#' Standard kinship estimator [kinship_std()] and the limit of the ROM version [kinship_std_limit()].
#' 
#' GCTA software, including history/update log.
#' https://yanglab.westlake.edu.cn/software/gcta/#Download
#'
#' @export
kinship_gcta_limit <- function( kinship ) {
    # make sure this is not missing
    if ( missing(kinship) )
        stop('`kinship` is required!')

    # validate input!
    popkin::validate_kinship(kinship)

    # extract number of individuals
    n_ind <- nrow(kinship)

    # formula requires a vector of 1's (even when weights are provided)
    u <- rep.int(1, n_ind)
    
    # weighted mean of each row
    mean_kinship_j <- rowMeans( kinship )
    
    # weighted average of entire matrix
    mean_kinship <- mean( mean_kinship_j )
    
    # return this "centered" matrix we expect given our bias calculation
    # inbr_diag fixed diagonal too! (but returns inbreeding instead of kinship)
    kinship_gcta_lim <- ( popkin::inbr_diag( kinship ) - tcrossprod(mean_kinship_j, u) - tcrossprod(u, mean_kinship_j) + mean_kinship ) / (1 - mean_kinship)

    # turn diagonal back to self-kinship
    diag( kinship_gcta_lim ) <- ( 1 + diag( kinship_gcta_lim ) ) / 2

    # done!
    return( kinship_gcta_lim )
}
