# a simple version to compare our more complex code against
# this will only work with regular matrix
# - missingness is handled
# - fixed loci will be excluded
inbr_gcta_basic <- function(
                            X,
                            mean_of_ratios = FALSE
                            ) {
    if ( missing( X ) )
        stop( 'Genotype matrix `X` is required!' )

    n <- ncol(X)
    m <- nrow(X)
    names_X <- colnames(X)

    # calculate MAF
    p_est <- rowMeans( X, na.rm = TRUE ) / 2

    # identify fixed loci, exclude
    # (this is true for loci to keep)
    indexes_not_fixed <- 0 < p_est & p_est < 1
    if ( any( !indexes_not_fixed ) ) {
        p_est <- p_est[ indexes_not_fixed ]
        X <- X[ indexes_not_fixed, ]
    }

    # direct formula from paper, in inbreeding form!
    # (probably less efficient and numerically stable than our implementation)
    # - to save memory, overwrite X (original X no longer needed)
    # - NA values in original X are NAs in this new X
    X <- X^2 - ( 1 + 2 * p_est ) * X + 2 * p_est^2
    
    # this is vector of estimates
    # note that estimates are actually self-kinship at this stage
    if ( mean_of_ratios ) {
        # normalize each locus, then sum excluding missing loci
        f_hat <- colMeans( X / ( p_est * ( 1 - p_est ) ), na.rm = TRUE ) / 2
    } else {
        # sum excluding missing loci
        # sum these too, but include for each individual only if its observation was not missing
        f_hat <- colSums( X, na.rm = TRUE ) / colSums( ( !is.na( X ) ) * ( p_est * ( 1 - p_est ) ) ) / 2
    }
    
    # transfer names from X to inbr if present
    if ( !is.null( names_X ) )
        names( f_hat ) <- names_X
    
    return( f_hat )
}
