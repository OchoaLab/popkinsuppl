#' Matrix of between-subpopulation Hudson FST estimates
#'
#' This function applies the `fst_hudson_k` formula to every pair of subpopulations in the dataset.
#' Since it is only applied to two subpopulations at the time, the values equal the original (non-generalized) "Hudson" FST estimator of Bhatia, Patterson, Sankararaman, and Price (2013).
#'
#' @param X The genotype matrix (BEDMatrix, regular R matrix, or function, same as `popkin`).
#' @param labs A vector of subpopulation assignments for every individual.
#' @param pops An optional vector of unique subpopulation labels, in the desired order for the matrix.
#' By default the unique labels in `labs` sorted alphabetically are used.
#' @param m The number of loci, required if `X` is a function (ignored otherwise).
#' In particular, `m` is obtained from `X` when it is a BEDMatrix or a regular R matrix.
#' @param loci_on_cols Determines the orientation of the genotype matrix (by default, `FALSE`, loci are along the rows).
#' If `X` is a BEDMatrix object, the input value is ignored (set automatically to `TRUE` internally).
#' @param mem_factor Proportion of available memory to use loading and processing genotypes.
#' Ignored if `mem_lim` is not `NA`.
#' @param mem_lim Memory limit in GB, used to break up genotype data into chunks for very large datasets.
#' Note memory usage is somewhat underestimated and is not controlled strictly.
#' Default in Linux and Windows is `mem_factor` times the free system memory, otherwise it is 1GB (OSX and other systems).
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
fst_hudson_pairwise <- function(
                                X,
                                labs,
                                pops = NULL,
                                m = NA,
                                loci_on_cols = FALSE,
                                mem_factor = 0.7,
                                mem_lim = NA
                                ) {
    if (missing(X))
        stop('Genotype matrix `X` is required!')
    if (missing(labs))
        stop('Subpopulation labels `labs` are required!')

    # default is to place subpops alphabetically in matrix
    if (is.null(pops))
        pops <- sort(unique(labs))
    
    k_subpops <- length(pops)

    # this is always given through labs
    n <- length(labs)

    # shared with WC estimator...
    obj <- preprocess_genotype_obj_fst(
        X = X,
        n = n,
        m = m,
        loci_on_cols = loci_on_cols
    )
    isFn <- obj$isFn
    m <- obj$m
    loci_on_cols <- obj$loci_on_cols

    # global computations shared by all chunks
    
    # general label processing to simplfy things
    out <- clean_labs(labs)
    k2is <- out$k2is
    k2n <- 2 * out$k2n - 1 # NOTE: in Hudson estimator, n_k is number of alleles, or here twice the number of individuals in the populations
    stopifnot(
        k_subpops == length( k2n )
    )

    if (k_subpops < 2)
        stop('There must be two or more subpopulations!')

    # construct a small matrix (between pops)
    # separate top and bottom matrices for before ratio
    # these contain running sums (will not save single locus statistics separately here)
    pwfst_hudson_tops <- matrix(0, nrow = k_subpops, ncol = k_subpops)
    pwfst_hudson_bots <- matrix(0, nrow = k_subpops, ncol = k_subpops)
    
    # calculate chunk size given available memory
    # use the same formula as for WC (the largest matrices are the same size)
    # TODO: this was almost entirely just copied from fst_hudson_k, but could optimize later...
    data <- popkin:::solve_m_mem_lim(
                         n = 1.5 * n + 2 * k_subpops,
                         m = m,
                         mat_m_n = 1, # left simple to not confound with above `n` hack
                         vec_m = 4, # 6,
                         mem = mem_lim,
                         mem_factor = mem_factor
                     )
    m_chunk <- data$m_chunk

        # navigate chunks
    i_chunk <- 1 # start of first chunk (needed for matrix inputs only; as opposed to function inputs)
    while(TRUE) { # start an infinite loop, break inside as needed
        # indexes to extract loci, and also so save to FstTs and FstBs vectors
        indexes_loci_chunk <- i_chunk : min(i_chunk + m_chunk - 1, m) # range of SNPs to extract in this chunk
        
        if (isFn) {
            # get next "m_chunk" SNPs
            Xi <- X( m_chunk )

            # stop when SNPs run out (only happens for functions X, not matrices)
            if (is.null(Xi))
                break
            
        } else {
            # this means all SNPs have been covered!
            if (i_chunk > m)
                break
            
            if (loci_on_cols) {
                Xi <- t(X[, indexes_loci_chunk, drop = FALSE]) # transpose for our usual setup
            } else {
                Xi <- X[indexes_loci_chunk, , drop = FALSE]
            }
        }

        # NOTE: in the other notes below I assume there aren't rows that are entirely NA
        # rowSums of all NAs is 0 (with na.rm = TRUE)
        # rowMeans of all NA is NaN (with na.rm = TRUE)

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
        
        # count number of non-NA individuals per SNP per pop k (then turn to allele counts with 2x-1)
        # NOTE: nothing here is NA, but there can be zeroes in the inner calls so the k2n matrix may have -1 cases
        k2n <- 2 * do.call(
                       cbind,
                       lapply(
                           k2is,
                           function(is)
                               rowSums( !is.na(Xi[, is, drop = FALSE]), na.rm = TRUE )
                       )
                   ) - 1

        ## # the bottoms are given by simple products of the allele frequencies
        ## # since k2ps may have NaNs, we have to create more careful versions without NAs
        ## k2psNN <- k2ps # copy
        ## k2psNN[ is.na(k2psNN) ] <- 0 # NAs become zeroes, so in p*q product they do not contribute
        ## k2qsNN <- 1 - k2ps # `q = 1 - p` version, copy first, also results in NAs (in the right places)
        ## k2qsNN[ is.na(k2qsNN) ] <- 0 # NAs become zeroes, so in p*q product they do not contribute (need NA q's to become zeroes too, otherwise p=0 gives q=1 != 0)
        ## # now this is non-NA entirely, and elements that would have been NA instead contribute zero to the running sum, as they should
        ## pq_mat <- crossprod( k2psNN, k2qsNN ) # k_subpops-squared matrix
        ## # add to running sum, add direct and transpose (as needed in formula)
        ## pwfst_hudson_bots <- pwfst_hudson_bots + pq_mat + t( pq_mat )

        # the tops have the squared allele frequency differences
        # plus small sample size corrections
        # to avoid confusion, I'll implement this first as a loop, then optimize if possible
        for (j in 2 : k_subpops) {
            # vector length m_chunk
            # some elements may be NaNs
            pj <- k2ps[, j]
            # small sample size correction terms
            # some elements may be NaNs too
            cj <- pj * ( 1 - pj ) / k2n[, j]
            # k < j, don't include self-comparison!
            for (k in 1 : (j-1)) {
                # vector length m_chunk
                pk <- k2ps[, k]
                ck <- pk * ( 1 - pk ) / k2n[, k]
                # the vector of values to add
                # elements to skip are precisely NA
                tjk <- ( pj - pk )^2 - cj - ck
                # turn into scalar, now never NA, store in one direction only (lower triangle)
                pwfst_hudson_tops[j, k] <- sum( tjk, na.rm = TRUE )

                # work on bottoms now
                # this is a vector of values to add
                # elements to skip are precisely NA
                bjk <- pj*(1-pk) + pk*(1-pj) 
                # turn into scalar, now never NA, store in one direction only (lower triangle)
                pwfst_hudson_bots[j, k] <- sum( bjk, na.rm = TRUE )
            }
        }

        # update starting point for next chunk! (overshoots at the end, that's ok)
        i_chunk <- i_chunk + m_chunk
    }

    # copy lower triangle to upper triangle
    # since the empty part is zeroes, this simple sum fixes this
    pwfst_hudson_tops <- pwfst_hudson_tops + t( pwfst_hudson_tops )
    pwfst_hudson_bots <- pwfst_hudson_bots + t( pwfst_hudson_bots )

    # perform final ratio
    pwfst_hudson <- pwfst_hudson_tops / pwfst_hudson_bots

    # the diagonal will be entirely 0/0 = NaN, fix that now
    diag( pwfst_hudson ) <- 0

    # add labels, can be helpful for manual inspection or some plots
    rownames( pwfst_hudson ) <- pops
    colnames( pwfst_hudson ) <- pops

    # return this!
    return( pwfst_hudson )
}

# kept slower, original version for validations only
fst_hudson_pairwise_hack <- function(X, labs, pops = NULL) {
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

