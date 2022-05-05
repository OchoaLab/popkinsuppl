# construct some random data to make sure all functions run on it and return objects in the expected dimensions

# dimensions of simulated data
n_ind <- 100
m_loci <- 1000
k_subpops <- 10
n_data <- n_ind * m_loci

# missingness rate
miss <- 0.1

# simulate ancestral allele frequencies
# uniform (0,1)
# it'll be ok if some of these are zero
p_anc <- runif(m_loci)

# simulate some binomial data
X <- rbinom(n_data, 2, p_anc)

# sprinkle random missingness
X[ sample(X, n_data * miss) ] <- NA

# turn into a matrix
X <- matrix(X, nrow = m_loci, ncol = n_ind)

# remove rows that are entirely missing and/or fixed
p_anc_hat <- rowMeans( X, na.rm = TRUE ) # MAFs
indexes <- !is.na( p_anc_hat ) & p_anc_hat > 0 & p_anc_hat < 1 # keep these
# apply filters
X <- X[ indexes, ]
# update this to match
m_loci <- nrow( X )

# create subpopulation labels
# k_subpops groups of equal size
labs <- ceiling( (1 : n_ind) / n_ind * k_subpops )

# create a more difficult case where there are lots of small subpopulations
# here each subpopulation nas one individual only, creating lots of subpops with NAs
labs1 <- 1 : n_ind


test_that("kinship_std ROM works", {
    # test ratio-of-means version
    kinship <- kinship_std(X)

    expect_true( is.numeric(kinship) )
    expect_true( !anyNA(kinship) )

    # test dimensions
    expect_true( is.matrix(kinship) )
    expect_equal( nrow(kinship), n_ind)
    expect_equal( ncol(kinship), n_ind)

    # repeat with want_M == TRUE
    obj <- kinship_std( X, want_M = TRUE )
    kinship <- obj$kinship
    M <- obj$M

    expect_true( is.numeric(kinship) )
    expect_true( !anyNA(kinship) )
    expect_true( is.numeric(M) )
    expect_true( !anyNA(M) )
    
    # test dimensions
    expect_true( is.matrix(kinship) )
    expect_equal( nrow(kinship), n_ind )
    expect_equal( ncol(kinship), n_ind )
    expect_true( is.matrix(M) )
    expect_equal( nrow(M), n_ind )
    expect_equal( ncol(M), n_ind )
    
    # M (pairwise sample sizes, so excluding NA pairs) must satisfy obvious range limits
    expect_true( min(M) >= 0 )
    expect_true( max(M) <= m_loci )

    # repeat with non-default m_chunk_max
    expect_silent(
        kinship2 <- kinship_std( X, m_chunk_max = 1 )
    )
    # outputs should have been the same
    expect_equal( kinship, kinship2 )
})

test_that("kinship_std MOR works", {
    # test mean-of-ratios version
    kinship <- kinship_std(X, mean_of_ratios = TRUE)
    
    expect_true( is.numeric(kinship) )
    expect_true( !anyNA(kinship) )
    
    # test dimensions
    expect_true( is.matrix(kinship) )
    expect_equal( nrow(kinship), n_ind)
    expect_equal( ncol(kinship), n_ind)

    # repeat with want_M == TRUE
    obj <- kinship_std( X, mean_of_ratios = TRUE, want_M = TRUE )
    kinship <- obj$kinship
    M <- obj$M

    expect_true( is.numeric(kinship) )
    expect_true( !anyNA(kinship) )
    expect_true( is.numeric(M) )
    expect_true( !anyNA(M) )
    
    # test dimensions
    expect_true( is.matrix(kinship) )
    expect_equal( nrow(kinship), n_ind )
    expect_equal( ncol(kinship), n_ind )
    expect_true( is.matrix(M) )
    expect_equal( nrow(M), n_ind )
    expect_equal( ncol(M), n_ind )
    
    # M (pairwise sample sizes, so excluding NA pairs) must satisfy obvious range limits
    expect_true( min(M) >= 0 )
    expect_true( max(M) <= m_loci )
    
    # repeat with non-default m_chunk_max
    expect_silent(
        kinship2 <- kinship_std( X, mean_of_ratios = TRUE, m_chunk_max = 1 )
    )
    # outputs should have been the same
    expect_equal( kinship, kinship2 )
})

test_that("fst_wc works", {
    # estimate FST using the Weir-Cockerham formula
    obj <- fst_wc(X, labs)

    # test overall object
    expect_true( is.list(obj) )
    expect_equal( length(obj), 4 )
    expect_equal( names(obj), c('fst', 'fst_loci', 'data', 'maf') )
    
    # the genome-wide FST estimate
    expect_true( is.numeric(obj$fst) )
    expect_equal( length(obj$fst), 1 )
    expect_true( !is.na(obj$fst) )
    
    # vector of per-locus FST estimates
    expect_true( is.numeric(obj$fst_loci) )
    expect_equal( length(obj$fst_loci), m_loci )
    expect_true( !anyNA( obj$fst_loci ) )
    
    # vector of MAFs
    expect_true( is.numeric(obj$maf) )
    expect_equal( length(obj$maf), m_loci )
    expect_true( !anyNA( obj$maf ) )

    # FST data matrix
    expect_true( is.numeric(obj$data) )
    expect_true( is.matrix(obj$data) )
    expect_equal( nrow(obj$data), m_loci )
    expect_equal( ncol(obj$data), 2 )
    expect_true( !anyNA( obj$data ) )
})

test_that("fst_wc with `FIT = TRUE` works", {
    # estimate FST using the Weir-Cockerham formula
    obj <- fst_wc(X, labs, FIT = TRUE)

    # test overall object
    expect_true( is.list(obj) )
    expect_equal( length(obj), 6 )
    expect_equal( names(obj), c('fst', 'fst_loci', 'data', 'maf', 'fit', 'fit_loci') )
    
    # the genome-wide FST estimate
    expect_true( is.numeric(obj$fst) )
    expect_equal( length(obj$fst), 1 )
    expect_true( !is.na(obj$fst) )
    # ditto FIT
    expect_true( is.numeric(obj$fit) )
    expect_equal( length(obj$fit), 1 )
    expect_true( !is.na(obj$fit) )
    
    # vector of per-locus FST estimates
    expect_true( is.numeric(obj$fst_loci) )
    expect_equal( length(obj$fst_loci), m_loci )
    expect_true( !anyNA( obj$fst_loci ) )
    # ditto FIT
    expect_true( is.numeric(obj$fit_loci) )
    expect_equal( length(obj$fit_loci), m_loci )
    expect_true( !anyNA( obj$fit_loci ) )
    
    # vector of MAFs
    expect_true( is.numeric(obj$maf) )
    expect_equal( length(obj$maf), m_loci )
    expect_true( !anyNA( obj$maf ) )

    # FST data matrix
    expect_true( is.numeric(obj$data) )
    expect_true( is.matrix(obj$data) )
    expect_equal( nrow(obj$data), m_loci )
    expect_equal( ncol(obj$data), 3 ) # this version has 3 columns!
    expect_true( !anyNA( obj$data ) )
})

## test_that("fst_wc works singleton subpops", {
##     # estimate FST using the Weir-Cockerham formula
##     obj <- fst_wc(X, labs1)
##
##     # test overall object
##     expect_true( is.list(obj) )
##     expect_equal( length(obj), 4 )
##     expect_equal( names(obj), c('fst', 'fst_loci', 'data', 'maf') )
##    
##     # the genome-wide FST estimate
##     expect_true( is.numeric(obj$fst) )
##     expect_equal( length(obj$fst), 1 )
##     expect_true( !is.na(obj$fst) )
##    
##     # vector of per-locus FST estimates
##     expect_true( is.numeric(obj$fst_loci) )
##     expect_equal( length(obj$fst_loci), m_loci )
##     expect_true( !anyNA( obj$fst_loci ) )
##    
##     # vector of MAFs
##     expect_true( is.numeric(obj$maf) )
##     expect_equal( length(obj$maf), m_loci )
##     expect_true( !anyNA( obj$maf ) )
##
##     # FST data matrix
##     expect_true( is.numeric(obj$data) )
##     expect_true( is.matrix(obj$data) )
##     expect_equal( nrow(obj$data), m_loci )
##     expect_equal( ncol(obj$data), 2 )
##     expect_true( !anyNA( obj$data ) )
## })

test_that("fst_wh works", {
    # estimate FST using the Weir-Hill formula
    obj <- fst_wh(X, labs)

    # test overall object
    expect_true( is.list(obj) )
    expect_equal( length(obj), 4 )
    expect_equal( names(obj), c('fst', 'fst_loci', 'data', 'maf') )
    
    # the genome-wide FST estimate
    expect_true( is.numeric(obj$fst) )
    expect_equal( length(obj$fst), 1 )
    expect_true( !is.na(obj$fst) )
    
    # vector of per-locus FST estimates
    expect_true( is.numeric(obj$fst_loci) )
    expect_equal( length(obj$fst_loci), m_loci )
    expect_true( !anyNA( obj$fst_loci ) )
    
    # vector of MAFs
    expect_true( is.numeric(obj$maf) )
    expect_equal( length(obj$maf), m_loci )
    expect_true( !anyNA( obj$maf ) )

    # FST data matrix
    expect_true( is.numeric(obj$data) )
    expect_true( is.matrix(obj$data) )
    expect_equal( nrow(obj$data), m_loci )
    expect_equal( ncol(obj$data), 2 )
    expect_true( !anyNA( obj$data ) )
})

test_that("fst_wh works singleton subpops", {
    # estimate FST using the Weir-Hill formula
    obj <- fst_wh(X, labs1)

    # test overall object
    expect_true( is.list(obj) )
    expect_equal( length(obj), 4 )
    expect_equal( names(obj), c('fst', 'fst_loci', 'data', 'maf') )
    
    # the genome-wide FST estimate
    expect_true( is.numeric(obj$fst) )
    expect_equal( length(obj$fst), 1 )
    expect_true( !is.na(obj$fst) )
    
    # vector of per-locus FST estimates
    expect_true( is.numeric(obj$fst_loci) )
    expect_equal( length(obj$fst_loci), m_loci )
    expect_true( !anyNA( obj$fst_loci ) )
    
    # vector of MAFs
    expect_true( is.numeric(obj$maf) )
    expect_equal( length(obj$maf), m_loci )
    expect_true( !anyNA( obj$maf ) )

    # FST data matrix
    expect_true( is.numeric(obj$data) )
    expect_true( is.matrix(obj$data) )
    expect_equal( nrow(obj$data), m_loci )
    expect_equal( ncol(obj$data), 2 )
    expect_true( !anyNA( obj$data ) )
})

test_that("fst_hudson_k works", {
    # estimate FST using the Weir-Cockerham formula
    obj <- fst_hudson_k(X, labs)

    # test overall object
    expect_true( is.list(obj) )
    expect_equal( length(obj), 3 )
    expect_equal( names(obj), c('fst', 'fst_loci', 'data') )
    
    # the genome-wide FST estimate
    expect_true( is.numeric(obj$fst) )
    expect_equal( length(obj$fst), 1 )
    expect_true( !is.na(obj$fst) )
    
    # vector of per-locus FST estimates
    expect_true( is.numeric(obj$fst_loci) )
    expect_equal( length(obj$fst_loci), m_loci )
    expect_true( !anyNA( obj$fst_loci ) )
    
    # FST data matrix
    expect_true( is.numeric(obj$data) )
    expect_true( is.matrix(obj$data) )
    expect_equal( nrow(obj$data), m_loci )
    expect_equal( ncol(obj$data), 2 )
    expect_true( !anyNA( obj$data ) )
})

test_that("fst_hudson_k works singleton subpops", {
    # estimate FST using the Weir-Cockerham formula
    obj <- fst_hudson_k(X, labs1)

    # test overall object
    expect_true( is.list(obj) )
    expect_equal( length(obj), 3 )
    expect_equal( names(obj), c('fst', 'fst_loci', 'data') )
    
    # the genome-wide FST estimate
    expect_true( is.numeric(obj$fst) )
    expect_equal( length(obj$fst), 1 )
    expect_true( !is.na(obj$fst) )
    
    # vector of per-locus FST estimates
    expect_true( is.numeric(obj$fst_loci) )
    expect_equal( length(obj$fst_loci), m_loci )
    expect_true( !anyNA( obj$fst_loci ) )
    
    # FST data matrix
    expect_true( is.numeric(obj$data) )
    expect_true( is.matrix(obj$data) )
    expect_equal( nrow(obj$data), m_loci )
    expect_equal( ncol(obj$data), 2 )
    expect_true( !anyNA( obj$data ) )
})

test_that("fst_hudson_pairwise works", {
    # estimated pairwise FST matrix using the "Hudson" formula
    fst_hudson_matrix <- fst_hudson_pairwise(X, labs)
    # compare to an older version that's also been validated, and which relies on fst_hudson_k internally so it's extra validated since fst_hudson_k was tested above
    fst_hudson_matrix_hack <- fst_hudson_pairwise_hack(X, labs)
    
    expect_true( is.numeric(fst_hudson_matrix) )

    # there should be no NAs
    expect_true( all( !is.na(fst_hudson_matrix) ) )

    # test dimensions
    expect_true( is.matrix(fst_hudson_matrix) )
    expect_equal( nrow(fst_hudson_matrix), k_subpops)
    expect_equal( ncol(fst_hudson_matrix), k_subpops)

    # diagonal must be zero
    expect_equal( range(diag(fst_hudson_matrix)), c(0, 0) )

    # compare to slower version
    expect_equal( fst_hudson_matrix, fst_hudson_matrix_hack )
})

test_that("fst_hudson_pairwise works singleton subpops", {
    # let's construct a smaller example (fewer individual pairs), otherwise this test is way too slow
    
    # dimensions of simulated data
    n_ind <- 10
    m_loci <- 100
    n_data <- n_ind * m_loci

    # missingness rate
    miss <- 0.1

    # simulate ancestral allele frequencies
    # uniform (0,1)
    # it'll be ok if some of these are zero
    p_anc <- runif(m_loci)

    # simulate some binomial data
    X <- rbinom(n_data, 2, p_anc)

    # sprinkle random missingness
    X[ sample(X, n_data * miss) ] <- NA

    # turn into a matrix
    X <- matrix(X, nrow = m_loci, ncol = n_ind)

    # create a more difficult case where there are lots of small subpopulations
    # here each subpopulation nas one individual only, creating lots of subpops with NAs
    labs1 <- 1 : n_ind
    
    # estimated pairwise FST matrix using the "Hudson" formula
    fst_hudson_matrix <- fst_hudson_pairwise(X, labs1)
    
    expect_true( is.numeric(fst_hudson_matrix) )

    # there should be no NAs
    expect_true( all( !is.na(fst_hudson_matrix) ) )

    # test dimensions
    expect_true( is.matrix(fst_hudson_matrix) )
    expect_equal( nrow(fst_hudson_matrix), n_ind )
    expect_equal( ncol(fst_hudson_matrix), n_ind )

    # diagonal must be zero
    expect_equal( range(diag(fst_hudson_matrix)), c(0, 0) )
})


###############################################

# tests restricted to data:
# - without missingness and
# - without fixed loci

# dimensions of simulated data
n_ind <- 100
m_loci <- 1000

# simulate ancestral allele frequencies
# uniform (0,1)
# it'll be ok if some of these are zero
p_anc <- runif(m_loci)

# simulate some binomial data
X <- rbinom(n_ind * m_loci, 2, p_anc)

# turn into a matrix
X <- matrix(X, nrow = m_loci, ncol = n_ind)

# find and remove fixed loci
p_anc_hat <- rowMeans(X) / 2
indexes <- 0 < p_anc_hat & p_anc_hat < 1
X <- X[indexes, ]

test_that("kinship_std ROM agrees with naive formula (requires no fixed loci)", {
    # test ratio-of-means version
    kinship <- kinship_std(X)

    # the "naive" formula
    p_anc_hat <- rowMeans(X) / 2
    # center, copying into new matrix
    X2 <- X - 2 * p_anc_hat
    # estimate p(1-p)
    p_q_hat <- 4 * p_anc_hat * (1 - p_anc_hat)
    # return "cross product" normalized as ratio of means
    kinship_naive <- crossprod(X2) / sum( p_q_hat )

    # finally, compare estimates!
    expect_equal( kinship, kinship_naive )
})

test_that("kinship_std MOR agrees with naive formula (requires no fixed loci)", {
    # test mean-of-ratios version
    kinship <- kinship_std(X, mean_of_ratios = TRUE)

    # the "naive" formula
    p_anc_hat <- rowMeans(X) / 2
    # estimate standard deviation
    stddev_hat <- 2 * sqrt( p_anc_hat * (1 - p_anc_hat) )
    # center, copying into new matrix
    X2 <- X - 2 * p_anc_hat
    # scale
    X2 <- X2 / stddev_hat
    # return "cross product", normalized
    kinship_naive <- crossprod(X2) / nrow(X2)
    
    # finally, compare estimates!
    expect_equal( kinship, kinship_naive )
})

###############################################

# tests restricted to data:
# - without missingness and
# - fixed loci allowed

# dimensions of simulated data
n_ind <- 100
m_loci <- 1000

# simulate ancestral allele frequencies
# uniform (0,1)
# it'll be ok if some of these are zero
p_anc <- runif(m_loci)

# simulate some binomial data
X <- rbinom(n_ind * m_loci, 2, p_anc)

# turn into a matrix
X <- matrix(X, nrow = m_loci, ncol = n_ind)

test_that("kinship_std ROM agrees with naive formula (admits fixed loci)", {
    # test ratio-of-means version
    kinship <- kinship_std(X)

    # the "naive" formula
    p_anc_hat <- rowMeans(X) / 2
    # estimate p(1-p)
    p_q_hat <- 4 * p_anc_hat * (1 - p_anc_hat)
    # copy data to edit
    X2 <- X
    # remove fixed loci (will make estimate fail)
    # these are the ones to keep
    i_keep <- 0 < p_q_hat
    if (any(!i_keep)) {
        # filter all data
        X2 <- X2[i_keep, ]
        p_anc_hat <- p_anc_hat[i_keep]
        p_q_hat <- p_q_hat[i_keep]
    }
    # center
    X2 <- X2 - 2 * p_anc_hat
    # return "cross product" normalized as ratio of means
    kinship_naive <- crossprod(X2) / sum( p_q_hat )

    # finally, compare estimates!
    expect_equal( kinship, kinship_naive )
})

test_that("kinship_std MOR agrees with naive formula (admits fixed loci)", {
    # test mean-of-ratios version
    kinship <- kinship_std(X, mean_of_ratios = TRUE)

    # the "naive" formula
    p_anc_hat <- rowMeans(X) / 2
    # estimate standard deviation
    stddev_hat <- 2 * sqrt( p_anc_hat * (1 - p_anc_hat) )
    # remove fixed loci (will make estimate fail)
    X2 <- X # copy in case edits are needed
    # these are the ones to keep
    i_keep <- 0 < stddev_hat
    if (any(stddev_hat <= 0)) {
        # filter all data
        X2 <- X2[i_keep, ]
        p_anc_hat <- p_anc_hat[i_keep]
        stddev_hat <- stddev_hat[i_keep]
    }
    # center
    X2 <- X2 - 2 * p_anc_hat
    # scale
    X2 <- X2 / stddev_hat
    # return "cross product", normalized
    kinship_naive <- crossprod(X2) / nrow(X2)
    
    # finally, compare estimates!
    expect_equal( kinship, kinship_naive )
})

test_that("kinship_std_limit works", {
    # construct a dummy kinship matrix
    kinship <- matrix(
        c(
            0.6, 0.1, 0,
            0.1, 0.6, 0.1,
            0, 0.1, 0.6
        ),
        nrow = 3
    )
    # this is its biased limit
    # (uniform weights)
    expect_silent(
        kinship_biased_limit <- kinship_std_limit(kinship)
    )
    
    # validate output!
    expect_silent( popkin::validate_kinship(kinship_biased_limit) )

    # match dims (already tested to be square by validate_kinship, so just check once)
    expect_equal( nrow(kinship), nrow( kinship_biased_limit ) )

    # zero overall mean
    expect_equal( mean( kinship_biased_limit ), 0 )

    # zero mean in every row
    # this max(abs(x)) boils it down to a single comparison
    expect_equal( max(abs(rowMeans( kinship_biased_limit ))), 0 )
})

test_that( "kinship_std agrees with kinship_std_limit( popkin(.) ) ROM and MOR", {
    # recent algebra suggests these should be the same!
    # test totally passes for data without missingness, but not with missingness
    expect_equal( kinship_std( X ), kinship_std_limit( popkin::popkin( X ) ) )
    
    # repeat for MOR version
    expect_equal( kinship_std( X, mean_of_ratios = TRUE ), kinship_std_limit( popkin::popkin( X, mean_of_ratios = TRUE ) ) )
})

test_that( "kinship_gcta_limit works", {
    # construct a dummy kinship matrix
    kinship <- matrix(
        c(
            0.6, 0.1, 0,
            0.1, 0.6, 0.1,
            0, 0.1, 0.6
        ),
        nrow = 3
    )
    # this is its biased limit  
    expect_silent(
        kinship_gcta_lim <- kinship_gcta_limit(kinship)
    )
    # compare directly to standard's limit
    expect_silent(
        kinship_std_lim <- kinship_std_limit(kinship)
    )
    # off-diagonal elements must agree!
    indexes <- lower.tri( kinship_gcta_lim )
    expect_equal( kinship_gcta_lim[ indexes ], kinship_std_lim[ indexes ] )
    
    # diagonal should have this form
    kinship_mean <- mean( kinship )
    kinship_mean_j <- rowMeans( kinship )
    diag_gcta_exp <- ( diag( kinship ) - kinship_mean_j ) / ( 1 - kinship_mean )
    expect_equal( diag( kinship_gcta_lim ), diag_gcta_exp )
    
    # run through generic validation regardless
    expect_silent( popkin::validate_kinship( kinship_gcta_lim ) )
})

test_that("kinship_wg_limit works", {
    # construct a dummy kinship matrix
    kinship <- matrix(
        c(
            0.6, 0.1, 0,
            0.1, 0.6, 0.1,
            0, 0.1, 0.6
        ),
        nrow = 3
    )
    # this is its biased WG limit
    kinship_wg_lim <- kinship_wg_limit(kinship)
    
    # validate output!
    expect_silent( popkin::validate_kinship( kinship_wg_lim ) )

    # match dims (already tested to be square by validate_kinship, so just check once)
    expect_equal( nrow( kinship ), nrow( kinship_wg_lim ) )

    # zero overall off-diagonal mean
    expect_equal( mean( kinship_wg_lim[ lower.tri( kinship_wg_lim ) ] ), 0 )
})

### tests that require BED files

if (suppressMessages(suppressWarnings(require(genio)))) {

    ### SNPRelate comparisons
    
    if (suppressMessages(suppressWarnings(require(SNPRelate)))) {

        test_that( "popkinsuppl agrees with SNPRelate", {

            ### SETUP

            # write BED files only if both packages are present
            # save data to temporary location
            # output path without extension
            file_out <- tempfile('X-popkinsuppl')
            genio::write_plink( file_out, X = X )

            # now write GDS file
            file_gds <- paste0( file_out, '.gds')
            SNPRelate::snpgdsBED2GDS(
                           paste0( file_out, '.bed'),
                           paste0( file_out, '.fam'),
                           paste0( file_out, '.bim'),
                           file_gds,
                           verbose = FALSE
                       )

            # load file, for all analyses
            genofile <- SNPRelate::snpgdsOpen( file_gds )

            ### WG FST (for subpops)

            # we found that our HudsonK calculations match what they call WH in their code (but is really WG FST for subpopulations).
            # here we repeat those tests

            # their version
            fst_wg <- SNPRelate::snpgdsFst(
                                     genofile,
                                     as.factor( labs ),
                                     method = 'W&H02',
                                     verbose = FALSE
                                 )
            
            # our version
            fst_hudsonk <- fst_hudson_k(
                X,
                labs
            )$fst

            expect_equal( fst_hudsonk, fst_wg$Fst )

            ### WG KINSHIP
            
            # get WG kinship estimate with their own package
            kinship_wg <- SNPRelate::snpgdsIndivBeta(
                                         genofile,
                                         with.id = FALSE,
                                         inbreeding = FALSE,
                                         verbose = FALSE
                                     )

            # get WG estimates via popkin and popkinsuppl
            # first get popkin estimates
            kinship_popkin <- popkin::popkin( X )
            # then just shift down with limit formula!
            kinship_wg_popkin <- kinship_wg_limit( kinship_popkin )
            # compare
            expect_equal( kinship_wg_popkin, kinship_wg )

            ### CLEANUP

            # close GDS file
            SNPRelate::snpgdsClose(genofile)

            # remove GDS file
            #file.remove( file_gds )
            genio:::delete_files_generic( file_out, 'gds' )
            # remove BED files
            genio::delete_files_plink( file_out )
        })
    }

}
