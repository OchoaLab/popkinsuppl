# used by FST approaches only
# in this case the function X is no longer supported (for better memory control and speed)
preprocess_genotype_obj_fst <- function(X, n, m = NA, loci_on_cols = FALSE, ind_keep = NULL) {
    # COPY FROM POPKIN
    # with some edits, not exact copy anymore
    # (n must not be missing in these applications)
    if (missing(X))
        stop('Genotype matrix `X` is required!')
    if (missing(n))
        stop('Number of individuals "n" is required!')
    
    # determine some behaviors depending on data type
    # first validate class and set key booleans
    isFn <- FALSE
    if (class(X) == 'function') {
        isFn <- TRUE
        if (is.na(m))
            stop('The number of loci `m` is required when `X` is a function!')
    } else if (class(X) == 'BEDMatrix') { # same as general matrix but transposed
        loci_on_cols <- TRUE # this is always imposed for this particular format!
    } else if (!is.matrix(X)) {
        stop('X has unsupported class: ', class(X))
    } 

    # extract dimensions from data (not possible for function version)
    if (!isFn) {
        # we overwrite m in this case, without checking for prior values
        if (loci_on_cols) {
            nX <- nrow(X)
            m <- ncol(X)
        } else {
            nX <- ncol(X)
            m <- nrow(X)
        }
        
        # if we're filtering individuals, then ind_keep is the one that should agree...
        if (!is.null(ind_keep)) {
            # this bizarre thing works if ind_keep is logical (nX<-sum(ind_keep)) or a vector of indexes (nX<-length(ind_keep)), or whichever other subsetting strategy is used (i.e. what about an input with negatives like -indexesKeep ?)
            nX <- length(
            ( 1 : nX )[ ind_keep ]
            )
        }

        # NOTE: for WC Fst disagreement between labs and X is fatal!
        if (n != nX)
            stop('Dimensions of labs and X disagree: ', n, ' != ', nX)
    } 

    # return values of interest
    list(
        isFn = isFn,
        m = m,
        loci_on_cols = loci_on_cols
    )
}
