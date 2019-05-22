preprocess_genotype_obj <- function(X, n, loci_on_cols = FALSE, indexes_keep = NULL) {
    # COPY FROM POPKIN
    # with some edits, not exact copy anymore
    # (n must not be missing in these applications)
    if (missing(n))
        stop('Required number of individuals "n" missing!')
    
    # determine some behaviors depending on data type
    # first validate class and set key booleans
    isFn <- FALSE
    if (class(X) == 'function') {
        isFn <- TRUE
    } else if (class(X) == 'BEDMatrix') { # same as general matrix but transposed
        loci_on_cols <- TRUE # this is always imposed for this particular format!
    } else if (!is.matrix(X)) {
        stop('X has unsupported class: ', class(X))
    } 

    # extract dimensions from data (not possible for function version)
    if (isFn) {
        m <- NA # have to define as NA to pass to getMemLimM below
    } else {
        if (loci_on_cols) {
            nX <- nrow(X)
            m <- ncol(X)
        } else {
            nX <- ncol(X)
            m <- nrow(X)
        }
        # if we're filtering, then indexes_keep is the one that should agree...
        if (!is.null(indexes_keep)) {
            # this bizarre thing works if indexes_keep is logical (nX<-sum(indexes_keep)) or a vector of indexes (nX<-length(indexes_keep)), or whichever other subsetting strategy is used (i.e. what about an input with negatives like -indexesKeep ?)
            nX <- length(
            ( 1 : nX )[ indexes_keep ]
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
