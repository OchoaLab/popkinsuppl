# constant
GB <- 1024 * 1024 * 1024

# kinship_std version of this function (constants are different compared to popkin)
get_mem_lim_m_kinship_std <- function(n, m = NA, mem = NA, factor = 0.7) {
    # NOTE: in estimating the chunk size, we don't know ahead of time if there are NAs or not!
    # so these calculations assume M in getMAInt is an n*n matrix
    # if there aren't any NAs, M is a scalar and we end up underestimating memory usage (better than the other way around)

    # n must be defined now, or this doesn't work!
    if (missing(n))
        stop('`n` is required!')
    
    # try to get total memory from the system if mem wasn't specified, so it works reasonably
    if (is.na(mem)) {
        mem <- popkin:::get_mem_lim(factor = factor) # infer a reasonable default from system!
    } else {
        # assuming mem is in GB, let's convert to bytes
        mem <- mem * GB
    }

    # estimating total memory usage in bytes
    # A = n*n*8+mo
    # M = n*n*8+mo
    # nu = 8
    # m = 8
    # Xi = m*n*4+mo # initial int version
    # Pi = m*8+ao
    # Xi = m*n*8+mo # centered version is type double
    # is.na(Xi) = m*n*4+mo # ignore (NEWEST) or count once (NEW) or appears twice (OLD)
    # !is.na(Xi) = m*n*8+mo # gets turned into double when input to tcrossprod
    # M = n*n*8+mo # introduce again as we update it in addition
    # Xi = m*n*8+mo # ignore (NEWEST) or introduce again as it is updated (replacing NAs with zeroes; OLD)
    # A = n*n*8+mo # introduce again as we update it in addition
    # Pi*(1-Pi) = 2*(m*8+ao) # two temporary arrays in computing this
    # nu = 8 # introduce again as it is updated
    # m = 8  # ditto
    
    data <- popkin:::solve_m_mem_lim(
                         mem = mem,
                         n = n,
                         m = m,
                         mat_m_n = 2.5,
                         mat_n_n = 4,
                         vec_m = 3,
                         vec_n = 0
                     )
    
    return( data$m_chunk ) # return desired value
}
