# WC version of this function (details are different compared to popkin)
get_mem_lim_m_WC <- function(n, r, m = NA, mem = NA, factor = 0.7) {
    # NOTE: in estimating the chunk size, we don't know ahead of time if there are NAs or not!
    # so these calculations assume M in getMAInt is an n*n matrix
    # if there aren't any NAs, M is a scalar and we end up underestimating memory usage (better than the other way around)

    # n must be defined now, or this doesn't work!
    if (missing(n))
        stop('`n` is required!')
    if (missing(r))
        stop('`r` is required!')
    
    # try to get total memory from the system if mem wasn't specified, so it works reasonably
    if (is.na(mem)) {
        mem <- popkin:::get_mem_lim(factor = factor) # infer a reasonable default from system!
    } else {
        # assuming mem is in GB, let's convert to bytes
        mem <- mem * GB
    }
    
    # estimating total memory usage in bytes
    # Pi = m*8+ao
    # k2ps = m*r*8+mo
    # Xi[, k2is[[k]] ] = m*(n/r)*8+mo # on average, let's be conservative: m*n*8+mo
    # (k2ps-Pi)^2 = m*r*8+mo # intermediate calculation, unnamed?
    # s2s = m*8+ao # final length
    # nus = m*8+ao # another vector
    # Xi == 1 = m*n*4+mo # intermediate matrix, one of the bigger ones!
    # hbars = m*8+ao
    # FstTsi = m*8+ao
    # FstBsi = m*8+ao

    # the general function doesn't have an "r" input, but since "n" only enters in m*n matrices and so does "r", we can hack it all to have an effective "n" that combines both
    # Have:
    # - 1.5 m*n (input int matrix, plus 1 subpop extract)
    # - 2 m*r
    # since there's no other n's this is equivalent to
    # n = 1.5 * n + 2 * r
    
    data <- popkin:::solve_m_mem_lim(
                         mem = mem,
                         n = 1.5 * n + 2 * r,
                         m = m,
                         mat_m_n = 1, # left simple to not confound with above `n` hack
                         vec_m = 6
                     )
    
    return( data$m_chunk ) # return desired value
}

