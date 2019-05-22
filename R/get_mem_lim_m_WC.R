# WC version of this function (details are different compared to popkin)
get_mem_lim_m_WC <- function(m = NA, n = NA, r = NA, mem = NA, factor = 0.7) {
    # NOTE: in estimating the chunk size, we don't know ahead of time if there are NAs or not!
    # so these calculations assume M in getMAInt is an n*n matrix
    # if there aren't any NAs, M is a scalar and we end up underestimating memory usage (better than the other way around)

    # n must be defined now, or this doesn't work!
    if (is.na(n))
        stop('`n` is required!')
    if (is.na(r))
        stop('`r` is required!')
    
    # try to get total memory from the system if mem wasn't specified, so it works reasonably
    if (is.na(mem)) {
        mem <- popkin:::get_mem_lim(factor = factor) # infer a reasonable default from system!
    } else {
        # assuming mem is in GB, let's convert to bytes
        mem <- mem*1024*1024*1024
    }

    # NOTE: some new experiments suggest matrix overhead is 200 rather than 40!
    # mo <- 200 # matrix overhead
    # ao <- 40  # array overhead
    mo8 <- 25 # mo/8 shortcut
    ao8 <- 5 # ao/8 shortcut

    # estimating total memory usage in bytes
    # Pi = m*8+ao
    # k2ps = m*r*8+mo
    # Xi[, k2is[[k]] ] = m*(n/r)*8+mo # on average
    # (k2ps-Pi)^2 = m*r*8+mo # intermediate calculation, unnamed?
    # s2s = m*8+ao # final length
    # nus = m*8+ao # another vector
    # Xi == 1 = m*n*4+mo # intermediate matrix, one of the bigger ones!
    # hbars = m*8+ao
    # FstTsi = m*8+ao
    # FstBsi = m*8+ao

    # total formula:
    # mem = 6*(m*8+ao) + 2*(m*r*8+mo) + m*(n/r)*8+mo + m*n*4+mo
    # mem = 6*m*8 + 6*ao + 2*m*r*8 + 2*mo + m*(n/r)*8 + mo + m*n*4 + mo
    # mem = 6*m*8 + 2*m*r*8 + m*(n/r)*8 + m*n*4 + 6*ao + 4*mo
    # mem = m*(6*8 + 2*r*8 + (n/r)*8 + n*4) + 6*ao + 4*mo
    # mem = 8*( m*(6 + 2*r + (n/r) + n/2) + 6*ao8 + 4*mo8 )
    
    # solve for m (all else fixed):
    # get maximum m (number of SNPs) given n and the memory requested
    mc <- (mem/8 - 6*ao8 - 4*mo8) / (6 + 2*r + n*(1/2+1)) # turned 1/r -> 1, will it work?
    # mc <- (mem/8 - 6*ao8 - 4*mo8) / (6 + 2*r + n*(1/2+1/r))

    # NOTE m may be missing if X is a function, so we can't make these simplifying decisions (to balance load) without m in that case...
    if (!is.na(m)) {
        if (m < mc) {
            mc <- m # use the smaller one
        } else {
            # should "redistribute" based on number of chunks, to lower memory even more per iteration
            mc <- ceiling( m/ceiling(m/mc) ) # this lowers mc even more, balances load better
        }
    }
    
    #memAct <- 8*( mc*(6 + 2*r + n*(1 + 1/2)) + 6*ao8 + 4*mo8 ) # turned 1/r -> 1, will it work?

    return(mc) # return desired value
}

