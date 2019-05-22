# kinship_std version of this function (constants are different compared to popkin)
get_mem_lim_m_kinship_std <- function(m = NA, n = NA, mem = NA, factor = 0.7) {
    # NOTE: in estimating the chunk size, we don't know ahead of time if there are NAs or not!
    # so these calculations assume M in getMAInt is an n*n matrix
    # if there aren't any NAs, M is a scalar and we end up underestimating memory usage (better than the other way around)

    # n must be defined now, or this doesn't work!
    if (is.na(n))
        stop('n must be defined!')
    
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

    # total formula:
    # mem = 4*(n*n*8+mo) + (m*n*4+mo) + 2*(m*n*8+mo) + 3*(m*8+ao) + 4*8
    # mem = 4*n*n*8 + 4*mo + m*n*4 + mo + 2*m*n*8 + 2*mo + 3*m*8 + 3*ao + 4*8
    # mem = 4*n*n*8 + m*n*4 + 2*m*n*8 + 3*m*8 + 7*mo + 3*ao + 4*8
    # mem = 4*n*n*8 + m*(n*4 + 2*n*8 + 3*8) + 7*mo + 3*ao + 4*8
    
    # OLD 3 still overestimates mem...
    # mem = 4*(n*n*8+mo) + (m*n*4+mo) + 3*(m*n*8+mo) + 3*(m*8+ao) + 4*8
    # mem = 4*n*n*8 + 4*mo + m*n*4 + mo + 3*m*n*8 + 3*mo + 3*m*8 + 3*ao + 4*8
    # mem = 4*n*n*8 + m*n*4 + 3*m*n*8 + 3*m*8 + 8*mo + 3*ao + 4*8
    # mem = 4*n*n*8 + m*(n*4 + 3*n*8 + 3*8) + 8*mo + 3*ao + 4*8

    # OLD 2 still overestimates mem...
    # mem = 4*(n*n*8+mo) + 2*(m*n*4+mo) + 3*(m*n*8+mo) + 3*(m*8+ao) + 4*8
    # mem = 4*n*n*8 + 4*mo + 2*m*n*4 + 2*mo + 3*m*n*8 + 3*mo + 3*m*8 + 3*ao + 4*8
    # mem = 4*n*n*8 + 2*m*n*4 + 3*m*n*8 + 3*m*8 + 9*mo + 3*ao + 4*8
    # mem = 4*n*n*8 + m*(2*n*4 + 3*n*8 + 3*8) + 9*mo + 3*ao + 4*8

    # OLD overestimates mem...
    # mem = 4*(n*n*8+mo) + 3*(m*n*4+mo) + 3*(m*n*8+mo) + 3*(m*8+ao) + 4*8
    # mem = 4*n*n*8 + 4*mo + 3*m*n*4 + 3*mo + 3*m*n*8 + 3*mo + 3*m*8 + 3*ao + 4*8
    # mem = 4*n*n*8 + 3*m*n*4 + 3*m*n*8 + 3*m*8 + 10*mo + 3*ao + 4*8
    # mem = 4*n*n*8 + m*(3*n*4 + 3*n*8 + 3*8) + 10*mo + 3*ao + 4*8

    # solve for m (all else fixed):
    # get maximum m (number of SNPs) given n and the memory requested
    mc <- (mem/8 - 4*n*n - 7*mo8 - 3*ao8 - 4) / ((2+1/2)*n + 3)
    # OLD 3 still overestimates mem...
    # mc <- (mem/8 - 4*n*n - 8*mo8 - 3*ao8 - 4) / ((3+1/2)*n + 3)
    # OLD 2 still overestimates mem...
    # mc <- (mem/8 - 4*n*n - 9*mo8 - 3*ao8 - 4) / (4*n + 3) # OLD denom: (n*(4+1/2) + 3)
    # OLD overestimates mem...
    # mc = (mem - 4*n*n*8 - 10*mo - 3*ao - 4*8) / (3*n*4 + 3*n*8 + 3*8)
    # mc <- (mem/8 - 4*n*n - 10*mo8 - 3*ao8 - 4) / (3*(n*3/2 + 1))

    # NOTE m may be missing if X is a function, so we can't make these simplifying decisions (to balance load) without m in that case...
    if (!is.na(m)) {
        if (m < mc) {
            mc <- m # use the smaller one
        } else {
            # should "redistribute" based on number of chunks, to lower memory even more per iteration
            mc <- ceiling( m / ceiling( m / mc ) ) # this lowers mc even more, balances load better
        }
    }
    
    # memAct <- 8*(4*n*n + mc*(n*(2+1/2) + 3) + 7*mo8 + 3*ao8 + 4)
    
    return(mc) # return desired value
}
