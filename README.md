# popkinsuppl

`popkinsuppl` is a supplementary package containing custom implementations of the existing kinship and Fst approaches benchmarked in Ochoa and Storey (2016).
The main package, `popkin`, provides our novel estimators and is available on CRAN.

Use of the functions in `popkinsuppl` is discouraged except for comparisons, which is why we are not planning on making it available on CRAN.
This package does have documentation and tests of the key functions.

## Installation

You can install `popkinsuppl` from the GitHub repository using `devtools`:
```R
install.packages("devtools") # if needed
library(devtools)
install_github('OchoaLab/popkinsuppl')
```

## Example

### Simulate random genotypes

First we need a genotype matrix `X` and, for some of the functions, subpopulation labels `labs`.
Here we construct simple data without population structure.
The key functions are only interesting under a structured population, but will work with any data in the correct format.

```R
# dimensions of simulated data
n_ind <- 100
m_loci <- 1000
n_data <- n_ind * m_loci

# missingness rate
miss <- 0.1

# simulate ancestral allele frequencies
# uniform (0,1)
p_anc <- runif(m_loci)

# simulate some binomial data
X <- rbinom(n_data, 2, p_anc)

# sprinkle random missingness
X[ sample(X, n_data * miss) ] <- NA

# turn into a matrix
X <- matrix(X, nrow = m_loci, ncol = n_ind)

# create fake subpopulation labels
# (note there are no actual subpopulations in genotype matrix)
# k_subpops groups of equal size
k_subpops <- 10
labs <- ceiling( (1 : n_ind) / k_subpops )
```

### Standard Kinship estimates

There are two forms of what we call the "Standard Kinship" estimator, both are implemented.
The "ratio-of-means" version is more robust, having known convergence properties.
However, the "mean-of-ratios" version is much more prevalent in practice and in the literature.

```R
library(popkinsuppl)

# estimate kinship matrices

# ... ratio-of-means version
kinship_rom <- kinship_std(X)

# ... mean-of-ratios version
kinship_mor <- kinship_std(X, mean_of_ratios = TRUE)
```

### FST estimates

This package implements the Weir-Cockerham and Generalized "Hudson" FST estimators.
It also provides a function to compute a matrix of pairwise Hudson FST estimates, where every entry compares a pair of subpopulations.

```R
# estimate FST using the Weir-Cockerham formula
fst_wc_obj <- fst_wc(X, labs)
# the genome-wide FST estimate
fst_wc_obj$fst
# vector of per-locus FST estimates
fst_wc_obj$fst_loci

# estimate FST using the "Hudson" formula
fst_hudson_k_obj <- fst_hudson_k(X, labs)
# the genome-wide FST estimate
fst_hudson_k_obj$fst
# vector of per-locus FST estimates
fst_hudson_k_obj$fst_loci

# estimated pairwise FST matrix using the "Hudson" formula
fst_hudson_matrix <- fst_hudson_subpops(X, labs)
```

