# popkinsuppl 1.0.0.9000 (2019-05-22)

* First public version on GitHub!

# popkinsuppl 1.0.1.9000 (2019-05-23)

* Changes shared by all FST methods (WC, Hudson K, and Hudson Pairwise).
  * Improved speed and memory usage
  * Added `m` option, required when genotypes are passed via a function (rare)
  * Renamed option `indexes_keep` to `ind_keep`, corrected documentation to say it filters individuals (incorrectly said it filtered loci).

# popkinsuppl 1.0.2.9000 (2019-05-30)

* Added `kinship_std_limit`

# popkinsuppl 1.0.3.9000 (2019-07-04)

* Fixed error in `kinship_std` description of `mean_of_ratios` parameter.
  Now it correctly says that `FALSE` is the default value.

# popkinsuppl 1.0.4.9000 (2019-07-24)

* Adapted package to use new `solve_m_mem_lim` from the `popkin` package, which generalizes some of the calculations that used to be part of the memory estimation/control code in `get_mem_lim_m_kinship_std` and `get_mem_lim_m_WC`.
* `kinship_std` now preserves individual names when present in input genotype matrix.

# popkinsuppl 1.0.5.9000 (2019-08-02)

* Internal changes to match updated `solve_m_mem_lim` from the `popkin` package (>= 1.2.6.9000).

# popkinsuppl 1.0.6.9000 (2019-12-17)

* Preemptively updated `class` usage now that matrices return a two-element array in R-devel
* Minor Roxygen-related updates

# popkinsuppl 1.0.7.9000 (2020-01-29)

* Implemented the Weir-Hill 2002 FST estimator (`fst_wh`).

# popkinsuppl 1.0.8.9000 (2020-05-07)

* `fst_wc` fixed a minor bug: sample size was half as it should have been, only noticeable in very small samples.
* Updates to `fst_wc`, `fst_wh`, `fst_hudson_k`, and `fst_hudson_pairwise`
  * Avoids NAs in output when entire subpopulations have missing data at one locus
  * Works correctly when there are singleton subpopulations (used to cause errors)
  * More clear error messages for fewer than two subpopulations
  * Expanded tests to cover these scenarios

# popkinsuppl 1.0.9.9000 (2020-05-21)

* Reimplemented `fst_hudson_pairwise` to make it much more efficient in large datasets
  * Original implementation simply called `fst_hudson_k` separately for every pair of subpopulations, which caused excessive and non-optimal file reading, which scaled increasingly poorly as the number of subpopulations increased.
  * New implementation directly constructs pairwise FST matrix, reading input genotype matrix only once.
  * The new implementation also better handles missing genotypes and is more numerically stable than before.
  * In the public Pacific dataset (Skoglund et al. 2016) it was verified to agree with the Hudson implementation in `EIGENSUITE` 7.2.1 (executable `smartpca`)

# popkinsuppl 1.0.10.9000 (2020-05-25)

* Reverted `fst_wc` sample size formula to what it was prior to version 1.0.8.9000.  Bias tests showed that the original sample size equation was best.

# popkinsuppl 1.0.11.9000 (2020-08-27)

* Function `fst_wc` has new option `FIT`, which if `TRUE` additionally includes calculations of the *total inbreeding* (FIT) next to the default *structural inbreeding* (FST) output list.

# popkinsuppl 1.0.12.9000 (2020-08-30)

* Debugged functions `fst_wc` and `fst_wh`, which previously propagated some matrix-vector products incorrectly.
  This bug resulted in small biases when subpopulation sample sizes were unequal (sample weighing was set incorrectly; there was no bug when all sample sizes were equal).
  The bug was introduced in version 1.0.8.9000 (commit b445e3e) and remained until this version.

# popkinsuppl 1.0.13.9000 (2020-09-24)

* Function `kinship_std` 
  - Added `want_M` option, which if `TRUE` returns a list containing the `kinship` matrix as well as the pairwise complete count matrix `M`.
  - Added `m_chunk_max` option (default 1000), which sets the maximum number of loci to process at the time.
	The new default behavior reduces memory usage a lot, especially on machines with large memory, without sacrificing speed.
	Original version would use a lot of memory just because it was available, which could be inconvenient when trying to run other processes, and did not result in increased speed, so it was unnecessary at best.
  - Fixed a rare bug that affected `mean_of_ratios = FALSE` version only, and only occured when `m_chunk` was smaller than `m_loci`.
    In original version `m_chunk` was as large as memory allowed, so the bug would have only occurred when `m_loci` was very large.

# popkinsuppl 1.0.14.9000 (2020-10-16)

* Extended `m_chunk_max` option (see `kinship_std` above) to functions: `fst_wc`, `fst_wh`, `fst_hudson_k`, and `fst_hudson_pairwise`.
  This reduces memory usage a lot without sacrificing speed.

# popkinsuppl 1.0.15.9000 (2021-02-23)

* Added function `kinship_wg_limit`, calculates the limit of the Weir-Goudet (WG) estimator.
  - Also obtains WG kinship estimates if input is the popkin kinship estimate instead of the true kinship matrix, which is validated in internal tests against the WG implementation from SNPRelate.

# popkinsuppl 1.0.16.9000 (2021-08-26)

- Added function `kinship_gcta_limit`.
- Updated DESCRIPTION main paragraph to mention kinship biased limit functions.
- Reformatted this `NEWS.md` slightly to improve its automatic parsing.
- Removed "LazyData: true" from DESCRIPTION (to avoid a new "NOTE" on CRAN).
