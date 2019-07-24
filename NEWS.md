# 2019-05-22 - popkinsuppl 1.0.0.9000

* First public version on GitHub!

# 2019-05-23 - popkinsuppl 1.0.1.9000

* Changes shared by all FST methods (WC, Hudson K, and Hudson Pairwise).
  * Improved speed and memory usage
  * Added `m` option, required when genotypes are passed via a function (rare)
  * Renamed option `indexes_keep` to `ind_keep`, corrected documentation to say it filters individuals (incorrectly said it filtered loci).

# 2019-05-30 - popkinsuppl 1.0.2.9000

* Added `kinship_std_limit`

# 2019-07-04 - popkinsuppl 1.0.3.9000

* Fixed error in `kinship_std` description of `mean_of_ratios` parameter.
  Now it correctly says that `FALSE` is the default value.

# 2019-07-24 - popkinsuppl 1.0.4.9000

* Adapted package to use new `solve_m_mem_lim` from the `popkin` package, which generalizes some of the calculations that used to be part of the memory estimation/control code in `get_mem_lim_m_kinship_std` and `get_mem_lim_m_WC`.

