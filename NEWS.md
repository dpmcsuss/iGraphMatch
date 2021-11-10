
# iGraphMatch 2.0.0

## Breaking Changes

* All `graph_match_` functions have been consolidated into one function called `gm` that takes a method argument indicating which of the methods to use. This will break pretty much any code but was instituted to improve flexibility and extensibility going forward. The returned values for these functions should be relatively consistent with the older values except that they'll have a lot more functionality through the `graphMatch` class.
* `get_perm` is now called `get_perm_mat` and is more flexible, taking a graph match object.
* Removed `match_report`, `edge_match_info`, `match_plot_igraph`, `match_plot_matrix`, `matched_adjs`. Instead use methods associated with `graphMatch` class such as `summary`, `plot`, `%*%` etc.
* `sample_correlated_rdpg` is now `sample_correlated_rdpg_pair`.
* `lapjv` and `lapmod` functions are now accessed solely through `do_lap`.
* Most start functions are removed and can be accessed through `init_start`.
* row_cor, row_diff, row_perm_stat have been removed and are now accessed through `best_matches`.
* `sample_*_with_junk` functions are removed. Instead use the normal `sample_*` functions with the `ncore` parameter.
* `gm_expand_when_stuck` is removed and can be accessed by specifying `ExpandWhenStuck = TRUE` when calling `gm` with `method = "percolation"`.

### List of all removed functions

* `bari_start`, `edge_match_info`, `get_perm`, `gm_indefinite`, `graph_match_ExpandWhenStuck`, `graph_match_FW`, `graph_match_IsoRank`, `graph_match_PATH`, `graph_match_Umeyama`, `graph_match_convex`, `graph_match_percolation`, `lapjv`, `lapmod`, `match_plot_igraph`, `match_plot_matrix`, `match_report`, `matched_adjs`, `rds_from_sim_start`, `rds_perm_bari_start`, `rds_sinkhorn_start`, `row_cor`, `row_diff`, `row_perm_stat`, `rperm`, `sample_correlated_gnp_pair_w_junk`,`sample_correlated_sbm_pair_w_junk`,

## New features

* Implemented `graphMatch` class for inspecting and using matches.
* `as.graphMatch` can convert data.frames and matrices to the `graphMatch` class.
* The `method` argument in `gm` can also take a function so that you can leverage the `graph_match` class and error checking of `gm`.
* `best_matches` will compute precision at k if given the true labels.


## Minor improvements and fixes

* Lots of small fixes.
* Convenience function `largest_cc` for finding largest connected component.
* Improved documentation
* The `check_graph` function is now exported to prepare graphs outside of `gm` functions.
* `matrix_list`s can now be named lists.


# iGraphMatch 1.0.1

* Added a `NEWS.md` file to track changes to the package.
* Fixed errors on Solaris.
