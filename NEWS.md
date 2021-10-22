
# iGraphMatch 2.0


* Implemented `graphMatch` class for inspecting and using matches.
* All `graph_match` functions have been consolidated into one function called `gm` that takes a method argument indicating which of the methods to use.
* The `method` argument in `gm` can also take a function so that you can leverage the `graph_match` class and error checking of `gm`.
* `get_perm` is now called `make_perm` and `get_perm` is a function that take a `graphMatch` object

# iGraphMatch 1.0.1

* Added a `NEWS.md` file to track changes to the package.
* Fixed errors on Solaris.
