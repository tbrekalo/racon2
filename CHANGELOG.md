# Changelog

Documents changes to the code base.

## racon -> racon2 refactor/cleanup

### dependencies
* deprecated
    * meson build support
        * meson was removed for lowering the const of maintance
    * dependency managment using git submodules
        * lowering the cost of maintance
        * less moving 3rd party parts in the directory
    * rvaser/thread_pool
        * switching to intel/oneapi/tbb for better scheduling
* added
    * cxxopts
        * more convinient command line parsing
    * tbrekalo/biosoup
        * uses compression for sequence data
    * tbb
        * cleaner interface and better scheduling

### Tests
    * removed unit tests for now
        * half of them were conected to CUDE which is no longer supported
        * the remaning half is broken due to the API changes
        * test relied on magic numbers to evaluate the correctnes of a polisher
            * this approach is too fragile
            * change the alignment strategy in POA or read aligner and they break
