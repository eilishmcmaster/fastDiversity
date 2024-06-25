fastDiversity: A simple and efficient population genetics analysis package.
=================================================================================

This is the repository for the latest version of the ```fastDiversity``` package. It contains the latest bug fixes, as well as the newest functionality.

This package was developed because although other packages exist for these methods, they require specific data structures, the writing of specific file types, and in some cases are no longer maintained. Here I reproduce the most essential functions for basic population genomic analysis in a format that is compatible with any SNP data in biallelic 0,1,2 (aa, aA, AA) format (e.g. DArTseq), and saves time compared to other methods. 


Installing fastDiversity
=======================


# Stable version
# ---------------
# 
# A stable version of the package can be installed from ```CRAN``` as follows:
# 
# ```s
# install.packages("fastDiversity")
# ```

Development version
-------------------

```s
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install the package
devtools::install("https://github.com/eilishmcmaster/fastDiversity")
```