
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Introduction

While single cell RNA-seq (scRNA-seq) gives us insights on biological
systems with unprecedented resolution, as tissue dissociation is
required for scRNA-seq, spatial context of gene expression is lost.
`novoSpaRc` is a way to reconstruct the spatial context of gene
expression by optimal transport. `novoSpaRc` is described in the paper
[Charting a Tissue from Single Cell
Transcriptomes](https://www.biorxiv.org/content/10.1101/456350v1). The
paper authors implemented this method in Python, which can be found
[here](https://github.com/rajewsky-lab/novosparc). This package is an R
implementation of this method.

In short, what `novoSpaRc` does is to find a probabilistic assignment of
cells to locations by aligning structural similarities between the
graphs generated for single cells in expression space and physical
space, with or without taking into account an in situ atlas that
quantifies expression of some landmark genes in spatial locations.

## Installation

This package is not yet on CRAN or Bioconductor. Please install it with

``` r
devtools::install_github("lambdamoses/novoSpaRc")
```

We also strongly recommend an optimized BLAS, as this package makes
heavy use of matrix multiplications, and the matrices can be large for
larger datasets. The default BLAS that comes with R is not optimized; an
optimized BLAS can speed up matrix multiplications several times. To use
an optimized BLAS from R, you can install [Microsoft R
Open](https://mran.microsoft.com/open), which uses the optimized [Intel
Math Kernel Library (MKL)](https://software.intel.com/en-us/mkl) on
Windows and Linux and Acceleration Framework on MacOS for BLAS. Other
optimized BLAS are [OpenBLAS](https://www.openblas.net) and
[ATLAS](http://math-atlas.sourceforge.net).

For MacOS users, you may change the version of BLAS used by R without
reinstalling R by adding a symbolic link from `libRblas.dylib` to the
desired BLAS. For example, for R installed from CRAN binary, here is how
to make R use the BLAS from Acceleration Framework, which comes with
MacOS:

``` bash
cd /Library/Frameworks/R.framework/Resources/lib
ln -sf /System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/libBLAS.dylib \
libRblas.dylib
```

For changing BLAS on Windows and Ubuntu, see [this blog
post](http://brettklamer.com/diversions/statistical/faster-blas-in-r/).
For Fedora, see [this
post](https://loveshack.fedorapeople.org/blas-subversion.html).
