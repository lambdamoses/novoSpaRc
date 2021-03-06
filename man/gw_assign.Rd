% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gw.R
\name{gw_assign}
\alias{gw_assign}
\title{Compute probabilistic assignment of cells to locations}
\usage{
gw_assign(D_cell, D_loc, D_cell_loc, alpha, p, q, loss_fun = "square",
  epsilon = 5e-04, maxiter = 1000, tol = sqrt(.Machine$double.eps),
  verbose = TRUE, check_every = 1)
}
\arguments{
\item{D_cell}{Graph-based distance matrix between pairs of cells, generated
from function \code{\link{calc_graph_dist}} in this package.}

\item{D_loc}{Graph-based distance matrix between pairs of
locations generated from function \code{\link{calc_graph_dist}} in this
package. Since the k-nearest neighbor graph is directed, the matrix is not
necessarily symmetric.}

\item{D_cell_loc}{Distance (e.g. Euclidean) matrix between each cell and each
location in gene expression space if in situ atlas is present. This argument
can be left missing if an in situ atlas (e.g. from seqFISH, MERFISH,
or STARmap).}

\item{alpha}{Weight to give to the in situ atlas if it's present.}

\item{p}{Numeric vector, marginal distribution over the cells, by default uniform.}

\item{q}{Numeric vector, Marginal distribution over the genes, by default uniform.}

\item{loss_fun}{A string that specifies loss function for difference between
cells' pairwise distance and locations' pairwise distance. Right now only
"square", meaning square loss, is supported.}

\item{epsilon}{Entropy regularization term.}

\item{maxiter}{Maximum number of iterations in projected gradient descent.}

\item{tol}{Tolerance. When \code{alpha != 1}, this function uses projected gradient
descent to find the probabilistic assignment of cells to locations. When the
change in objective value is less than a tolerance, then the grdient descent
iteration will stop.}

\item{verbose}{Whether to display progress during projected gradient descent.}

\item{check_every}{Check change in objective value once in how many iterations
in gradient descent. For instance, if 10 is passed to this argument, then
check once every 10 iterations. You may not want to check every single
iteration since checking involves computing the Frobenius norm of a potentially
large matrix, which can be time consuming. However, if you do not check
frequently enough, then this function will keep on running after the actual
change in objective is less than the tolerance until the next check.}
}
\value{
A matrix with cells in rows and locations in columns.
}
\description{
This function uses the graph-based distances between cells and between
locations, and optionally a distance matrix (e.g. Euclidean) bewteen cells
and locations in gene expression space if an in situ atlas is present, to
compute a probabilistic assignment of cells to locations. See the Methods
section of the \href{https://www.biorxiv.org/content/biorxiv/early/2018/10/30/456350.full.pdf}{novoSpaRc paper}.
}
\details{
R uses BLAS for linear algebra, and the default BLAS that
comes with R is not optimized. Using an optimized version of BLAS, such as
OpenBLAS and Intel MKL (as in Microsoft R Open), can speed up matrix
multiplication several times compared to default BLAS. For the projected
gradient descent, a matrix multiplication is done in each iteration, and the
matrix can be large if your dataset has many cells. Therefore, an optimized
BLAS is strongly recommended. See vignette for how to make R use optimized
BLAS.
}
