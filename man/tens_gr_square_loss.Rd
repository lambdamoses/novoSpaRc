% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gw.R
\name{tens_gr_square_loss}
\alias{tens_gr_square_loss}
\title{Calculate square loss term in gradient of objective function}
\usage{
tens_gr_square_loss(D_cell, tD_loc, T_tmp)
}
\arguments{
\item{D_cell}{Graph-based distance matrix between pairs of cells, generated
from function \code{\link{calc_graph_dist}} in this package.}

\item{tD_loc}{Transpose of the graph-based distance matrix between pairs of
locations generated from function \code{\link{calc_graph_dist}} in this
package. Since the k-nearest neighbor graph is directed, the matrix is not
necessarily symmetric.}

\item{T_tmp}{A matrix for probabilistic assignment of cells (in rows) to
locations (in columns). This matrix will be updated in gradient descent.}
}
\value{
A numeric matrix with cells in rows and locations in columns.
}
\description{
Projected gradient descent is used to numerically solve the optimization
problem of finding a probabilistic assignment of cells to locations. The
gradient of the loss function has a term \eqn{L(C_1, C_2) \otimes T}. This
function calculates that term when \eqn{L(a,b) = \frac 1 2 (a - b)^2}, namely
the square loss.
}
