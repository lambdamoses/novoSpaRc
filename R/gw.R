#' Calculate square loss term in gradient of objective function
#'
#' Projected gradient descent is used to numerically solve the optimization
#' problem of finding a probabilistic assignment of cells to locations. The
#' gradient of the loss function has a term \eqn{L(C_1, C_2) \otimes T}. This
#' function calculates that term when \eqn{L(a,b) = \frac 1 2 (a - b)^2}, namely
#' the square loss.
#'
#' @param D_cell Graph-based distance matrix between pairs of cells, generated
#' from function \code{\link{calc_graph_dist}} in this package.
#' @param D_loc Graph-based distance matrix between pairs of locations,
#' generated from function \code{\link{calc_graph_dist}} in this package.
#' @param T_tmp A matrix for probabilistic assignment of cells (in rows) to
#' locations (in columns). This matrix will be updated in gradient descent.
#'
#' @return A numeric matrix with cells in rows and locations in columns.
#' @importFrom Rfast mat.mult transpose
tens_gr_square_loss <- function(D_cell, D_loc, T_tmp) {
  # Rfast::mat.mult for parallel, faster matrix multiplication
  out <- -mat.mult(mat.mult(D_cell, T_tmp), transpose(D_loc))
  out - min(out)
}

#' Compute probabilistic assignment of cells to locations
#'
#' This function uses the graph-based distances between cells and between
#' locations, and optionally a distance matrix (e.g. Euclidean) bewteen cells
#' and locations in gene expression space if an in situ atlas is present, to
#' compute a probabilistic assignment of cells to locations. See the Methods
#' section of the \href{https://www.biorxiv.org/content/biorxiv/early/2018/10/30/456350.full.pdf}{novoSpaRc paper}.
#'
#' @inheritParams tens_gr_square_loss
#' @param D_cell_loc Distance (e.g. Euclidean) matrix between each cell and each
#' location in gene expression space if in situ atlas is present. This argument
#' can be left missing if an in situ atlas (e.g. from seqFISH, MERFISH,
#' or STARmap).
#' @param alpha Weight to give to the in situ atlas if it's present.
#' @param p Numeric vector, marginal distribution over the cells, by default uniform.
#' @param q Numeric vector, Marginal distribution over the genes, by default uniform.
#' @param loss_fun A string that specifies loss function for difference between
#' cells' pairwise distance and locations' pairwise distance. Right now only
#' "square", meaning square loss, is supported.
#' @param epsilon Entropy regularization term.
#' @param tol Tolerance. When \code{alpha != 1}, this function uses gradient
#' descent to find the probabilistic assignment of cells to locations. When the
#' change in objective value is less than a tolerance, then the grdient descent
#' iteration will stop.
#' @param maxiter Maximum number of iterations in gradient descent.
#' @param verbose Whether to display progress during gradient descent.
#'
#' @return A matrix with cells in rows and locations in columns.
#' @importFrom Barycenter Sinkhorn
gw_assign <- function(D_cell, D_loc, D_cell_loc, alpha, p, q,
                      loss_fun = "square", epsilon = 5e-4, maxiter = 1000,
                      tol = sqrt(.Machine$double.eps), verbose = TRUE) {
  if (epsilon < 0) stop("epsilon must not be negative.")
  # Initialize
  T_tmp <- tcrossprod(p, q)
  cpt <- 1
  err <- 1
  # Normalize D_cell_loc if present
  if (missing(D_cell_loc)) {
    if (alpha > 0) {
      stop("Weight is given to D_cell_loc, which is missing")
    } else Dcl_norm <- 0
  } else {
    if (alpha == 0) {
      warning("D_cell_loc is ignored because alpha = 0.")
      Dcl_norm <- 0
    } else {
      Dcl_norm <- D_cell_loc / max(D_cell_loc)
    }
  }
  p <- matrix(p, ncol = 1)
  q <- matrix(q, ncol = 1)
  if (alpha == 1) {
    T_tmp <- Sinkhorn(p, q, T_tmp, epsilon, maxIter = 1000)$Transportplan
  } else {
    # Gradient descent
    if (loss_fun != "square") {
      message("Only square loss is currently supported. Using square loss.")
    }
    while (err > tol && cpt < maxiter) {
      Tprev <- T_tmp
      tens <- tens_gr_square_loss(D_cell, D_loc, T_tmp)
      tens_all <- (1 - alpha) * tens + alpha * Dcl_norm
      T_tmp <- Sinkhorn(p, q, tens_all, epsilon, maxIter = 1000)$Transportplan
      if (cpt %% 10 == 0) {
        err <- norm(T_tmp - Tprev, "F")
        if (verbose) {
          cat("Iteration", cpt, ", Err:", err, "\n", sep = " ")
        }
      }
      cpt <- cpt + 1
    }
    if (err > tol) {
      warning("Error is greater than tolerance when maxiter is exhausted.")
    }
  }
  return(T_tmp)
}

# To do: Unit test on toy example
