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
#' @param tD_loc Transpose of the graph-based distance matrix between pairs of
#' locations generated from function \code{\link{calc_graph_dist}} in this
#' package. Since the k-nearest neighbor graph is directed, the matrix is not
#' necessarily symmetric.
#' @param T_tmp A matrix for probabilistic assignment of cells (in rows) to
#' locations (in columns). This matrix will be updated in gradient descent.
#'
#' @return A numeric matrix with cells in rows and locations in columns.
tens_gr_square_loss <- function(D_cell, tD_loc, T_tmp) {
  out <- -D_cell %*% T_tmp %*% tD_loc
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
#' @param tol Tolerance. When \code{alpha != 1}, this function uses projected gradient
#' descent to find the probabilistic assignment of cells to locations. When the
#' change in objective value is less than a tolerance, then the grdient descent
#' iteration will stop.
#' @param maxiter Maximum number of iterations in projected gradient descent.
#' @param verbose Whether to display progress during projected gradient descent.
#' @param check_every Check change in objective value once in how many iterations
#' in gradient descent. For instance, if 10 is passed to this argument, then
#' check once every 10 iterations. You may not want to check every single
#' iteration since checking involves computing the Frobenius norm of a potentially
#' large matrix, which can be time consuming. However, if you do not check
#' frequently enough, then this function will keep on running after the actual
#' change in objective is less than the tolerance until the next check.
#' @return A matrix with cells in rows and locations in columns.
#' @importFrom Barycenter Sinkhorn
gw_assign <- function(D_cell, D_loc, D_cell_loc, alpha, p, q,
                      loss_fun = "square", epsilon = 5e-4, maxiter = 1000,
                      tol = sqrt(.Machine$double.eps), verbose = TRUE,
                      check_every = 10) {
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
    D_loc <- t(D_loc) # Don't transpose in every single iteration
    while (err > tol && cpt < maxiter) {
      Tprev <- T_tmp
      tens <- tens_gr_square_loss(D_cell, D_loc, T_tmp)
      tens_all <- (1 - alpha) * tens + alpha * Dcl_norm
      T_tmp <- Sinkhorn(p, q, tens_all, epsilon, maxIter = 1000)$Transportplan
      if (cpt %% check_every == 0) {
        err <- norm(T_tmp - Tprev, "F")
        if (verbose) {
          cat("Iteration ", cpt, ", delta: ", err, "\n", sep = "")
        }
      }
      cpt <- cpt + 1
    }
    if (err > tol) {
      warning("Change in objective is greater than tolerance when maxiter is exhausted.")
    }
  }
  return(T_tmp)
}

# To do: Unit test on toy example
