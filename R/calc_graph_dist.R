#' Compute graph-based distance among cells or locations
#'
#' Since Euclidean distance and Pearson correlation cannot capture the true
#' geometry of non-linear low dimensional manifolds, a graph-based distance
#' is used instead in \code{novoSpaRc}. This function first computes a k-nearest
#' neighbor graph among cells or locations. Then it infers the shortest pairwise
#' path lengths on the graph for cells and locations, resulting in a graph-based
#' distance matrix, which is then used for the optimal transport reconstruction
#' of locations of gene expression.
#'
#' Whlie the Python implementation of this package uses the Floyd Warshall
#' algorithm to find the shortest path between vertices in the graph, this
#' function uses the Johnson algorithm, which is more efficient for sparse
#' graphs. Let \eqn{V} denote the number of vertices in the graph, and \eqn{E} number of
#' edges. The Floyd Warshall algorithm has complexity \eqn{O(V^3)}, while the
#' Johnson algorithm has complexity \eqn{O(V^2 \log(V) + VE)}. We expect k-nearest
#' neighbor graphs to be sparse since k is usually much smaller than the number
#' of vertices, so the number of edges is much smaller than in the complete
#' graph, which is \eqn{V(V-1)} in directed graphs.
#'
#' The \code{BPPARAM} argument is used for parallel computing in k-nearrest
#' neighbor search. For instance, use \code{BPPARAM = MulticoreParam(3)} for
#' using 3 threads in shared memory computing.
#'
#' The \code{BNINDEX} argument is for precomputed index information for
#' different algorithms to find k-nearests neighbors. Using a pre-computed index
#' will save when multiple KNN search are performed on the same x.
#' If \code{BNINDEX} is to be specified, then x should be missing.
#'
#' The \code{BNPARAM} argument is used for setting parameters for KNN search
#' algorithms, such as the kind of distance metric used.
#'
#' Only one of \code{BNINDEX} and \code{BNPARAM} is needed to determine the
#' algorithm used, and if both are supplied, they must specify the same algorithm.
#' If both are missing, then the KmKNN algorithm will be used.
#'
#' @inheritParams BiocNeighbors::findKNN
#' @param X Numeric matrix (can be sparse) with genes in rows and cells in
#' columns. If genes are in columns, then set \code{transposed = TRUE}.
#' @param k Number of nearest neighbor when constructing k-nearest neighbor graph.
#' @param \dots Further arguments to pass to \code{\link[BiocNeighbors]{findKNN}}.
#' @param transposed Logical, whether the matrix has cells in rows rather than
#' in columns. Defaults to \code{FALSE}.
#' @param BPPARAM An object of \code{\link[BiocParallel]{BiocParallelParam}}
#' class for parallelization. See details.
#' @return A dense square numeric matrix with n cells columns and rows. The
#' entry at ith row and jth column represents the normalized shortest path
#' length between vertex i and vertex j.
#' @importFrom graph graphBAM
#' @importFrom BiocNeighbors findKNN
#' @importFrom RBGL johnson.all.pairs.sp
#' @importFrom BiocParallel SerialParam
#' @importFrom Matrix t
#' @export
calc_graph_dist <- function(X, k, BNINDEX, BNPARAM, BPPARAM = SerialParam(),
                            transposed = FALSE, ...) {
  if (!transposed) X <- Matrix::t(X)
  fknn_args <- c(k = k, BPPARAM = BPPARAM, list(...))
  if (missing(BNINDEX)) {
    fknn_args$X <- X
    if (!missing(BNPARAM)) fknn_args$BNPARAM <- BNPARAM
  } else {
    fknn_args$BNINDEX <- BNINDEX
    if (!missing(BNPARAM)) fknn_args$BNPARAM <- BNPARAM
  }
  fknn_args$get.distance <- FALSE
  knn <- do.call(findKNN, fknn_args)
  # Convert to graph
  g <- graphBAM(data.frame(from = as.vector(row(knn$index)),
                           to = as.vector(knn$index),
                           weight = 1),
                edgemode = "directed")
  # Shortest path
  sp <- johnson.all.pairs.sp(g)
  # Normalize
  sp_max <- max(sp[!is.infinite(sp)])
  sp[is.infinite(sp)] <- sp_max
  sp <- (sp - mean(sp)) / sp_max
  return(sp)
}
