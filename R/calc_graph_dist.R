#' Internal function for graph-based distance
#'
#' This function is internal and is called by various S4 methods after method
#' specific pre-processing of data.
#'
#' @inheritParams BiocNeighbors::findKNN
#' @param \dots Further arguments to pass to \code{\link[BiocNeighbors]{findKNN}}.
#' @importFrom graph graphBAM
#' @importFrom BiocNeighbors findKNN
#' @importFrom RBGL johnson.all.pairs.sp
#' @importFrom methods setMethod
#' @importFrom BiocParallel SerialParam
#'

.calc_graph_dist <- function(X, k, BNINDEX, BNPARAM, BPPARAM = SerialParam(), ...) {
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
#' different algorithms to find k-nearests neighbors. Use this argument to
#' change the algorithm. Using a pre-computed index will save when multiple KNN
#' search are performed on the same X. If \code{BNINDEX} is specified, then
#' X does not need to be specified and any value specified for X will be ignored.
#'
#' The \code{BNPARAM} argument is used for setting parameters for KNN search
#' algorithms, such as the kind of distance metric used.
#'
#' Only one of \code{BNINDEX} and \code{BNPARAM} is needed to determine the
#' algorithm used, and if both are supplied, they must specify the same algorithm.
#' If both are missing, then the KmKNN algorithm will be used.
#'
#' @inheritParams .calc_graph_dist
#' @rdname calc_graph_dist
#' @param x A SingleCellExperiment object, a \code{seurat} object, or a matrix
#' containing expression values for each gene (row) in each cell (column). The
#' matrix can be a sparse matrix (\code{\link[Matrix]{dgCMatrix}} or other
#' sparse matrix classes from the \code{Matrix} package). The data in this matrix
#' should be normalized. If the cells are in rows, then set
#' \code{transposed = TRUE} when calling this function.
#' @return A dense square numeric matrix with n cells columns and rows. The
#' entry at ith row and jth column represents the normalized shortest path
#' length between vertex i and vertex j.
#' @export
setGeneric("calc_graph_dist", function(x, k, BNINDEX, BNPARAM,
                                       BPPARAM = SerialParam(),
                                       ...) {
  standardGeneric("calc_graph_dist")
})

#' @rdname calc_graph_dist
#' @param transposed Logical, whether the matrix has cells in rows rather than
#' in columns.
#' @param n.pcs Number of principal components to use if KNN search is to be
#' done in PCA space. If \code{NA}, which is the default, the full matrix as
#' specified in x will be used for KNN search. If a positive integer, then
#' the number specified will be the number of top principal components used.
#' @param irlba.args Named list of arguments to be passed to
#' \code{\link[irlba]{prcomp_irlba}}, such as whether to scale and center the
#' data prior to PCA.
#' @export
setMethod("calc_graph_dist", "ANY",
          function(x, k, BNINDEX, BNPARAM, BPPARAM = SerialParam(),
                   transposed = FALSE,
                   n.pcs = NA,
                   irlba.args = list(), ...) {
            if (!transposed) x <- t(x)
            if (!is.na(n.pcs)) {
              if (n.pcs < 0) {
                stop("n.pcs must be NA or a positive integer.")
              }
              irlba.args$x <- x
              irlba.args$retx <- TRUE
              irlba.args$n <- n.pcs
              x_use <- do.call(prcomp_irlba, irlba.args)$x
              out <- .calc_graph_dist(x_use, k, BNINDEX = BNINDEX, BNPARAM = BNPARAM,
                                      BPPARAM = BPPARAM, ...)
            } else {
              out <- .calc_graph_dist(x, k, BNINDEX = BNINDEX, BNPARAM = BNPARAM,
                                      BPPARAM = BPPARAM, ...)
            }
            return(out)
          })

#' @rdname calc_graph_dist
#' @param assay.use A string specifying which assay to use, defaults to
#' \code{logcounts}, namely log1p normalized data.
#' @param use.dimred The low dimensional representation of the data to use for
#' KNN search. Should be a string to use to access dimension reductions in
#' \code{\link[SingleCellExperiment]{reducedDim}}, If \code{NA}, as default, the
#' full data as specified in \code{assay.use} will be used. This argument can
#' also be a numeric index of the position of the desired dimension reduction
#' result. If not \code{NA}, then \code{assay.use} will be ignored and the
#' low dimensional representation will be used for KNN search.
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment assay
#' @export
setMethod("calc_graph_dist", "SingleCellExperiment",
          function(x, k, BNINDEX, BNPARAM, BPPARAM = SerialParam(),
                   assay.use = "logcounts",
                   use.dimred = NA, ...) {
            if (!is.na(use.dimred)) {
              out <- .calc_graph_dist(reducedDim(x, use.dimred), k,
                                      BNINDEX, BNPARAM, BPPARAM, ...)
            } else {
              out <- .calc_graph_dist(assay(x, i = assay.use), k,
                                      BNINDEX, BNPARAM, BPPARAM, ...)
            }
            return(out)
          })

#' @rdname calc_graph_dist
#' @importFrom Seurat GetAssayData GetDimReduction
#' @importClassesFrom Seurat seurat
#' @inheritParams Seurat::GetAssayData
#' @param reduction.type Type of dimension reduction to use for KNN search. If
#' \code{NA}, then the full data as specified by \code{assay.type} and \code{slot}
#' will be used. Otherwise \code{assay.type} and \code{slot} will be ignored,
#' and the dimension reduction specified by \code{reduction.type} and \code{slot.dr}
#' will be used instead.
#' @param slot.dr A string specifying the slot within the dimension reduction to
#' use for KNN search, defaults to \code{"cell.embeddings"}.
#' @export
setMethod("calc_graph_dist", "seurat",
          function(x, k, BNINDEX, BNPARAM, BPPARAM = SerialParam(),
                   assay.type = "RNA", slot = "data",
                   reduction.type = NA, slot.dr = "cell.embeddings",
                   ...) {
            if (!is.na(reduction.type)) {
              out <- .calc_graph_dist(GetDimReduction(x, reduction.type, slot.dr),
                                      k, BNINDEX, BNPARAM, BPPARAM, ...)
            } else {
              out <- .calc_graph_dist(GetAssayData(x, assay.type, slot),
                                      k, BNINDEX, BNPARAM, BPPARAM, ...)
            }
            return(out)
          })
