#' Predict gene expression at locations
#'
#' With the probabilistic assignment of cells to locations from optimal transport,
#' we can predict gene expression at each location based on the gene count matrix.
#' Let \eqn{X} denote the gene count matrix with genes in rows and cells in
#' columns, and \eqn{T} denote the probabilistic assignment of cells to locations,
#' with cells in rows and locations in columns. Gene expression at each location
#' can be predicted by \eqn{XT}. This function will use this to predict gene
#' expression. It will also scale the prediction so the mean gene expression
#' across locations will match the mean gene expression among cells in scRNA-seq.
#'
#' @inheritParams calc_graph_dist
#' @param cell_loc Matrix that probabistically assigns cells to locations, made
#' from the function \code{\link{gw_assign}}.
#' @param scale Logical, whether predictions should be scaled so their means are
#' the same as in single cell RNA-seq data supplied in \code{x}.
#' Defaults to \code{TRUE}.
#' @return A numeric matrix with genes in rows and locations in columns.
#' @importFrom Matrix rowMeans
#' @export
predict_expr_loc <- function(X, cell_loc, transposed = FALSE,
                             scale = TRUE) {
  if (transposed) X <- t(X)
  out <- X %*% cell_loc
  if (scale) {
    means <- Matrix::rowMeans(X)
    means_pred <- Matrix::rowMeans(out)
    out <- diag(means / means_pred) %*% out
  }
  if (!is.null(rownames(X))) rownames(out) <- rownames(X)
  return(out)
}

#' Expression of 84 landmark genes in the Drosophila embryo
#'
#' The [Berkeley Drosophila Transcription Network Project (BDTNP)](http://bdtnp.lbl.gov:8080/Fly-Net/)
#' quantified expression of landmark genes at single cell resolution in
#' _Drosophila melanogaster_ embryos at different stages of development. This
#' dataset has BDTNP data for 84 landmark genes at late stage 5. The order of
#' cells in this dataset corresponds to the order of cells in the \code{\link{locations}}
#' dataset.
#'
#' @format A numeric matrix with genes in columns and cells in rows, 84 columns
#' and 3039 rows.
#' @source \url{http://bdtnp.lbl.gov:8080/Fly-Net/bioimaging.jsp?w=releases}
"bdtnp"

#' Locations of cells in the Drosophila embryo
#'
#' This dataset has the x, y, and z spatial coordinates for cells in an
#' archetypical late stage 5 Drosophila embryo. The order of cells here is the
#' same as in the \code{\link{bdtnp}} dataset. This is only half of the embryo,
#' which is symmetric about the y axis here.
#'
#' @format A data frame with 3 columns, xcoord, ycoord, and zcoord, and 3039 rows.
#' @source \url{http://bdtnp.lbl.gov:8080/Fly-Net/bioimaging.jsp?w=releases}
"locations"
