#' If something is a matrix
#'
#' To allow users to use sparse matrices from the Matrix package and easily test
#' if an object is a matrix, either base R matrix or Matrix matrix.
#'
#' @param X Object to be tested.
#' @return Logical, whether X is a matrix of some sort.
#' @importFrom methods is
isMatrix <- function(X) {
  is.matrix(X) | is(X, "Matrix")
}

#' Plot spatial gene expressions
#'
#' This function can be used to plot expression of one or more genes in spatial
#' locations. This function can be used to visualize predictions from this
#' package, and gene expression in general as well if locations are known.
#' Plotly-based interactive graphics can be generated, but by default, ggplot2
#' will be used.
#'
#' @param gene_expressions Numeric vector or matrix with genes in rows and
#' locations in columns. Matrices with genes in columns instead are also permitted,
#' with \code{transposed = TRUE}.
#' @param locations Data frame or numeric matrix with 1 to 3 columns specifying
#' locations, assumed to have x in the first column, y in the second, and z in
#' the 3rd.
#' @param symmetry A string for the column name in \code{locations} that specifies
#' an axis about which symmetry is present. This is useful if due to symmetry,
#' only half of the locations are in \code{locations} as the other half would be
#' the same anyway. However, if all locations are present in \code{locations}
#' even though the actual tissue does have symmetry, this argument can be left
#' missing.
#' @param transposed Logical, whether the matrix in \code{gene_expressions} has
#' genes in columns and locations in rows. Defaults to \code{FALSE}.
#' @param interactive Logical, whether the plot returned should be interactive.
#' If \code{TRUE}, then a Plotly-based html widget will be returned. Otherwise
#' a \code{ggplot} object. Defaults to \code{FALSE}. When there are 3 columns
#' in \code{locations}, this argument will be ignored and an interactive plot
#' will always be returned.
#' @param n_col Number of columns in multi-faceted plot when multiple genes are
#' plotted at once.
#' @param pt_size Numeric, point size.
#' @param alpha Numeric, transparency of points.
#'
#' @return If \code{interactive = FALSE} and the plot is less than 3D, a \code{ggplot2}
#' object. Otherwise a Plotly object.
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap scale_color_viridis_c coord_equal
#' @importFrom plotly plot_ly toWebGL add_markers subplot layout
#' @importFrom dplyr filter
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @export
plot_spatial_expression <- function(gene_expressions, locations, symmetry,
                                    transposed = FALSE, interactive = FALSE,
                                    n_col = 2, pt_size = 1, alpha = 1) {
  if (!is.atomic(gene_expressions) && !isMatrix(gene_expressions)) {
    stop("gene_expressions must be a numeric vector or matrix.")
  }
  if (is.matrix(gene_expressions) && !is.numeric(gene_expressions)) {
    stop("gene_expression must be numeric")
  }
  if (!isMatrix(locations) && !is.data.frame(locations)) {
    stop("locations must be a matrix or a data frame.")
  }
  if (isMatrix(locations)) {
    if (!is.numeric(locations)) {
      stop("location matrix must be numeric.")
    }
    locations <- as.data.frame(locations)
  }
  if (!transposed) gene_expressions <- t(gene_expressions)
  if (isMatrix(gene_expressions) && is.null(colnames(gene_expressions))) {
    stop("Genes must be named in gene_expressions.")
  }
  n <- if (isMatrix(gene_expressions)) nrow(gene_expressions) else length(gene_expressions)
  if (n != nrow(locations)) {
    stop("Number of locations must be the same in gene_expressions and locations.")
  }

  # Rename columns of locations into xyz
  ndims <- ncol(locations)
  if (ndims > 3) {
    warning("locations has more than 3 dimensions. The first 3 columns are used.")
    ndims <- 3
    locations <- locations[,1:3]
  }
  if (!missing(symmetry) && !symmetry %in% names(locations)) {
    stop("symmetry does not match any column in locations.")
  }
  nms <- letters[24:26]
  if (!missing(symmetry)) symmetry <- nms[match(symmetry, names(locations))]
  for (i in seq_len(ndims)) {
    names(locations)[i] <- nms[i]
  }

  # Make the data frame for plotting
  if (isMatrix(gene_expressions)) {
    ngenes <- ncol(gene_expressions)
    gene_expressions <- as.data.frame(gene_expressions)
    df <- cbind(locations, gene_expressions)
  } else {
    ngenes <- 1
    df <- locations
    df$value <- gene_expressions
  }
  if (!missing(symmetry)) {
    df2 <- df
    df2[[symmetry]] <- -df[[symmetry]]
    df <- rbind(df, df2)
  }
  # Avoid no visible binding for global variagle x, y, z, gene, value
  x <- y <- z <- value <- gene <- NULL
  # 2D case
  if (ndims < 3) {
    if (ngenes > 1) {
      df <- df %>%
        gather(key = "gene", "value", -x, -y)
    }
    p <- ggplot(df, aes(x, y, color = value)) +
      geom_point(size = pt_size, alpha = alpha) +
      scale_color_viridis_c() +
      coord_equal()
    if (ngenes > 1) {
      p <- p +
        facet_wrap(~gene, ncol = n_col)
    }
    if (interactive) p <- toWebGL(p)
  } else {
    # 3D case
    if (!(ngenes > 1)) {
      p <- plot_ly(df, x = ~x, y = ~y, z = ~z, color = ~value) %>%
        add_markers(marker = list(size = 3 * pt_size, sizemode = "diameter",
                                  opacity = alpha))
    } else {
      df <- df %>%
        gather(key = "gene", "value", -x, -y, -z)
      genes_use <- colnames(gene_expressions)
      plts <- list()
      for (i in seq_along(genes_use)) {
        plts[[i]] <- df %>%
          filter(gene == genes_use[i]) %>%
          plot_ly(x = ~x, y = ~y, z = ~z, color = ~value,
                  scene = paste0("scene", i)) %>%
          add_markers(marker = list(size = 3 * pt_size, sizemode = "diameter",
                                    opacity = alpha))
      }
      # Partition canvas
      n_row <- ceiling(ngenes / n_col)
      part_x <- (0:n_col) / n_col
      part_y <- (n_row:0) / n_row
      scenes <- list()
      annots <- list()
      plt_n <- 1
      for (i in seq_len(n_row)) {
        for (j in seq_len(n_col)) {
          scenes[[plt_n]] <- list(domain = list(x = c(part_x[j], part_x[j + 1]),
                                                y = c(part_y[i + 1], part_y[i])))
          annots[[plt_n]] <- list(x = part_x[j] + 0.5/n_col,
                                  y = part_y[i] - 0.1/n_row,
                                  text = genes_use[plt_n],
                                  showarrow = FALSE,
                                  xref = 'paper', yref = 'paper')
          plt_n <- plt_n + 1
        }
      }
      names(scenes) <- c("scene", paste0("scene", 2:ngenes))
      scenes$p <- subplot(plts)
      scenes$annotations <- annots
      p <- do.call(layout, scenes)
    }
  }
  return(p)
}
