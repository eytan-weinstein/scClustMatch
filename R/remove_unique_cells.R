#' remove_unique_cells
#'
#' Removes any cells from two given clusterings which are not shared between
#' those clusterings. This allows a similarity score between these clusterings
#' to be computed over a set of cells common to both clusterings.
#'
#' @param first_clustering The S4 \code{clustering} class storing information
#' related to the user's Reference Genome #1 clustering.
#'
#' @param second_clustering The S4 \code{clustering} class storing information
#' related to a given Reference Genome #2 clustering.
#'
#' @return A list of length 2; the first element is a list of non-unique cell
#' IDs by cluster in \code{first_clustering}; the second a list of a non-unique
#' cell IDs by cluster in \code{second_clustering}.
#'
#' @import methods
#' @import Seurat

remove_unique_cells <- function(first_clustering, second_clustering) {
    cells_by_cluster_first <- get_cells_by_cluster(first_clustering) # nolint
    cells_by_cluster_second <- get_cells_by_cluster(second_clustering) # nolint
    common_cells <-
    intersect(unlist(cells_by_cluster_first), unlist(cells_by_cluster_second))

    # Removing non-unique cells from the user's Reference Genome #1 clustering

    cells_by_cluster_first <- lapply(seq_len(length(cells_by_cluster_first)),
    function(i)
        intersect(cells_by_cluster_first[[i]], common_cells)
    )

    # Removing non-unique cells from the given Reference Genome #2 clustering

    cells_by_cluster_second <- lapply(seq_len(length(cells_by_cluster_second)),
    function(i)
        intersect(cells_by_cluster_second[[i]], common_cells)
    )

    return(list(cells_by_cluster_first, cells_by_cluster_second))
}