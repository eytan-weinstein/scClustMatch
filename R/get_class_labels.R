#' get_class_labels
#'
#' Assigns class labels to the set of cells common to the user's Reference
#' Genome #1 clustering and the given Reference Genome #2 clustering; denoting
#' every cell in that set as belonging to cluster i in the user's Reference
#' Genome #1 clustering, and cluster j in the given Reference Genome #2
#' clustering.
#'
#' @param first_clustering The S4 \code{clustering} class storing information
#' related to the user's Reference Genome #1 clustering.
#'
#' @param second_clustering The S4 \code{clustering} class storing information
#' related to a given Reference Genome #2 clustering.
#'
#' @return A list of length 2; the first element is a vector of integer class
#' labels for \code{first_clustering}; the second a vector of integer class
#' labels for \code{second_clustering}.
#'
#' @import methods
#' @import Seurat

get_class_labels <- function(first_clustering, second_clustering) {

    # Obtaining set of cells common to the two clusterings

    cells_by_cluster <- remove_unique_cells(first_clustering, second_clustering) # nolint

    # Class labels for the user's Reference Genome #1 clustering

    first_labels <- unlist(lapply(seq_len(length(cells_by_cluster[[1]])),
    function(i)
        rep(i - 1, length((cells_by_cluster[[1]])[[i]]))
    ))

    # Class labels for the given Reference Genome #2 clustering; same order as
    # above

    second_labels <- sapply(seq_len(length(unlist(cells_by_cluster[[1]]))),
    function(i)
        grep(unlist(cells_by_cluster[[1]])[i], cells_by_cluster[[2]]) - 1
    )

    return(list(first_labels, second_labels))
}