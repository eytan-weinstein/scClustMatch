#' compute_similarity
#'
#' Computes the similarity between two different clusterings of the same
#' reads--one using Reference Genome #1 and the other using Reference Genome
#' #2--as an adjusted Rand index (ARI) similarity score.
#'
#' @param first_clustering The S4 \code{clustering} class storing information
#' related to the user's Reference Genome #1 clustering.
#'
#' @param second_clustering The S4 \code{clustering} class storing information
#' related to a given Reference Genome #2 clustering.
#'
#' @return The adjusted Rand index (ARI) similarity score between
#' \code{first_clustering} and \code{second_clustering}.
#'
#' @import methods
#' @import Seurat
#' @import aricode
#'
#' @export

compute_similarity <-
function(first_clustering, second_clustering) {
    if (!(class(first_clustering) == "clustering")) {
        stop("Invalid first_clustering.")
    }
    if (!(class(second_clustering) == "clustering")) {
        stop("Invalid second_clustering.")
    }
    class_labels <- get_class_labels(first_clustering, second_clustering) #nolint
    return(aricode::ARI(class_labels[[1]], class_labels[[2]]))
}