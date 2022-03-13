#' find_matching_clusters
#'
#' Finds concordant ("matching") clusters betwen the user's Reference Genome #1
#' clustering and their Reference Genome #2 clustering by maximum weight
#' bipartite graph matching. See the README for this package for more details
#' regarding the mathematics of this algorithm.
#'
#' @param first_clustering The S4 \code{clustering} class storing information
#' related to the user's Reference Genome #1 clustering.
#'
#' @param second_clustering The S4 \code{clustering} class storing information
#' related to the user's Reference Genome #2 clustering.
#'
#' @return A vector with a length equal to the number of clusters in
#' \code{first_clustering}. Each element in this vector is the integer label
#' for the cluster in \code{second_clustering} which is most concordant with
#' the cluster labelled i - 1 in \code{first_clustering}, where i is the index
#' of the given vector element.
#'
#' @import methods
#' @import Seurat
#' @import igraph

find_matching_clusters <- function(first_clustering, second_clustering) {
    biadjacency_matrix <-
    compute_biadjacency_matrix(first_clustering, second_clustering)
    bipartite_graph <-
    igraph::graph.incidence(biadjacency_matrix, weighted = TRUE)
    bipartite_matching <- igraph::max_bipartite_match(bipartite_graph)
    matching_cluster_labels <-
    sapply(seq_len(get_number_of_clusters(first_clustering)),
    function(i)
        as.numeric(bipartite_matching$matching[i])
    )
    return(matching_cluster_labels)
}