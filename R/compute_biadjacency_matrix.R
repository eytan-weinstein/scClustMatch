#' compute_biadjacency_matrix
#'
#' Computes a biadjacency matrix between the clusters in the user's Reference
#' Genome #1 clustering and those in their Reference Genome #2 clustering. Each
#' entry i,j in that matrix is the Jaccard similarity score between cluster i in
#' the user's Reference Genome #1 clustering and cluster j in their Reference
#' Genome #2 clustering. Jaccard similarity is computed on the basis of cell IDs
#' shared between clusters. Cluster labels are replaced to be Seurat-compliant.
#'
#' @param first_clustering The S4 \code{clustering} class storing information
#' related to the user's Reference Genome #1 clustering.
#'
#' @param second_clustering The S4 \code{clustering} class storing information
#' related to the user's Reference Genome #2 clustering.
#'
#' @return A biadjacency matrix between the clusters in \code{first_clustering}
#' and those in \code{second_clustering}.
#'
#' @import methods
#' @import Seurat

compute_biadjacency_matrix <- function(first_clustering, second_clustering) {
    cells_by_cluster <- remove_unique_cells(first_clustering, second_clustering) # nolint
    similarity_scores <- c()
    for (cluster_1 in cells_by_cluster[[1]]) {
        for (cluster_2 in cells_by_cluster[[2]]) {
            similarity_scores <- c(similarity_scores,

            # Calculating Jaccard similarity between each pair of clusters

            ((length(intersect(cluster_1, cluster_2))) /
            (length(union(cluster_1, cluster_2)))))
        }
    }
    biadjacency_matrix <- matrix(similarity_scores,
    byrow = TRUE,
    nrow = get_number_of_clusters(first_clustering),
    ncol = get_number_of_clusters(second_clustering))

    # Replacing cluster labels to be Seurat-compliant

    rownames(biadjacency_matrix) <- c(0:(nrow(biadjacency_matrix) - 1))
    colnames(biadjacency_matrix) <- c(0:(ncol(biadjacency_matrix) - 1))
    return(biadjacency_matrix)
}