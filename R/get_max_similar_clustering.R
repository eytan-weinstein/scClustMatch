#' get_max_similar_clustering
#'
#' Clusters the reads aligned to Reference Genome #2 iteratively at increasing
#' increments of \code{increment} resolutions; starting at
#' \code{init_resolution} and stopping only when the adjusted Rand index (ARI)
#' similarity measure between the user's Reference Genome #1 clustering and the
#' given Reference Genome #2 clustering has been globally maximized with an
#' error of +/- \code{increment} resolutions. See the README for this package
#' for more details regarding the mathematics of this algorithm.
#'
#' @param first_clustering The S4 \code{clustering} class storing information
#' related to the user's Reference Genome #1 clustering.
#'
#' @param nearest_neighbors The nearest neighbors graph construction computed
#' from the alignment of reads to Reference Genome #2.
#'
#' @param dims The dimensions of reduction with which to cluster those reads.
#'
#' @param init_resolution The initial resolution at which to cluster those
#' reads.
#'
#' @param increment The increment at which to cluster iteratively at increasing
#' resolutions beginning with \code{init_resolution}.
#'
#' @return An S4 \code{clustering} class storing information related to that
#' Reference Genome #2 clustering which is maximally similar to the user's
#' \code{first_clustering} with an error of +/- \code{increment} resolutions.
#'
#' @import methods
#' @import Seurat
#' @import aricode

get_max_similar_clustering <-
function(first_clustering, nearest_neighbors, dims, init_resolution,
increment) {
    max_similar_clustering <- NULL
    max_similarity_score <- -2
    curr_resolution <- init_resolution
    repeat {
        clustering <- Seurat::FindClusters(nearest_neighbors,
        resolution = curr_resolution)
        clustering <- Seurat::RunUMAP(clustering, dims = 1:dims)
        clustering <- Seurat::RunTSNE(clustering)
        clustering <- methods::new("clustering",
        Seurat_object = clustering,
        clustering_resolution = curr_resolution, dims = dims)
        curr_similarity_score <-
        compute_similarity(first_clustering, clustering)
        if (curr_similarity_score < max_similarity_score) {
            break
        }
        max_similarity_score <- curr_similarity_score
        max_similar_clustering <- clustering
        curr_resolution <- curr_resolution + increment
    }
    return(max_similar_clustering)
}