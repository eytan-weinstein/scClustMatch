#' get_second_clustering
#'
#' Generates an S4 \code{clustering} class storing information related to
#' the clustering of reads aligned to the user's Reference Genome #2. The
#' resolution of this clustering is optimized within an error of +/- 0.01 so as
#' to maximize its similarity to the user's Reference Genome #1 clustering.
#'
#' @param second_directory The file path to the user's local directory
#' containing the (1) barcodes/cell IDs, (2) features/gene IDs, and (3) count
#' matrix outputted by the alignment of reads to Reference Genome #2.
#'
#' @param first_clustering The S4 \code{clustering} class storing information
#' related to the user's Reference Genome #1 clustering.
#'
#' @return An S4 \code{clustering} class storing information related to the
#' user's Reference Genome #2 clustering.
#'
#' @import methods
#' @import Seurat
#' @import aricode
#'
#' @export

get_second_clustering <- function(second_directory, first_clustering) {
    single_cell_data <- list.files(second_directory)
    if (!(("barcodes.tsv.gz" %in% single_cell_data) &
    ("features.tsv.gz" %in% single_cell_data) &
    ("matrix.mtx.gz" %in% single_cell_data) &
    dir.exists(second_directory))) {
        stop("Invalid second_directory. Ensure that your directory contains each
        of: \"barcodes.tsv.gz\", \"features.tsv.gz\", and \"matrix.mtx.gz\" as
        outputs from some valid read alignment pipeline.")
    }
    if (!(class(first_clustering) == "clustering")) {
        stop("Invalid first_clustering.")
    }

    # Inheriting dimensions of reduction from the user's Reference Genome #1
    # clustering

    dims <- get_dims(first_clustering) # nolint

    # Clustering with the SCTransform pipeline in Seurat
    # (see https://satijalab.org/seurat/articles/sctransform_vignette.html)

    second_clustering <- Seurat::Read10X(second_directory)
    second_clustering <- Seurat::CreateSeuratObject(counts = second_clustering)
    second_clustering <- Seurat::SCTransform(second_clustering)
    second_clustering <- Seurat::RunPCA(second_clustering)
    second_clustering <- Seurat::FindNeighbors(second_clustering, dims = 1:dims)

    # Finding the optimal resolution for the user's Reference Genome #2
    # clustering within an error of +/- 0.1 resolutions

    optimal_resolution_large_error <-
    get_clustering_resolution(get_max_similar_clustering(first_clustering,
    second_clustering, dims, 0, 0.1))

    # Finding the optimally resolved Reference #2 clustering within an error of
    # +/- 0.01 resolutions

    second_clustering <- get_max_similar_clustering(first_clustering,
    second_clustering, dims, optimal_resolution_large_error - 0.09, 0.01)
    return(second_clustering)
}