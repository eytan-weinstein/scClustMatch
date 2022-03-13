#' get_first_clustering
#'
#' Generates an S4 \code{clustering} class storing information related to
#' the clustering of reads aligned to the user's Reference Genome #1.
#'
#' @param first_directory The file path to the user's local directory
#' containing the (1) barcodes/cell IDs, (2) features/gene IDs, and (3) count
#' matrix outputted by the alignment of reads to Reference Genome #1.
#'
#' @param reference_resolution The resolution at which to cluster the reads.
#'
#' @param dims (optional) The dimensions of reduction with which to cluster the
#' reads. If not selected by the user, this parameter is set to 30 by default.
#'
#' @return An S4 \code{clustering} class storing information related to the
#' user's Reference Genome #1 clustering.
#'
#' @import methods
#' @import Seurat
#'
#' @export

get_first_clustering <- function(first_directory, reference_resolution, dims) {
    single_cell_data <- list.files(first_directory)
    if (!(("barcodes.tsv.gz" %in% single_cell_data) &
    ("features.tsv.gz" %in% single_cell_data) &
    ("matrix.mtx.gz" %in% single_cell_data) &
    dir.exists(first_directory))) {
        stop("Invalid first_directory. Ensure that your directory contains each
        of: \"barcodes.tsv.gz\", \"features.tsv.gz\", and \"matrix.mtx.gz\" as
        outputs from some valid read alignment pipeline.")
    }
    if (reference_resolution < 0 | reference_resolution > 10) {
        stop("Invalid reference_resolution. The resolution for this clustering
        must be some real number between 0 and 10.")
    }
    if (missing(dims)) {
        dims <- 30
    }
    if (dims < 1 | dims > 100 | dims %% 1 != 0) {
        stop("Invalid dims. The dimensions of reduction for this clustering
        must be some integer between 1 and 100.")
    }

    # Clustering with the SCTransform pipeline in Seurat
    # (see https://satijalab.org/seurat/articles/sctransform_vignette.html)

    first_clustering <- Seurat::Read10X(first_directory)
    first_clustering <-
    Seurat::CreateSeuratObject(counts = first_clustering)
    first_clustering <- Seurat::SCTransform(first_clustering)
    first_clustering <- Seurat::RunPCA(first_clustering)
    first_clustering <-
    Seurat::FindNeighbors(first_clustering, dims = 1:dims)
    first_clustering <- Seurat::FindClusters(first_clustering,
    resolution = reference_resolution)
    first_clustering <- Seurat::RunUMAP(first_clustering, dims = 1:dims)
    first_clustering <- Seurat::RunTSNE(first_clustering)
    first_clustering <- methods::new("clustering",
    Seurat_object = first_clustering,
    clustering_resolution = reference_resolution, dims = dims)
    return(first_clustering)
}