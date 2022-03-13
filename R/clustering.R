#' clustering
#'
#' An S4 class storing information related to a clustering of scRNAseq data.
#'
#' @slot Seurat_object The Seurat object output of a clustering pipeline run
#' in Seurat. This output--commonly saved locally using the Seurat::saveRDS
#' command--should bind the data produced by running all steps in a standard
#' clustering pipeline up to and including non-linear dimensionality reduction
#' (see "https://satijalab.org/seurat/articles/pbmc3k_tutorial.html").
#'
#' @slot clustering_resolution The resolution parameter used in that pipeline.
#' This must be some real number between 0 and 10 (inclusive).
#'
#' @slot dims The dimensions of reduction used in that pipeline. This must be
#' some integer between 1 and 100 (inclusive).
#'
#' @import methods
#' @import Seurat
#'
#' @export

clustering <- setClass(
    "clustering",
    representation(Seurat_object = "Seurat",
    clustering_resolution = "numeric",
    dims = "numeric"),
    validity = function(object) {
        if (object@clustering_resolution < 0 |
        object@clustering_resolution > 10) {
            stop("Invalid clustering_resolution. The resolution for this
            clustering must be some real number between 0 and 10 (inclusive).")
        }
        if (object@dims < 1 | object@dims > 100 |
        object@dims %% 1 != 0) {
            stop("Invalid dims. The dimensions of reduction for this clustering
            must be some integer between 1 and 100 (inclusive).")
        }
        return(TRUE)
    }
)


#' get_Seurat_object
#'
#' Gets the Seurat object associated with this \code{clustering}.
#'
#' @param object An S4 \code{clustering} class.
#'
#' @return The \code{Seurat_object} slot stored by this class.
#'
#' @import methods
#'
#' @export

setGeneric(
    "get_Seurat_object",
    function(object) standardGeneric("get_Seurat_object")
)

#' @rdname get_Seurat_object

setMethod(
    "get_Seurat_object",
    signature(object = "clustering"),
    function(object) {
        return(object@Seurat_object)
    }
)


#' get_clustering_resolution
#'
#' Gets the clustering resolution used to resolve this \code{clustering}.
#'
#' @param object An S4 \code{clustering} class.
#'
#' @return The \code{clustering_resolution} slot stored by this class.
#'
#' @import methods
#'
#' @export

setGeneric(
    "get_clustering_resolution",
    function(object) standardGeneric("get_clustering_resolution")
)

#' @rdname get_clustering_resolution

setMethod(
    "get_clustering_resolution",
    signature(object = "clustering"),
    function(object) {
        return(object@clustering_resolution)
    }
)


#' get_dims
#'
#' Gets the dimensions of reduction used in this \code{clustering}.
#'
#' @param object An S4 \code{clustering} class.
#'
#' @return The \code{dims} slot stored by this class.
#'
#' @import methods

setGeneric(
    "get_dims",
    function(object) standardGeneric("get_dims")
)

#' @rdname get_dims

setMethod(
    "get_dims",
    signature(object = "clustering"),
    function(object) {
        return(object@dims)
    }
)


#' get_number_of_clusters
#'
#' Gets the number of clusters in this \code{clustering}.
#'
#' @param object An S4 \code{clustering} class.
#'
#' @return The number of clusters in this \code{clustering}.
#'
#' @import methods
#' @import Seurat

setGeneric(
    "get_number_of_clusters",
    function(object) standardGeneric("get_number_of_clusters")
)

#' @rdname get_number_of_clusters

setMethod(
    "get_number_of_clusters",
    signature(object = "clustering"),
    function(object) {
        cluster_number <- 0
        repeat {
            test_cluster_number <- tryCatch(
                Seurat::WhichCells(object = object@Seurat_object,
                ident = cluster_number),
                error = function(e) e
            )
            if (inherits(test_cluster_number, "error")) {
                break
            }
            cluster_number <- cluster_number + 1
        }
        return(cluster_number)
    }
)


#' get_cells_by_cluster
#'
#' Returns a list of cells belonging to each cluster in this \code{clustering}.
#'
#' @param object An S4 \code{clustering} class.
#'
#' @return A list of lists; each sublist is a list of cell IDs belonging to the
#' cluster numbered i - 1 in this \code{clustering}, where i is the sublist's
#' index in the list.
#'
#' @import methods
#' @import Seurat

setGeneric(
    "get_cells_by_cluster",
    function(object) standardGeneric("get_cells_by_cluster")
)

#' @rdname get_cells_by_cluster

setMethod(
    "get_cells_by_cluster",
    signature(object = "clustering"),
    function(object) {
        cells_by_cluster <- lapply(0:(get_number_of_clusters(object) - 1),
        function(i)
            Seurat::WhichCells(object@Seurat_object, ident = i)
        )
        return(cells_by_cluster)
    }
)