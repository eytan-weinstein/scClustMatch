#' compare_clusterings
#'
#' Generates an HTML summary report comparing the user's Reference Genome #1
#' clustering to their Reference Genome #2 clustering.
#'
#' @param first_clustering The S4 \code{clustering} class storing information
#' related to the user's Reference Genome #1 clustering.
#'
#' @param second_clustering The S4 \code{clustering} class storing information
#' related to the user's Reference Genome #2 clustering.
#'
#' @return An HTML summary report comparing \code{first_clustering} and
#' \code{second_clustering}. See the README for this package for more details
#' regarding the figures generated here.
#'
#' @import methods
#' @import utils
#' @import Seurat
#' @import aricode
#' @import igraph
#' @import rmarkdown
#' @import ggplot2
#' @import gplots
#' @import ggalluvial
#' @import ggrepel
#' @import knitr
#'
#' @export

compare_clusterings <- function(first_clustering, second_clustering) {
    if (!(class(first_clustering) == "clustering")) {
        stop("Invalid Reference Genome #1 clustering.")
    }
    if (!(class(second_clustering) == "clustering")) {
        stop("Invalid Reference Genome #2 clustering.")
    }
    utils::browseURL(rmarkdown::render(system.file("Rmd",
    "compare_clusterings.Rmd",
    package = "scClustMatch"),
    params = list(first_clustering = first_clustering,
    second_clustering = second_clustering)))
}