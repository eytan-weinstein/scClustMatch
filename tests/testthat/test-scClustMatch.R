library(scClustMatch)
library(rdrop2)

### Run this code to access the data for the example provided in the README

scClustMatch_token <- readRDS(system.file("testdata", "scClustMatch_token.rds", # nolint
package = "scClustMatch"))

woodchuck_MONAX5_directory <- tempdir() # nolint
rdrop2::drop_download(path = "woodchuck_MONAX5/barcodes.tsv.gz",
local_path = woodchuck_MONAX5_directory, dtoken = scClustMatch_token)
rdrop2::drop_download(path = "woodchuck_MONAX5/features.tsv.gz",
local_path = woodchuck_MONAX5_directory, dtoken = scClustMatch_token)
rdrop2::drop_download(path = "woodchuck_MONAX5/matrix.mtx.gz",
local_path = woodchuck_MONAX5_directory, dtoken = scClustMatch_token)
woodchuck_MONAX5_clustering <- get_first_clustering(woodchuck_MONAX5_directory, 0.2) # nolint

woodchuck_10_directory <- tempdir()
rdrop2::drop_download(path = "woodchuck_10/barcodes.tsv.gz",
local_path = woodchuck_10_directory, dtoken = scClustMatch_token,
overwrite = TRUE)
rdrop2::drop_download(path = "woodchuck_10/features.tsv.gz",
local_path = woodchuck_10_directory, dtoken = scClustMatch_token,
overwrite = TRUE)
rdrop2::drop_download(path = "woodchuck_10/matrix.mtx.gz",
local_path = woodchuck_10_directory, dtoken = scClustMatch_token,
overwrite = TRUE)
woodchuck_10_clustering <- get_second_clustering(woodchuck_10_directory,
woodchuck_MONAX5_clustering) # nolint

###

test_that("Clustering resolution set properly; first_clustering", {
  expect_equal(get_clustering_resolution(woodchuck_MONAX5_clustering), 0.2)
})

test_that("PCA dimensions set properly; first_clustering", {
  expect_equal(get_dims(woodchuck_MONAX5_clustering), 30)
})

test_that("Correct number of clusters assigned; first_clustering", {
  expect_equal(get_number_of_clusters(woodchuck_MONAX5_clustering),
  length(get_cells_by_cluster(woodchuck_MONAX5_clustering)))
})

test_that("Correct number of cells identified; first_clustering", {
  expect_equal(length(colnames(get_Seurat_object(woodchuck_MONAX5_clustering)@assays$RNA@counts)), # nolint
  length(unlist(get_cells_by_cluster(woodchuck_MONAX5_clustering))))
})

test_that("Clustering resolution selected properly; second_clustering", {
  expect_equal(round(compute_similarity(woodchuck_MONAX5_clustering,
  woodchuck_10_clustering), 2), 0.75)
})

test_that("Clustering resolution set properly; second_clustering", {
  expect_equal(get_clustering_resolution(woodchuck_10_clustering), 0.27)
})

test_that("PCA dimensions set properly; second_clustering", {
  expect_equal(get_dims(woodchuck_10_clustering), 30)
})

test_that("Correct number of clusters assigned; second_clustering", {
  expect_equal(get_number_of_clusters(woodchuck_10_clustering),
  length(get_cells_by_cluster(woodchuck_10_clustering)))
})

test_that("Correct number of cells identified; second_clustering", {
  expect_equal(length(colnames(get_Seurat_object(woodchuck_10_clustering)@assays$RNA@counts)), # nolint
  length(unlist(get_cells_by_cluster(woodchuck_10_clustering))))
})

test_that("Correct number of cells removed from first_clustering", {
  expect_equal(length(unlist(remove_unique_cells(woodchuck_MONAX5_clustering,
  woodchuck_10_clustering)[[1]])),
  length(colnames(get_Seurat_object(woodchuck_MONAX5_clustering)@assays$RNA@counts)) - # nolint
  length(setdiff(colnames(get_Seurat_object(woodchuck_MONAX5_clustering)@assays$RNA@counts), # nolint
  colnames(get_Seurat_object(woodchuck_10_clustering)@assays$RNA@counts)))
  )
})

test_that("Correct number of cells removed from second_clustering", {
  expect_equal(length(unlist(remove_unique_cells(woodchuck_MONAX5_clustering,
  woodchuck_10_clustering)[[2]])),
  length(colnames(get_Seurat_object(woodchuck_10_clustering)@assays$RNA@counts)) - # nolint
  length(setdiff(colnames(get_Seurat_object(woodchuck_10_clustering)@assays$RNA@counts), # nolint
  colnames(get_Seurat_object(woodchuck_MONAX5_clustering)@assays$RNA@counts)))
  )
})

test_that("Clusters matched correctly", {
  expect_equal(find_matching_clusters(woodchuck_MONAX5_clustering,
  woodchuck_10_clustering), c(1, 2, 0, 3, 4, 5, 6, 7, 8, 10, 9))
})