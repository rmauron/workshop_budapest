#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: setup
#| include: false

.libPaths("/home/raphael.mauron/miniconda3/envs/VisiumHD/lib/R/library") # set project library

.libPaths() # check project library

lapply(.libPaths()[1], list.files) # check library

rm(list = ls()) # empty object loaded
#
#
#
#
#
#| label: load packages

package.list <- c("semla", "Seurat", "viridis", "dplyr", "ggplot2", "gtools", "patchwork",
                  "RColorBrewer", "harmony", "tidyverse")
invisible(lapply(package.list, function(element) {
  library(element, character.only = TRUE)
}))
rm(package.list)
#
#
#
#
#
# setwd("/home/raphael.mauron/projects/hdca_whole_embryo/code/")
# getwd()
# savepath <- "../results/"
#
#
#
#
#
set.seed(1)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
data_root_directory <- file.path("./data/spaceranger/V10F24-105_C1")

samples <- Sys.glob(paths = file.path(data_root_directory, "filtered_feature_bc_matrix.h5"))
imgs <- Sys.glob(paths = file.path(data_root_directory, "spatial", "tissue_hires_image.png"))
spotfiles <- Sys.glob(paths = file.path(data_root_directory, "spatial", "tissue_positions_list.csv"))
json <- Sys.glob(paths = file.path(data_root_directory, "spatial", "scalefactors_json.json"))
#
#
#
#
#
#
infoTable <- tibble(samples,
                    imgs,
                    spotfiles,
                    json)
#
#
#
#
#
#
#
#
#
#
object <- ReadVisiumData(infoTable)
object <- LoadImages(object)
object
#
#
#
#
#
#
#

ImagePlot(object)
#
#
#
#
#
MapFeatures(
  object = object,
  slot = "counts",
  features = "MYH6",
  pt_size = 3,
  image_use = "raw"
)
#
#
#
#
#
#
#
#
#
#
#
#
object
#
#
#
#
#
#
#
#
object <- NormalizeData(object = object)
object <- ScaleData(object = object)
#
#
#
#
#
p1 <- MapFeatures(
  object = object,
  slot = "counts", # Here we use the non-normalized & unscaled layer
  features = "MYH6",
  pt_size = 3,
  image_use = "raw"
)

p2 <- MapFeatures(
  object = object,
  slot = "data", # Here we use the normalized and sclaed layer
  features = "MYH6",
  pt_size = 3,
  image_use = "raw"
)

p <- p1 | p2
p
#
#
#
#
#
#
object <- FindVariableFeatures(object = object, nfeatures = 3000)
object <- RunPCA(object)
#
#
#
#
#
#
#
data_root_directory <- file.path("./data/spaceranger/*")
metadata <- read.csv("./data/metadata.csv")

samples <- Sys.glob(paths = file.path(data_root_directory, "filtered_feature_bc_matrix.h5"))
imgs <- Sys.glob(paths = file.path(data_root_directory, "spatial", "tissue_hires_image.png"))
spotfiles <- Sys.glob(paths = file.path(data_root_directory, "spatial", "tissue_positions_list.csv"))
json <- Sys.glob(paths = file.path(data_root_directory, "spatial", "scalefactors_json.json"))
#
#
#
infoTable <- tibble(samples,
                    imgs,
                    spotfiles,
                    json,
                    slide_id = list.files(path = "./data/spaceranger/"),
                    sample_id = list.files(path = "./data/spaceranger/"),
                    organ = metadata$organ,
                    organism = metadata$organism,
                    age = metadata$pcw_age,
                    sex = metadata$sex)

rm(data_root_directory, samples, imgs, spotfiles, json)
#
#
#
#
#
object_3sect <- ReadVisiumData(infoTable)
object_3sect <- LoadImages(object_3sect, image_height = 1931)
object_3sect

rm(infoTable)
gc()
#
#
#
#
#
#| eval: false

ImagePlot(object_3sect)
#
#
#
MapFeatures(
  object = object_3sect,
  slot = "counts",
  features = "MYH6",
  pt_size = 3,
  image_use = "raw"
)
#
#
#

#
#
#
#
