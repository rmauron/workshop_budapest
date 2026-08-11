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
undesirable_genes <- grep(pattern = c("^MT-|^HB|^RP|MALAT1"), x = rownames(object), value = T)

undesirable_genes <- setdiff(undesirable_genes, c("HBEGF", "HBP1")) # we want to keep that Growth factor and transcription factor

undesirable_genes
#
#
#
#
#
#
genes_to_keep <- setdiff(rownames(object), undesirable_genes)
length(rownames(object)) - length(genes_to_keep)
#
#
#
#
#
#
object <- SubsetSTData(object, features = genes_to_keep)
object
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
#| eval: false
#| 
library(rtracklayer)

gtf <- import("/srv/home/10x_references/refdata-gex-GRCh38-2020-A/genes/genes.gtf")

colnames(mcols(gtf))

unique(gtf$type)
# [1] gene           transcript     exon           CDS            start_codon   
# [6] stop_codon     UTR            Selenocysteine

annot <- gtf[gtf$type %in% c("gene")] |>
  mcols() |>
  as.data.frame() |>
  dplyr::select(gene_id, gene_name, gene_type) |>
  dplyr::distinct()

saveRDS(annot, "./data/GRCh38-2020-A_gene_biotypes.rds")
#
#
#
#
annot <- readRDS("./data/GRCh38-2020-A_gene_biotypes.rds")
head(annot)
genes_query <- rownames(object)

#
#
#
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
