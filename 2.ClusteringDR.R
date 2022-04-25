rm(list = ls())

library(readxl)
library(CATALYST)
library(flowCore)
library(tidyverse)
library(diffcyt)

### Set PrimaryDirectory
dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

### Set 'input' directory
setwd(PrimaryDirectory)
if(!dir.exists("data"))
  dir.create("data")
setwd("data/")
InputDirectory <- getwd()
if(!dir.exists("fcs"))
  dir.create("fcs")
setwd("fcs/")
fcsDir <- getwd()
fcsDir
setwd(PrimaryDirectory)

### Set 'metadata' directory
setwd(PrimaryDirectory)
if(!dir.exists("metadata"))
  dir.create("metadata")
setwd("metadata/")
MetaDirectory <- getwd()
setwd(PrimaryDirectory)

### Create output directory
if(!dir.exists("output"))
  dir.create("output", showWarnings = FALSE)
setwd("output")
OutputDirectory <- getwd()
setwd(PrimaryDirectory)

sce <- readRDS(file.path(OutputDirectory, "sce_Tcell.rds"))

pbMDS(sce, dims = c(1,2), fun = "median", features = type_markers(sce),
      size_by = TRUE, by = "sample_id")

seed <- 123456
set.seed(seed)
sce <- cluster(sce, features = "type", xdim = 10, ydim = 10, maxK = 20, 
               verbose = TRUE, seed = seed)
delta_area(sce)
# Run dimensionality reduction
n_cells <- 1000
n_events <- min(n_cells(sce))
sce <- runDR(sce, dr =  "UMAP", cells = n_cells, features = "type")

plotAbundances(sce, k = "meta8", by = "cluster_id", group_by = "condition")
ggsave(filename = file.path(OutputDirectory, "ExpVsNExp_Abundance.pdf"))
pdf(filename = file.path(OutputDirectory, "ExpVsNExp_Heatmap.pdf"))
plotExprHeatmap(sce, features = type_markers(sce), k = "meta8", by = "cluster_id",
                fun = "mean", scale = "last", bars = TRUE, perc = TRUE)
dev.off()
CATALYST::plotDR(sce, dr = "UMAP", color_by = "meta8", facet_by = "condition") +
  geom_density2d(binwidth = 0.006, colour = "black")
ggsave(filename = file.path(OutputDirectory, "UMAP_ExpVsNExp.pdf"))

annotation_table <- as.data.frame(cbind(c(1:8), c(1:8)))

colnames(annotation_table) <- c("meta8", "Clusters")
annotation_table$Clusters <- factor(annotation_table$Clusters)
sce <- mergeClusters(sce, k = "meta8", 
                     table = annotation_table, id = "cluster_annotation", overwrite = TRUE)
sce$cluster_annotation <- cluster_ids(sce, "cluster_annotation")


FDR_cutoff <- 0.05
ei <- sce@metadata$experiment_info
plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
ggsave(filename = file.path(OutputDirectory, "ExpVsNExp_Abundance.pdf"))
pdf(file.path(OutputDirectory, "ExpVsNExp_Heatmap.pdf"))
plotExprHeatmap(sce, features = type_markers(sce), k = "cluster_annotation", 
                by = "cluster_id", scale = "last", bars = TRUE, perc = TRUE)
dev.off()
# DA using edgeR
design <- createDesignMatrix(ei,
                             cols_design = c("condition"))
contrast <- createContrast(c(0, 1))

nrow(contrast) == ncol(design)

out_DA <- diffcyt(sce,
                  experiment_info = ei, design = design, contrast = contrast,
                  analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE,
                  transform = FALSE, normalize = FALSE, min_cells = 50, min_samples = 2)
da <- rowData(out_DA$res)
plotDiffHeatmap(sce, da, top_n = 8, all = TRUE, fdr = FDR_cutoff)
