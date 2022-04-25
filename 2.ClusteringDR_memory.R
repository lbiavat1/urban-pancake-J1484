rm(list = ls())

library(readxl)
library(CATALYST)
library(flowCore)
library(tidyverse)
library(diffcyt)
library(ggpubr)

### Set PrimaryDirectory
dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

### Set 'input' directory
setwd(PrimaryDirectory)
if(!dir.exists("memory_data"))
  dir.create("memory_data")
setwd("memory_data/")
InputDirectory <- getwd()
if(!dir.exists("fcs"))
  dir.create("fcs")
setwd("fcs/")
fcsDir <- getwd()
fcsDir
setwd(PrimaryDirectory)

### Set 'metadata' directory
setwd(PrimaryDirectory)
if(!dir.exists("memory_metadata"))
  dir.create("memory_metadata")
setwd("memory_metadata/")
MetaDirectory <- getwd()
setwd(PrimaryDirectory)

### Create output directory
if(!dir.exists("memory_output"))
  dir.create("memory_output", showWarnings = FALSE)
setwd("memory_output")
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

plotAbundances(sce, k = "meta12", by = "cluster_id", group_by = "condition")

plotExprHeatmap(sce, features = type_markers(sce), k = "meta12", by = "cluster_id",
                fun = "mean", scale = "last", bars = TRUE, perc = TRUE)
CATALYST::plotDR(sce, dr = "UMAP", color_by = "meta12", facet_by = "condition") +
  geom_density2d(binwidth = 0.006, colour = "black")
ggsave(file.path(OutputDirectory, "UMAP_BMvsPB.pdf"), plot = last_plot())

CATALYST::plotDR(sce, dr = "UMAP", color_by = c("CD69", "TIGIT", "PD1"), facet_by = "condition") +
  geom_density2d(binwidth = 0.006, colour = "black")
ggsave(file.path(OutputDirectory, "UMAP_Dysfxn.pdf"), plot = last_plot())

annotation_table <- as.data.frame(cbind(c(1:12), c(1:12)))

colnames(annotation_table) <- c("meta12", "Clusters")
annotation_table$Clusters <- factor(annotation_table$Clusters)
sce <- mergeClusters(sce, k = "meta12", 
                     table = annotation_table, id = "cluster_annotation", overwrite = TRUE)
sce$cluster_annotation <- cluster_ids(sce, "cluster_annotation")


FDR_cutoff <- 0.05
ei <- sce@metadata$experiment_info
plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "condition")
ggsave(file.path(OutputDirectory, "BMvsPB_Abundance.pdf"), plot = last_plot())
pdf(file.path(OutputDirectory, "BMvsPB_Heatmap.pdf"))
plotExprHeatmap(sce, features = type_markers(sce), k = "cluster_annotation", 
                by = "cluster_id", scale = "last", bars = TRUE, perc = TRUE)
dev.off()


# DA using edgeR
design <- createDesignMatrix(ei,
                             cols_design = c("condition", "patient_id"))
contrast <- createContrast(c(0, 1, rep(0,10)))

nrow(contrast) == ncol(design)

out_DA <- diffcyt(sce,
                  experiment_info = ei, design = design, contrast = contrast,
                  analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = FALSE,
                  transform = FALSE, normalize = FALSE, min_cells = 50, min_samples = 8)
da <- rowData(out_DA$res)
plotDiffHeatmap(sce, da, top_n = 12, all = TRUE, fdr = FDR_cutoff)


