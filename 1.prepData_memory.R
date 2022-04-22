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


# prep and read in metadata

# # Step 1
#
# filenames <- list.files(InputDirectory, pattern = ".csv")
# write_csv(as.data.frame(filenames), file = file.path(MetaDirectory, "sample.details.csv"))

# filenames <- read_csv(file.path(MetaDirectory, "sample.details.csv"))
# filenames
# exp_data <- read_excel(file.path(MetaDirectory, "J1484_SampleFlow_LAST-EXP.xlsx"))
# exp_data

# # Step 2
# 
# metadata <- left_join(filenames, exp_data, by = c("Position")) %>%
#   write_csv(file = file.path(MetaDirectory, "metadata.csv"))

# # Step 3
# 
# metadata <- read_csv(file = file.path(MetaDirectory, "metadata.csv"))
# 
# patient_id <- as_tibble(unique(metadata$MRN)) %>% rename(., MRN = value) %>%
#   rownames_to_column(.) %>%
#   mutate(PT_ID = paste("PT", rowname, sep = "_")) %>%
#   select(MRN, PT_ID)
# 
# metadata <- metadata %>% mutate(Disease = "AML") %>%
#   rename(., Filename = filenames) %>%
#   rename(., Sample = "Barcode [Parent]") %>%
#   rename(., origin = 'Specimen Source') %>%
#   left_join(., patient_id, by = c("MRN")) %>%
#   mutate(Tissue = ifelse(grepl("PBL",origin), "PB", "BM")) %>%
#   select(Filename, Sample, PT_ID, MRN, Disease, Tissue, Expanded) %>%
#   arrange(PT_ID)
# #
# #
# metadata %>% write_csv(file = file.path(MetaDirectory, "tidy_metadata.csv"))



metadata <- read_csv(file = file.path(MetaDirectory, "tidy_metadata.csv"))
metadata

sample_md <- metadata
sample_md

CSVfiles <- sample_md$Filename

# convert csv files to fcs
csvTofcs <- function(file.names, dest){
  # create an empty list to start
  DataList <- list()

  for(file in file.names){
    tmp <- read_csv(file.path(file))
    file <- gsub(".csv", "", file)
    DataList[[file]] <- tmp
  }
  rm(tmp)

  filenames <- names(DataList)
  head(DataList)

  # convert csv to fcs

  for(i in c(1:length(filenames))){
    data_subset <- DataList[i]
    data_subset <- data.table::rbindlist(as.list(data_subset))
    file_name <- names(DataList)[i]

    metadata <- data.frame(name = dimnames(data_subset)[[2]], desc = "")

    # create FCS file metadata
    # metadata$range <- apply(apply(data_subset, 2, range), 2, diff)
    metadata$minRange <- apply(data_subset, 2, min)
    metadata$maxRange <- apply(data_subset, 2, max)


    # data as matrix by exprs
    data_subset.ff <- new("flowFrame", exprs = as.matrix(data_subset),
                          parameters = AnnotatedDataFrame(metadata))

    head(data_subset.ff)
    write.FCS(data_subset.ff, paste0(dest, "/", file_name, ".fcs"), what = "numeric")
  }
}
setwd(InputDirectory)
# csvTofcs(CSVfiles, fcsDir)

# select samples to keep for analysis

sample_md <- sample_md 

sample_md

fcsFiles <- list.files(path = fcsDir, pattern = ".fcs")
fcsFiles %in% gsub(".csv", ".fcs", sample_md$Filename)
fcsToLoad <- fcsFiles[fcsFiles %in% gsub(".csv", ".fcs", sample_md$Filename)]
# read fcs files as flowSet and add $CYT keyword
fs <- read.flowSet(files = fcsToLoad, path = fcsDir, truncate_max_range = FALSE)
fs

# create tibble sample_md
sample_md

# prep sample_md for SCE (prepData)
sample_md <- sample_md %>% mutate(Filename = gsub(".csv", ".fcs", Filename)) %>%
  rename(file_name = Filename) %>%
  rename(patient_id = PT_ID) %>%
  mutate(condition = Tissue) %>%
  mutate(sample_id = paste(Disease, patient_id, Tissue, sep = "-")) %>%
  select(file_name, patient_id, condition, sample_id)
sample_md

# create tibble panel_md
fcs_colname <- colnames(fs)
fcs_colname
antigen <- fcs_colname
antigen[10:35] <- sapply(fcs_colname[10:35], function(.) unlist(str_split(., " :: "))[2])
fluorochrome <- fcs_colname
fluorochrome[10:35] <- sapply(fcs_colname[10:35], function(.) unlist(str_split(., " :: "))[1])
marker_class <- fcs_colname
marker_class[ c(1:10, 12, 33:36)] <- "none"
marker_class[c(11, 13:32)] <- "type"
panel_md <- as_tibble(cbind(fcs_colname, antigen, fluorochrome, marker_class))
as.data.frame(panel_md)

# use CATALYST::prepData() to create SCE object from flowSet
sce <- prepData(fs, panel_md, sample_md, FACS = TRUE)
assay(sce, "exprs") <- assay(sce, "counts")

subsetSCE <- function(x, n_cells){
  cs <- split(seq_len(ncol(x)), x$sample_id)
  cs <- unlist(lapply(cs, function(.) sample(., min(n_cells, length(.)))))
  x <- x[, cs]
  return(x)
}

sub_sce <- subsetSCE(sce, 15000)


pbMDS(sub_sce, dims = c(1,2), fun = "median", features = type_markers(sce),
      size_by = TRUE, by = "sample_id")

sce <- sub_sce

p <- plotExprs(sce, features = NULL, color_by = "condition")
p$facet$params$ncol <- 9
p

n_events <- min(n_cells(sce))
n_events
n_cells(sce)
plotCounts(sce, group_by = "sample_id", color_by = "condition")

plotNRS(sce, features = type_markers(sce), color_by = "condition")



saveRDS(sce, file = file.path(OutputDirectory, "sce_Tcell.rds"))

