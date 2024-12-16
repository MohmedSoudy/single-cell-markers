# Use case on Peripheral Blood Mononuclear Cell (PBMC) Analysis

## Background
The ifnb dataset is a single-cell RNA sequencing dataset that captures immune responses in PBMCs (Peripheral Blood Mononuclear Cells). 

These cells consist of:
  
- Lymphocytes: T cells, B cells, and NK cells.
- Monocytes: A subset of white blood cells.
- Dendritic Cells: Present in smaller quantities.

## Prerequisites

Ensure the following R libraries are installed:
- `Seurat`
- `SeuratData`
- `openxlsx`
- `ggplot2`

## Required Libraries
```r
library(Seurat)
library(SeuratData)
library(openxlsx)
library(ggplot2)
```

## Analysis Script
```r
# Load IFNB dataset
ifnb <- LoadData("ifnb")
ifnb <- UpdateSeuratObject(ifnb)

# Modify cell annotations
ifnb@meta.data$seurat_annotations <- as.character(ifnb@meta.data$seurat_annotations)
ifnb@meta.data$seurat_annotations[ifnb@meta.data$seurat_annotations == "CD8 T"] <- "CD8+ T cell"
ifnb <- SetIdent(ifnb, value = "seurat_annotations")

# Download cell markers database
cell_markers <- read.xlsx("http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download_files/file/Cell_marker_Human.xlsx")

# Filter markers for peripheral blood normal cells
cell_markers_blood <- cell_markers[cell_markers$tissue_type == "Peripheral blood" & cell_markers$cell_type == "Normal cell",]

# Define target cell type
target_cell_type <- "CD8+ T cell"

# Extract unique markers for target cell type
cell_markers_mt <- unique(cell_markers_blood$marker[cell_markers_blood$cell_name == target_cell_type])
source("https://raw.githubusercontent.com/MohmedSoudy/single-cell-markers/refs/heads/main/get_relevant_markers.R")
# Get relevant markers (function not shown in original script)
relevant_markers <- get_relevant_markers(ifnb, "RNA", cell_markers_mt, target_cell_type)

# Create dot plot
p <- DotPlot(ifnb, relevant_markers) + scale_size(range = c(1,8))

# Save plot
ggsave(filename = "CD8 T cell.jpeg", plot = p, width = 8, height = 5, dpi = 300)
```

![](https://raw.githubusercontent.com/MohmedSoudy/single-cell-markers/refs/heads/main/CD8%20T%20cell.jpeg)


## References
- [CellMarker Database](http://bio-bigdata.hrbmu.edu.cn/CellMarker/)
