library(Seurat)
library(ggplot2)

#Details
#Fetch Data:
  
#The function retrieves expression data for the specified features from the Seurat object.
#Each cell is associated with its corresponding identity (cell type).

#Calculate Statistics:
  
#For each cell type, the function calculates:
#AvgExp: Average expression of each feature, computed as the mean of expm1-transformed values.
#PctExp: Percentage of cells in the group expressing the feature, defined as the percentage of cells with expression values > 0.

#Filter Features:
  
#For each feature, the function examines its average expression across cell types.
#If the highest average expression occurs in the target cell type and exceeds 1.4 times the second-highest average expression, the feature is retained as a relevant marker.

#Return Markers:
 
#The function returns a list of features (markers) that meet the criteria.
#Usage: get_relevant_markers(object, assay, features, target_cell_type)
get_relevant_markers <- function(object, assay, features, target_cell_type){
  cells <- unlist(x = CellsByIdentities(object = object, cells = colnames(object[[assay]]), 
                                        idents = Idents(object)))
  
  # Fetch data for selected features
  data.features <- FetchData(object = object, vars = features, 
                             cells = cells)
  
  # Add cell identities to the data frame
  data.features$id <- Idents(object = object)[cells]
  
  # Ensure `id` is a factor
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  
  # Initialize result data frame
  result.df <- data.frame()
  
  # Iterate over unique identities to calculate avg.exp and pct.exp
  for (ident in levels(data.features$id)) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(100 * mean(x > 0))
    })
    
    # Combine into a data frame with identifiers
    ident.df <- data.frame(
      Feature = names(avg.exp),
      CellType = ident,
      AvgExp = avg.exp,
      PctExp = pct.exp
    )
    
    # Append to the result data frame
    result.df <- rbind(result.df, ident.df)
  }
  
  markers <- c()
  for (feature in features){
    # Filter results_df for the target cell type
    cell_type_data <- result.df[result.df$Feature == feature, ]
    
    cell_type_data <- cell_type_data[order(cell_type_data$AvgExp, decreasing = T),]
    if(nrow(cell_type_data) > 1){
      if (cell_type_data$CellType[1] == target_cell_type){
        
        second_highest_avg_exp <- sort(cell_type_data$AvgExp, decreasing = TRUE)[2]
        threshold <- second_highest_avg_exp * 1.4
        # Find features with AvgExp greater than the threshold
        if (cell_type_data$AvgExp[1]> threshold){
          markers <- c(markers, cell_type_data$Feature[1])
        }
      }
    }

  }
  return(markers)
}


object <- readRDS("Seuratobject.rds")
object <- UpdateSeuratObject(object)
assay <- "RNA"

markers <- read.csv("cell-markers-database.csv")
target_cell_type <- "Astrocyte"


markers <- markers[markers$cell_type == target_cell_type,]

features <- markers$marker

cell_markers <- get_relevant_markers(object, "RNA", features, target_cell_type)

DotPlot(object, cell_markers) + scale_size(range = c(3,8))
