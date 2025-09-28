https://satijalab.org/seurat/articles/spatial_vignette


## Identification of Spatially Variable Features
modify the codes
```
SpatiallyVariableFeatures_workaround <- function(object, assay = "SCT", selection.method = "moransi") {
  
  # Check if object is a Seurat object
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }
  
  # Check if assay is a valid assay
  if (!assay %in% names(object@assays)) {
    stop("assay must be a valid assay")
  }
  
  # Extract meta.features from the specified object and assay
  data <- object@assays[[assay]]@meta.features
  
  # Select columns starting with the provided selection.method (e.g., 'moransi')
  moransi_cols <- grep(paste0("^", selection.method), colnames(data), value = TRUE)
  
  # Filter rows where ".spatially.variable" is TRUE AND NOT NA
  # This is the crucial bug fix (avoiding NAs)
  filtered_data <- data[
    data[[paste0(selection.method, ".spatially.variable")]] & 
      (!is.na(data[[paste0(selection.method, ".spatially.variable")]])), 
    moransi_cols
  ]
  
  # Sort filtered data by ".spatially.variable.rank" column in ascending order
  sorted_data <- filtered_data[order(filtered_data[[paste0(selection.method, ".spatially.variable.rank")]]), ]
  
  # Return row names (the gene names) of the sorted data frame
  return(rownames(sorted_data))
}

top.features <- head(SpatiallyVariableFeatures_workaround(brain, selection.method = "moransi"), 6)
SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))
```
