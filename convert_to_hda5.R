
library(Seurat)
library(hdf5r)

rna_files <- list(
  "E11" = "E11_5_filtered_gene_bc_matrices_h5.h5", 
  "E12" = "E12_5_filtered_gene_bc_matrices_h5.h5", 
  "E13" = "E13_5_filtered_gene_bc_matrices_h5.h5", 
  "E14" = "E14_5_filtered_gene_bc_matrices_h5.h5", 
  "E15" = "E15_5_filtered_gene_bc_matrices_h5.h5", 
  "E16" = "E16_5_filtered_gene_bc_matrices_h5.h5")
rna_list <- lapply(rna_files, Read10X_h5)

E11_rna <- rna_list$E11
E12_rna <- rna_list$E12
E13_rna <- rna_list$E13
E14_rna <- rna_list$E14
E15_rna <- rna_list$E15
E16_rna <- rna_list$E16

library(SingleR)
library(anndata)

E11_seurat <- CreateSeuratObject(E11_rna, project = "E11.5")
E12_seurat <- CreateSeuratObject(E12_rna, project = "E12.5")
E13_seurat <- CreateSeuratObject(E13_rna, project = "E13.5")
E14_seurat <- CreateSeuratObject(E14_rna, project = "E14.5")
E15_seurat <- CreateSeuratObject(E15_rna, project = "E15.5")
E16_seurat <- CreateSeuratObject(E16_rna, project = "E16.5")

combined <- merge(E11_seurat, y = c(E12_seurat, E13_seurat, E14_seurat, E15_seurat, E16_seurat),
                  add.cell.ids = c("E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E16.5"))

# convert to h5ad for MapMyCells
# Get list of all layers
layers <- Layers(combined)

# Extract counts from all layers and combine
counts_list <- lapply(layers[grep("counts", layers)], function(layer) {
  GetAssayData(combined, layer = layer)
})

# Combine all count matrices
counts <- do.call(cbind, counts_list)
metadata <- combined@meta.data

# Create anndata object
adata <- AnnData(
  X = t(counts),
  obs = metadata,
  var = data.frame(gene = rownames(counts), row.names = rownames(counts))
)

write_h5ad(adata, "merged_data.h5ad")