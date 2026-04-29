# ============================================================
# Title: Figure 7 - BPDCN Bone Marrow Integration (xV-seq)
# Author: Alina Ulezko Antonova
# Date: 2024-05-27
# Description: 
#   1. Process individual BPDCN patient Seurat objects.
#   2. Combine patients and integrate using Harmony.
#   3. Calculate mutation metadata (Mutation_Overall, Mutation_Counts).
#   4. Recreate Figure 7 plots.
# ============================================================

library(Seurat)
library(harmony)
library(ggplot2)
library(ggrepel)
library(tidyverse)

# ----------------------------
# 1. Setup and Load
# ----------------------------
base_dir <- "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/bpdcn_xvseq/seurat_objects"
plot_dir <- "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/plots"
setwd(base_dir)

# Helper to process raw patient files
process_patient <- function(file_path) {
  obj <- readRDS(file_path)
  obj <- SCTransform(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
  return(obj)
}

# Define files
patient_files <- c(
  "Pt1Dx_Seurat_Final.rds", "Pt1Rem_Seurat_Final.rds", "Pt5Dx_Seurat_Final.rds",
  "Pt9Dx_Seurat_Final.rds", "Pt10Dx_Seurat_Final.rds", "Pt10Rel_Seurat_Final.rds",
  "Pt12Dx_Seurat_Final.rds", "Pt12Rel_Seurat_Final.rds", "Pt14Dx_Seurat_Final.rds",
  "Pt15Dx_Seurat_Final.rds", "Pt16Dx_Seurat_Final.rds"
)

# Process and save objects
patient_objs <- lapply(patient_files, function(f) {
  obj <- process_patient(f)
  saveRDS(obj, file = gsub("_Seurat_Final.rds", "_proc.rds", f))
  return(obj)
})

# ----------------------------
# 2. Integrate with Harmony
# ----------------------------
pts <- merge(patient_objs[[1]], y = patient_objs[2:length(patient_objs)])
pts[["percent.mt"]] <- PercentageFeatureSet(pts, pattern = "^MT-")
pts <- SCTransform(pts, vars.to.regress = "percent.mt", verbose = FALSE)
pts <- RunPCA(pts, verbose = FALSE)
pts <- RunHarmony(pts, assay = "SCT", group.by.vars = "orig.ident", max.iter.harmony = 20)
pts <- RunUMAP(pts, reduction = "harmony", dims = 1:20)

# ----------------------------
# 3. Add Mutation Metadata
# ----------------------------
# Identify mutation status columns
metadata <- pts@meta.data
mutation_cols <- grep("TET2|RAB9A|MAP4K5|ASXL1|ACAP2|CWF19L2|DOLPP1|NOTCH1|NFIC|MALAT1|IDH2|EZH2|ETV6", 
                      colnames(metadata), value = TRUE)

# Mutation_Overall: "mutant" if any slot is "mutant"
metadata$Mutation_Overall <- ifelse(rowSums(metadata[, mutation_cols] == "mutant", na.rm = TRUE) > 0, 
                                   "mutant", "wildtype")

# Mutation_Counts: Count total distinct mutation slots flagged as "mutant"
mutation_counts <- rowSums(metadata[, mutation_cols] == "mutant", na.rm = TRUE)
metadata$Mutation_Counts <- paste("Count", mutation_counts)

pts@meta.data <- metadata

# Save integrated object
saveRDS(pts, "../Single-cell_BPDCN/BPDCN_cases_withMut.rds")

# ----------------------------
# 4. Generate Figures
# ----------------------------
# Example: UMAP mutation status
p_mut <- DimPlot(pts, group.by = "Mutation_Overall", cols = c("#5C3566FF", "#FCE94FFF"), shuffle = TRUE) + theme_void()
ggsave(file.path(plot_dir, "bpdcn_umap_mut.png"), p_mut, width = 6, height = 5, dpi = 700)

# Example: Cycling pDC signature
pts <- AddModuleScore(pts, features = list(c("CLEC4C", "CUX2")), name = "Cycling_pDC_sig")
p_sig <- FeaturePlot(pts, features = "Cycling_pDC_sig1") + theme_void()
ggsave(file.path(plot_dir, "bpdcn_sig.png"), p_sig, width = 7, height = 5, dpi = 700)
