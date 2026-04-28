# ============================================================
# Title: Cross-species integration of mouse pDC and human pDC
# Author: Alina Ulezko Antonova
# Date: 2024-05-04
# Steps this code follows:
#   - Convert mouse pDC genes to human orthologues
#   - Rebuild a mouse pDC Seurat object using human gene symbols
#   - Integrate mouse and human pDC using Seurat CCA
#   - Save processed objects and publication-quality figures
# ============================================================

# ----------------------------
# Reproducibility settings
# ----------------------------
set.seed(1234)

# ----------------------------
# Libraries
# ----------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(biomaRt)
})

# ----------------------------
# File paths
# ----------------------------
mouse_dcs_rds <- "/storage1/fs1/mcolonna/Active/Patrick/Analysis/003/DCs.combined.rds" 
human_pdc_rds <- "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/objects/all_pDC_SCT_Clean_030524.rds"
mouse_pred_rds <- "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/objects/Mouse_pDC_with_Predictions.rds"

out_mouse_human_rds <- "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/objects/mouse_pDC_hGenes_050424.rds"
out_cca_rds <- "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/objects/pdc_cca_crossSpecies_050324.rds"

plot_dir <- "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/plots"

# ----------------------------
# Sanity checks
# ----------------------------
stopifnot(file.exists(mouse_dcs_rds))
stopifnot(file.exists(human_pdc_rds))
stopifnot(dir.exists(plot_dir))

# ----------------------------
# Helper: mouse -> human orthologues
# ----------------------------
convert_mouse_to_human_genes <- function(mouse_genes) {
  human <- useMart(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl",
    host = "https://dec2021.archive.ensembl.org/"
  )

  mouse <- useMart(
    biomart = "ensembl",
    dataset = "mmusculus_gene_ensembl",
    host = "https://dec2021.archive.ensembl.org/"
  )

  map <- getLDS(
    attributes = c("mgi_symbol"),
    filters = "mgi_symbol",
    values = mouse_genes,
    mart = mouse,
    attributesL = c("hgnc_symbol"),
    martL = human,
    uniqueRows = TRUE
  )

  colnames(map) <- c("Mouse", "Human")
  map <- unique(map)
  map <- map[map$Mouse != "" & map$Human != "", ]
  map
}

# ----------------------------
# Load mouse dendritic cells
# ----------------------------
dcs <- readRDS(mouse_dcs_rds)

# Inspect broad clusters
DimPlot(dcs, label = TRUE)

# Subset pDC and pDC-like populations
pdc <- subset(dcs, idents = c("0", "7"))
pdc <- RunUMAP(pdc, dims = 1:20)
DimPlot(pdc)

# ----------------------------
# Map mouse genes to human orthologues
# ----------------------------
mouse_genes <- rownames(pdc)
gene_map <- convert_mouse_to_human_genes(mouse_genes)

# Keep only mapped genes
counts_mouse <- as.data.frame(GetAssayData(pdc, assay = "RNA", slot = "counts"))
counts_mouse$Mouse <- rownames(counts_mouse)

merged_counts <- merge(gene_map, counts_mouse, by = "Mouse", all.x = FALSE, all.y = FALSE)

# Resolve duplicate human gene symbols reproducibly
merged_counts$Human <- make.unique(merged_counts$Human)
rownames(merged_counts) <- merged_counts$Human
merged_counts$Human <- NULL
merged_counts$Mouse <- NULL

# Remove any duplicated rows after merging
merged_counts <- merged_counts[!duplicated(rownames(merged_counts)), ]

# ----------------------------
# Create mouse object in human gene space
# ----------------------------
mouse_pdc_human <- CreateSeuratObject(counts = t(as.matrix(merged_counts)))

mouse_pdc_human[["percent.mt"]] <- PercentageFeatureSet(mouse_pdc_human, pattern = "^MT-")
mouse_pdc_human <- SCTransform(mouse_pdc_human, vars.to.regress = "percent.mt", verbose = FALSE)
mouse_pdc_human <- RunPCA(mouse_pdc_human, verbose = FALSE)
mouse_pdc_human <- FindNeighbors(mouse_pdc_human, dims = 1:30)
mouse_pdc_human <- FindClusters(mouse_pdc_human, resolution = 0.1)
mouse_pdc_human <- RunUMAP(mouse_pdc_human, dims = 1:30)

# Add metadata from original pDC object
mouse_meta <- pdc[[]]
common_cells <- intersect(colnames(mouse_pdc_human), rownames(mouse_meta))
mouse_pdc_human <- AddMetaData(mouse_pdc_human, metadata = mouse_meta[common_cells, , drop = FALSE])

# Optional annotation for plotting
mouse_pdc_human$species <- "Mouse"

# Save processed mouse object
saveRDS(mouse_pdc_human, out_mouse_human_rds)

# ----------------------------
# Load human pDC object
# ----------------------------
human_pdc <- readRDS(human_pdc_rds)
Idents(human_pdc) <- human_pdc@meta.data$annotation_coarse

# Downsample to match mouse cell number
human_pdc_sub <- subset(human_pdc, downsample = 11188)
human_pdc_sub$species <- "Human"

# ----------------------------
# Harmonize genes before integration
# ----------------------------
common_genes <- intersect(rownames(mouse_pdc_human), rownames(human_pdc_sub))
stopifnot(length(common_genes) > 0)

mouse_pdc_human <- subset(mouse_pdc_human, features = common_genes)
human_pdc_sub <- subset(human_pdc_sub, features = common_genes)

# ----------------------------
# Prepare metadata for plotting
# ----------------------------
mouse_pdc_human$annotation_coarse_species <- paste("Mouse", mouse_pdc_human$annotation_coarse, sep = "_")
human_pdc_sub$annotation_coarse_species <- paste("Human", human_pdc_sub$annotation_coarse, sep = "_")

if ("Tissue" %in% colnames(mouse_pdc_human@meta.data)) {
  mouse_pdc_human$tissue_species <- paste("Mouse", mouse_pdc_human$Tissue, sep = "_")
}
if ("organ" %in% colnames(human_pdc_sub@meta.data)) {
  human_pdc_sub$tissue_species <- paste("Human", human_pdc_sub$organ, sep = "_")
}

# ----------------------------
# Cross-species CCA integration
# ----------------------------
pdc_cca <- RunCCA(object1 = human_pdc_sub, object2 = mouse_pdc_human)
pdc_cca <- FindNeighbors(pdc_cca, reduction = "cca", dims = 1:20)
pdc_cca <- FindClusters(pdc_cca, reduction = "cca")
pdc_cca <- RunUMAP(pdc_cca, reduction = "cca", dims = 1:20)

# ----------------------------
# Publishable plots
# ----------------------------
p_species <- DimPlot(
  pdc_cca,
  group.by = "species",
  cols = c("Human" = "#FF3D7FFF", "Mouse" = "#088BBEFF")
) & NoAxes() & NoLegend()

ggsave(
  filename = file.path(plot_dir, "pdc_cca_species.png"),
  plot = p_species,
  dpi = 700,
  width = 6,
  height = 4
)

p_markers <- FeaturePlot(
  pdc_cca,
  features = c("TLR7", "CX3CR1", "MKI67"),
  order = TRUE,
  split.by = "species"
) & NoAxes()

ggsave(
  filename = file.path(plot_dir, "pdc_cca_markers.png"),
  plot = p_markers,
  dpi = 700,
  width = 6,
  height = 8
)

# ----------------------------
# Save integrated object
# ----------------------------
saveRDS(pdc_cca, out_cca_rds)

# ----------------------------
# Optional inspecition plots
# ----------------------------
Idents(pdc_cca) <- pdc_cca$species
pdc_cca_mouse <- subset(pdc_cca, idents = "Mouse")
DimPlot(pdc_cca_mouse, group.by = "annotation_coarse", cols = "red")
