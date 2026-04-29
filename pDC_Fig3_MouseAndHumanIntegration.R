# ============================================================
# Title: Cross-species integration of mouse pDC and human pDC
# Author: Alina Ulezko Antonova
# Date: 2024-05-04
# Steps this code follows:
#   - Convert mouse pDC genes to human orthologues
#   - Rebuild a mouse pDC Seurat object using human gene symbols
#   - Integrate mouse and human pDC using Seurat CCA
#   - Save processed objects and publication-quality figures
#
# Naming convention used here:
#   - mGenes = object in mouse gene symbols
#   - hGenes = object in human gene symbols
#   - orthologous = restricted to shared genes used for integration
# ============================================================

# ----------------------------
# Libraries used
# ----------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(biomaRt)
})

# ----------------------------
# Define input and output paths
# ----------------------------
mouse_dcs_rds <- "..."   # Processed object from GSE314567 (mouse tissue DC data)
human_pdc_rds <- "..."   # Processed object from Figure 1D (see code in pDC_Fig1B_HumanSeurat.R)

out_mouse_pdc_hGenes_rds <- "...rds" #path for your rds of mouse DCs using human genes
out_pdc_cca_rds <- "..." #path for your rds of CCA-integrated mouse and human pDC data
plot_dir <- "plots"

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

stopifnot(file.exists(mouse_dcs_rds))
stopifnot(file.exists(human_pdc_rds))

# ----------------------------
# Function to convert mouse -> human orthologues
# ----------------------------
convertMouseGeneList <- function(x) {
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

  genesV2 <- getLDS(
    attributes = c("mgi_symbol"),
    filters = "mgi_symbol",
    values = x,
    mart = mouse,
    attributesL = c("hgnc_symbol"),
    martL = human,
    uniqueRows = TRUE
  )

  genePairs <- data.frame(
    Mouse = genesV2[, 1],
    Human = genesV2[, 2],
    stringsAsFactors = FALSE
  )

  genePairs <- unique(genePairs)
  genePairs <- genePairs[genePairs$Mouse != "" & genePairs$Human != "", ]

  return(genePairs)
}

# ----------------------------
# Load mouse DC object with mouse genes
# ----------------------------
mouse_DCs_mGenes <- readRDS(mouse_dcs_rds)

# You can inspect clusters broadly here
DimPlot(mouse_DCs_mGenes, label = TRUE)

# ----------------------------
# Subset mouse pDC and pDC-like (AS DC) clusters with mouse genes
# ----------------------------
FeaturePlot(mouse_DCs_mGenes, genes =c("Siglech", "Cx3cr1", "Klf4", "Sirpa", "Bst2"))
mouse_pDCs_mGenes <- subset(mouse_DCs_mGenes, idents = c("0", "7"))
mouse_pDCs_mGenes <- RunUMAP(mouse_pDCs_mGenes, dims = 1:20)
DimPlot(mouse_pDCs_mGenes) 
#I made sure these are really pDC & pDC-like by doing re-integration of these cells with initial object.
#Because cycling pDC express most pDC genes, I assumed that I would not loose cycling pDCs by doing the analysis this way.

# ----------------------------
# Convert mouse genes to human orthologues
# ----------------------------
musGenes10 <- rownames(mouse_pDCs_mGenes)
mouse_to_human_orthologues <- convertMouseGeneList(musGenes10)

# Extract raw count matrix
counts_mouse <- GetAssayData(mouse_pDCs_mGenes, assay = "RNA", slot = "counts")
counts_mouse <- as.matrix(counts_mouse)
counts_mouse_df <- as.data.frame(counts_mouse)
counts_mouse_df$Mouse <- rownames(counts_mouse_df)

# Merge orthologue table with count matrix
mouse_counts_humanGeneSymbols <- merge(
  mouse_to_human_orthologues,
  counts_mouse_df,
  by = "Mouse",
  all = FALSE
)

# Reproducible handling of duplicate human gene symbols
mouse_counts_humanGeneSymbols$Human <- make.unique(mouse_counts_humanGeneSymbols$Human)

rownames(mouse_counts_humanGeneSymbols) <- mouse_counts_humanGeneSymbols$Human
mouse_counts_humanGeneSymbols$Human <- NULL
mouse_counts_humanGeneSymbols$Mouse <- NULL

mouse_counts_humanGeneSymbols <- mouse_counts_humanGeneSymbols[
  !duplicated(rownames(mouse_counts_humanGeneSymbols)), ,
  drop = FALSE
]

# ----------------------------
# Create mouse pDC object using human gene symbols
# ----------------------------
mouse_pDC_hGenes <- CreateSeuratObject(counts = as.matrix(mouse_counts_humanGeneSymbols))

mouse_pDC_hGenes[["percent.mt"]] <- PercentageFeatureSet(mouse_pDC_hGenes, pattern = "^MT-")
mouse_pDC_hGenes <- SCTransform(mouse_pDC_hGenes, vars.to.regress = "percent.mt", verbose = FALSE)
mouse_pDC_hGenes <- RunPCA(mouse_pDC_hGenes, verbose = FALSE)
mouse_pDC_hGenes <- FindNeighbors(mouse_pDC_hGenes, dims = 1:30)
mouse_pDC_hGenes <- FindClusters(mouse_pDC_hGenes, resolution = 0.1)
mouse_pDC_hGenes <- RunUMAP(mouse_pDC_hGenes, dims = 1:30)

# Transfer metadata from the original mouse pDC object (because now we have 2 different objects)
meta <- mouse_pDCs_mGenes[[]]
meta <- meta[colnames(mouse_pDC_hGenes), , drop = FALSE]
mouse_pDC_hGenes <- AddMetaData(mouse_pDC_hGenes, metadata = meta)

#Do a sanity check of expression & cell distribution
FeaturePlot(mouse_pDC_hGenes, features = c("VEGFB", "CX3CR1", "TCF4"), order = TRUE)
DimPlot(mouse_pDC_hGenes, label = TRUE)

if ("annotation_coarse" %in% colnames(mouse_pDC_hGenes@meta.data)) {
  mouse_pDC_hGenes@meta.data$annotation_coarse_species <- paste(
    "Mouse",
    mouse_pDC_hGenes@meta.data$annotation_coarse,
    sep = "_"
  )
}

mouse_pDC_hGenes@meta.data$species <- "Mouse"

if ("Tissue" %in% colnames(mouse_pDC_hGenes@meta.data)) {
  mouse_pDC_hGenes@meta.data$tissue_species <- paste(
    "Mouse",
    mouse_pDC_hGenes@meta.data$Tissue,
    sep = "_"
  )
}

saveRDS(mouse_pDC_hGenes, out_mouse_pdc_hGenes_rds) #these paths are defined at the top

# Mouse-only plot
p1 <- DimPlot(
  mouse_pDC_hGenes,
  group.by = "orig.ident",
  cols = "#4E9A06FF"
) & NoAxes() & NoLegend()

ggsave(
  filename = file.path(plot_dir, "p1_Mouse.png"),
  plot = p1,
  dpi = 700,
  width = 4,
  height = 4
)

# ----------------------------
# Load human pDC object with human genes
# ----------------------------
human_pDC_hGenes <- readRDS(human_pdc_rds)
Idents(human_pDC_hGenes) <- human_pDC_hGenes@meta.data$annotation_coarse #annotation_coarse is the cell type

# Downsample human object to match mouse cell number
n_mouse <- ncol(mouse_pDC_hGenes)
human_pDC_hGenes_downs <- subset(human_pDC_hGenes, downsample = n_mouse)

# ----------------------------
# Harmonize genes between human and mouse objects
# ----------------------------
common_genes <- intersect(rownames(human_pDC_hGenes_downs), rownames(mouse_pDC_hGenes))
stopifnot(length(common_genes) > 0)

human_pDC_hGenes_orthologous <- subset(human_pDC_hGenes_downs, features = common_genes)
mouse_pDC_hGenes <- subset(mouse_pDC_hGenes, features = common_genes)

# ----------------------------
# Plot input objects before CCA
# ----------------------------
p1 <- DimPlot(
  human_pDC_hGenes_orthologous,
  cols = c("#088BBEFF", "#FF3D7FFF", "#1BB6AFFF")
) & NoAxes() & NoLegend()

ggsave(
  filename = file.path(plot_dir, "p1_CCA.png"),
  plot = p1,
  dpi = 700,
  width = 4,
  height = 4
)

p2 <- DimPlot(
  mouse_pDC_hGenes,
  cols = c("#088BBEFF", "#088BBEFF", "#088BBEFF")
) & NoAxes() & NoLegend()

ggsave(
  filename = file.path(plot_dir, "p2_CCA.png"),
  plot = p2,
  dpi = 700,
  width = 4,
  height = 4
)

# ----------------------------
# Prepare metadata for cross-species CCA
# ----------------------------
human_pDC_hGenes_orthologous@meta.data$species <- "Human"
mouse_pDC_hGenes@meta.data$species <- "Mouse"

if ("annotation_coarse" %in% colnames(human_pDC_hGenes_orthologous@meta.data)) {
  human_pDC_hGenes_orthologous@meta.data$annotation_coarse_species <- paste(
    "Human",
    human_pDC_hGenes_orthologous@meta.data$annotation_coarse,
    sep = "_"
  )
}

if ("annotation_coarse" %in% colnames(mouse_pDC_hGenes@meta.data)) {
  mouse_pDC_hGenes@meta.data$annotation_coarse_species <- paste(
    "Mouse",
    mouse_pDC_hGenes@meta.data$annotation_coarse,
    sep = "_"
  )
}

if ("Tissue" %in% colnames(mouse_pDC_hGenes@meta.data)) {
  mouse_pDC_hGenes@meta.data$tissue_species <- paste(
    "Mouse",
    mouse_pDC_hGenes@meta.data$Tissue,
    sep = "_"
  )
}

if ("organ" %in% colnames(human_pDC_hGenes_orthologous@meta.data)) {
  human_pDC_hGenes_orthologous@meta.data$tissue_species <- paste(
    "Human",
    human_pDC_hGenes_orthologous@meta.data$organ,
    sep = "_"
  )
}

# ----------------------------
# Cross-species CCA integration
# ----------------------------
pdc_cca <- RunCCA(
  object1 = human_pDC_hGenes_orthologous,
  object2 = mouse_pDC_hGenes
)

pdc_cca <- FindNeighbors(pdc_cca, reduction = "cca", dims = 1:20)
pdc_cca <- FindClusters(pdc_cca, reduction = "cca")
pdc_cca <- RunUMAP(pdc_cca, reduction = "cca", dims = 1:20)

DimPlot(pdc_cca, group.by = "species")
DimPlot(pdc_cca, group.by = "annotation_coarse")
DimPlot(pdc_cca, split.by = "annotation_coarse_species", label = TRUE)

# ----------------------------
# Publication-quality CCA plots
# ----------------------------
p3 <- DimPlot(
  pdc_cca,
  split.by = "species",
  group.by = "species",
  cols = c("#FF3D7FFF", "#088BBEFF")
) & NoAxes() & NoLegend()

ggsave(
  filename = file.path(plot_dir, "p3_CCA_species.png"),
  plot = p3,
  dpi = 700,
  width = 6,
  height = 4
)

p4 <- FeaturePlot(
  pdc_cca,
  features = c("TLR7", "CX3CR1", "CD83", "MKI67"),
  order = TRUE,
  split.by = "species"
) & NoAxes()

ggsave(
  filename = file.path(plot_dir, "p4_CCA_markers.png"),
  plot = p4,
  dpi = 700,
  width = 6,
  height = 8
)

# ----------------------------
# Save integrated object
# ----------------------------
saveRDS(pdc_cca, out_pdc_cca_rds)
