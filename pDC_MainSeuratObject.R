#This is the code to reproduce the main pDC object used in Figure 1.
#Alina Ulezko Antonova

# Load packages ------------------------------------------------------------

library(Seurat)
library(plyr)
library(ggplot2)
library(harmony)

# Helper functions ---------------------------------------------------------

process_standard_sample <- function(path, project, sample_origin, age, sex, organ, patientID) {
  obj <- Seurat::Read10X(path)
  obj <- CreateSeuratObject(obj, project = project)

  obj$sample_origin <- sample_origin
  obj$age <- age
  obj$sex <- sex
  obj$organ <- organ
  obj$patientID <- patientID

  return(obj)
}

process_hto_sample <- function(path, use_kmeans = FALSE) {
  matrix <- Read10X(path)
  matrixGE <- matrix$`Gene Expression`
  matrixAC <- matrix$`Antibody Capture`
  joint.bcs <- intersect(colnames(matrixGE), colnames(matrixAC))
  matrixGE <- matrixGE[, joint.bcs]
  matrixAC <- as.matrix(matrixAC[, joint.bcs])
  rownames(matrixAC)

  hashtag <- CreateSeuratObject(counts = matrixGE)
  hashtag <- NormalizeData(hashtag)
  hashtag <- FindVariableFeatures(hashtag, selection.method = "mean.var.plot")
  hashtag <- ScaleData(hashtag, features = VariableFeatures(hashtag))
  hashtag[["HTO"]] <- CreateAssayObject(counts = matrixAC)
  hashtag <- NormalizeData(hashtag, assay = "HTO", normalization.method = "CLR")

  if (use_kmeans) {
    hashtag <- HTODemux(hashtag, assay = "HTO", positive.quantile = 0.99, kfunc = "kmeans")
  } else {
    hashtag <- HTODemux(hashtag, assay = "HTO", positive.quantile = 0.99)
  }

  table(hashtag$HTO_classification.global)

  Idents(hashtag) <- "HTO_maxID"
  RidgePlot(hashtag, assay = "HTO", features = rownames(hashtag[["HTO"]])[1:3], ncol = 3)

  Idents(hashtag) <- "HTO_classification.global"
  VlnPlot(hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

  hashtag.subset <- subset(hashtag, idents = "Negative", invert = TRUE)

  DefaultAssay(hashtag.subset) <- "HTO"
  hashtag.subset <- ScaleData(
    hashtag.subset,
    features = rownames(hashtag.subset),
    verbose = FALSE
  )
  hashtag.subset <- RunPCA(
    hashtag.subset,
    features = rownames(hashtag.subset),
    approx = FALSE
  )
  hashtag.subset <- RunTSNE(hashtag.subset, dims = 1:8, perplexity = 100)
  DimPlot(hashtag.subset)

  singlet <- subset(hashtag, idents = "Singlet")
  Idents(singlet) <- singlet@meta.data$HTO_classification
  singlet[["percent.mt"]] <- PercentageFeatureSet(singlet, pattern = "^MT-")

  VlnPlot(singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  VlnPlot(singlet, features = c("percent.mt"), ncol = 1, y.max = 20)
  VlnPlot(singlet, features = c("nFeature_RNA"), ncol = 1)

  return(singlet)
}

# Standard samples ---------------------------------------------------------

pbmc1 <- process_standard_sample(
  path = "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/samples/PBMC1/",
  project = "pbmc1",
  sample_origin = "pbmc1",
  age = "adult",
  sex = "Female",
  organ = "PBMC",
  patientID = "01"
)

pbmc3 <- process_standard_sample(
  path = "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/samples/PBMC3/",
  project = "pbmc3",
  sample_origin = "pbmc3",
  age = "adult",
  sex = "Male",
  organ = "PBMC",
  patientID = "02"
)

thymus1 <- process_standard_sample(
  path = "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/samples/THYMUS1/",
  project = "thymus1",
  sample_origin = "thymus1",
  age = "child",
  sex = "Female",
  organ = "Thymus",
  patientID = "03"
)

thymus2 <- process_standard_sample(
  path = "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/samples/THYMUS2/",
  project = "thymus2",
  sample_origin = "thymus2",
  age = "child",
  sex = "Male",
  organ = "Thymus",
  patientID = "04"
)

tonsil1 <- process_standard_sample(
  path = "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/samples/TONSIL1/",
  project = "tonsil1",
  sample_origin = "tonsil1",
  age = "child",
  sex = "Female",
  organ = "Tonsil",
  patientID = "05"
)

tonsil2 <- process_standard_sample(
  path = "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/samples/TONSIL2/",
  project = "tonsil2",
  sample_origin = "tonsil2",
  age = "child",
  sex = "Male",
  organ = "Tonsil",
  patientID = "06"
)

tonsil3 <- process_standard_sample(
  path = "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/samples/TONSIL3/",
  project = "tonsil3",
  sample_origin = "tonsil3",
  age = "child",
  sex = "Female",
  organ = "Tonsil",
  patientID = "07"
)

# Hashtagged samples -------------------------------------------------------

# sample 1910_11_M7WO (MGI4595) - works, but kfunc in HTODemux() has to be kmeans
sample191011.singlet <- process_hto_sample(
  path = "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/samples/Colonna_fbc_MGI4595_10X/MGI4595_MCAM-AA-1910_11_M7WO_pDC/outs/filtered_feature_bc_matrix/",
  use_kmeans = TRUE
)
sample191011.singlet$sample_origin <- sample191011.singlet@meta.data$HTO_classification
sample191011.singlet$sex <- "Male"
sample191011.singlet$age <- "7WO"
sample191011.singlet$patientID <- "191011"
sample191011.singlet$organ <- plyr::mapvalues(
  x = sample191011.singlet$sample_origin,
  from = c("LN191011-TotalSeqB", "PBMC191011-TotalSeqB", "Thymus191011-TotalSeqB"),
  to = c("LN", "PBMC", "Thymus")
)

# sample 1910_12_M2WO (MGI4658)
sample191012.singlet <- process_hto_sample(
  path = "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/samples/Colonna_fbc_MGI4658_10X/MGI4658_MCAM-AA-1910_12_M2WO_pDC/outs/filtered_feature_bc_matrix/"
)
sample191012.singlet$sample_origin <- sample191012.singlet@meta.data$HTO_classification
sample191012.singlet$sex <- "Male"
sample191012.singlet$age <- "2WO"
sample191012.singlet$patientID <- "191012"
sample191012.singlet$organ <- plyr::mapvalues(
  x = sample191012.singlet$sample_origin,
  from = c("LN191012-TotalSeqB", "PBMC191012-TotalSeqB", "Thymus191012-TotalSeqB"),
  to = c("LN", "PBMC", "Thymus")
)

# sample 1910_13_M6MO (MGI4669)
sample191013.singlet <- process_hto_sample(
  path = "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/samples/Colonna_fbc_MGI4669_10X/MGI4669_MCAM-AA-1910_13_6MO_pDC/outs/filtered_feature_bc_matrix/"
)
sample191013.singlet$sample_origin <- sample191013.singlet@meta.data$HTO_classification
sample191013.singlet$sex <- "Male"
sample191013.singlet$age <- "6MO"
sample191013.singlet$patientID <- "191013"
sample191013.singlet$organ <- plyr::mapvalues(
  x = sample191013.singlet$sample_origin,
  from = c("PBMC191013-TotalSeqB", "Thymus191013-TotalSeqB"),
  to = c("PBMC", "Thymus")
)

# sample 1910_14_F2WO (SR003091)
sample191014.singlet <- process_hto_sample(
  path = "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/samples/Colonna_fbc_SR003091_10X/SR003091_1910_14_F2WO_pDC/outs/while_waiting/"
)
sample191014.singlet$sample_origin <- sample191014.singlet@meta.data$HTO_classification
sample191014.singlet$sex <- "Female"
sample191014.singlet$age <- "2WO"
sample191014.singlet$patientID <- "191014"
sample191014.singlet$organ <- plyr::mapvalues(
  x = sample191014.singlet$sample_origin,
  from = c("LN191014-TotalSeqB", "PBMC191014-TotalSeqB", "Thymus191014-TotalSeqB"),
  to = c("LN", "PBMC", "Thymus")
)

# sample 1910_15_F4MO (SR003160) 
sample191015.singlet <- process_hto_sample(
  path = "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/samples/Colonna_fbc_SR003160_10X/SR003160_1910_15_F4MO_pDC/outs/filtered_feature_bc_matrix/"
)
sample191015.singlet$sample_origin <- sample191015.singlet@meta.data$HTO_classification
sample191015.singlet$sex <- "Female"
sample191015.singlet$age <- "4MO"
sample191015.singlet$patientID <- "191015"
sample191015.singlet$organ <- plyr::mapvalues(
  x = sample191015.singlet$sample_origin,
  from = c("PBMC191015-TotalSeqB", "Thymus191015-TotalSeqB"),
  to = c("PBMC", "Thymus")
)

# sample 1910_16_F8DO (SR003223)
sample191016.singlet <- process_hto_sample(
  path = "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/samples/Colonna_fbc_SR003223_10X/SR003223_1910_16_F8DO_pDC/outs/filtered_feature_bc_matrix/"
)
sample191016.singlet$sample_origin <- sample191016.singlet@meta.data$HTO_classification
sample191016.singlet$sex <- "Female"
sample191016.singlet$age <- "8DO"
sample191016.singlet$patientID <- "191016"
sample191016.singlet$organ <- plyr::mapvalues(
  x = sample191016.singlet$sample_origin,
  from = c("LN191016-TotalSeqB", "PBMC191016-TotalSeqB", "Thymus191016-TotalSeqB"),
  to = c("LN", "PBMC", "Thymus")
)

# Merge all samples --------------------------------------------------------

all_pDC <- merge(
  pbmc1,
  y = c(
    pbmc3, thymus1, thymus2, tonsil1, tonsil2, tonsil3,
    sample191011.singlet, sample191012.singlet, sample191013.singlet,
    sample191014.singlet, sample191015.singlet, sample191016.singlet
  )
)

# Subset data and QC -------------------------------------------------------

all_pDC[["percent.mt"]] <- PercentageFeatureSet(all_pDC, pattern = "^MT-")
p1 <- VlnPlot(all_pDC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/plots/qc.png", p1, width = 25, height = 10, dpi = 700)
VlnPlot(all_pDC, features = c("percent.mt"), ncol = 1, y.max = 20)
VlnPlot(all_pDC, features = c("nFeature_RNA"), ncol = 1)
saveRDS(all_pDC, "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/objects/all_pDC_030524.rds")
# all_pDC <- readRDS("/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/objects/all_pDC_022224.rds")

# SCTransform, Harmony, clustering -----------------------------------------

all_pDC1 <- subset(all_pDC, subset = nFeature_RNA > 1500 & nFeature_RNA < 6000 & percent.mt < 7.5)
all_pDC1 <- SCTransform(all_pDC1, vars.to.regress = "percent.mt")
all_pDC1 <- RunPCA(all_pDC1)
all_pDC1 <- RunHarmony(all_pDC1, assay = "SCT", group.by.vars = "patientID")
all_pDC1 <- FindNeighbors(all_pDC1, dims = 1:20, reduction = "harmony")
all_pDC1 <- FindClusters(all_pDC1, resolution = 0.3, reduction = "harmony")
all_pDC1 <- RunUMAP(all_pDC1, dims = 1:20, reduction = "harmony")
DimPlot(all_pDC1, split.by = "organ")
DimPlot(all_pDC1, split.by = "sex")
DimPlot(all_pDC1, split.by = "patientID")
DimPlot(all_pDC1, label = TRUE)
saveRDS(all_pDC1, "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/objects/all_pDC_SCT_030524.rds")

# Remove contaminating clusters -----------------------

# Cluster 13 (n = 245) was excluded because it was detected only in sample 1910-11,
# suggesting a sample-specific population rather than a reproducible pDC state.
# Cluster 14 (n = 66) was excluded because it expressed T-cell-associated genes suggesting
# that this is a pDC-T cell doublet.
# These two clusters were therefore removed prior to downstream pDC-focused analyses.

all_pDC2 <- subset(all_pDC1, idents = c("13", "14"), invert = TRUE)
all_pDC2 <- FindNeighbors(all_pDC2, dims = 1:20, reduction = "harmony")
all_pDC2 <- FindClusters(all_pDC2, resolution = 0.3, reduction = "harmony")
all_pDC2 <- RunUMAP(all_pDC2, dims = 1:20, reduction = "harmony")
DimPlot(all_pDC2, label = TRUE)
DimPlot(all_pDC2, split.by = "patientID")
DimPlot(all_pDC2, group.by = "organ")
DimPlot(all_pDC2, split.by = "sample_origin", ncol = 4)

# Add annotations ----------------------------------------------------------
# Initially called pDC "Bona fide pDC" and AS DC "pDC-like"
all_pDC2$annotation_coarse <- plyr::mapvalues(
  x = all_pDC2$seurat_clusters,
  from = c("12", "10", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "11"),
  to = c(
    "Cycling pDC", "pDC-like",
    "Bona fide pDC", "Bona fide pDC", "Bona fide pDC", "Bona fide pDC", "Bona fide pDC",
    "Bona fide pDC", "Bona fide pDC", "Bona fide pDC", "Bona fide pDC", "Bona fide pDC", "Bona fide pDC"
  )
)

all_pDC2$age_group <- plyr::mapvalues(
  x = all_pDC2$sample_origin,
  from = c(
    "LN191011-TotalSeqB", "LN191012-TotalSeqB", "LN191014-TotalSeqB", "LN191016-TotalSeqB", "pbmc1", "PBMC191011-TotalSeqB",
    "PBMC191012-TotalSeqB", "PBMC191013-TotalSeqB", "PBMC191014-TotalSeqB", "PBMC191015-TotalSeqB", "PBMC191016-TotalSeqB", "pbmc3",
    "thymus1", "Thymus191011-TotalSeqB", "Thymus191012-TotalSeqB", "Thymus191013-TotalSeqB", "Thymus191014-TotalSeqB", "Thymus191015-TotalSeqB",
    "Thymus191016-TotalSeqB", "thymus2", "tonsil1", "tonsil2", "tonsil3"
  ),
  to = c(
    "Child", "Child", "Child", "Child", "Adult", "Child",
    "Child", "Child", "Child", "Child", "Child", "Adult",
    "Child", "Child", "Child", "Child", "Child", "Child",
    "Child", "Child", "Child", "Child", "Child"
  )
)

# Save final object --------------------------------------------------------

saveRDS(all_pDC2, "/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/objects/all_pDC_SCT_Clean.rds")
all_pDC2 <- readRDS("/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/objects/all_pDC_SCT_Clean.rds")
