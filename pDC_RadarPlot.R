#### Radar Plot (Figure 1D) ####
# Alina Ulezko Antonova

library(triwise)
library(limma)
library(Biobase)
library(Seurat)
library(ggplot2)

#Set working directory
setwd("/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/")

#Load object
all_pDC2 <- readRDS("/objects/all_pDC2.rds")
#Set identity to cell type (pDC/AS_DC/Cycling_pDC; in my case, this is in the metadata slot "annotation_coarse")
Idents(all_pDC2) <- all_pDC2@meta.data$annotation_coarse
#Calculate average expression of every gene
pdc_avg_coarse <- AverageExpression(all_pDC2, assay ="SCT", return.seurat = T)
pdc_avg_coarse_counts <- pdc_avg_coarse@assays$SCT@counts
pdc_avg_coarse_counts_df <- as.data.frame(pdc_avg_coarse_counts)
#Manually annotate which genes you would like to be displayed in the plot
genes_pDC_coarse <- c("AXL", "SIGLEC6", 
                      "FGL2", "COTL1", "HLA-DRB1", "KLF4", 
                      "CD63", "SPI1", "IL1B", "IFI30", "S100A10",
                      "LST1", "ANXA2", "VIM", "VSIR", "PPP1R14A", "LTK",
                      "NFKBIA", "AREG", "LILRA4", "KLF6", 
                      "IRF4", "PLXNA4",
                      "GZMB", "NPC2", "PTCRA", "IRF7", "IGKC", "PPP1R14B",
                      "IL3RA", "JCHAIN", "CLEC4C", "TCF4", "IGLC2", "MX1", "IFI16", "TRDC",
                      "IGLL1", "MKI67", "TOP2A", "TUBA1B", "HMGB1", "COTL1A", "SEMA4A")
barycoords_pdc_avg_coarse <- transformBarycentric(as.matrix(pdc_avg_coarse_counts_df))
k <- plotDotplot(barycoords_pdc_avg_coarse, Goi=genes_pDC_coarse)
#Save as you prefer (png/svg/pdf)
ggsave("/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/plots/dotplot_highlighted.png",k, dpi=700, width=6, height=6)
ggsave("/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/plots/dotplot_highlighted.svg",k, dpi=700, width=8, height=8)
ggsave("/storage1/fs1/mcolonna/Active/Alina/Own_Analysis/pDC_Project/plots/dotplot_highlighted.pdf",k, dpi=700, width=8, height=8)
#Something useful you can do is have an interactive dotplot pulled in the viewer and you can search your gene of interest to see where it is located at. 
#This is useful to manually annotate the genes, makes a cleaner plot.
dotplot = interactiveDotplot(as.matrix(pdc_avg_coarse_counts_df))
dotplot


dotplot <- interactiveDotplot(as.matrix(pdc_avg_coarse_counts_df))
dotplot
