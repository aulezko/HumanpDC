# HumanpDC

This repository contains the code associated with the publication:
**<u>Age- and tissue-dependent diversity of human plasmacytoid dendritic cells uncovers a cycling subset dominant in early life and cancer<u>**, by *Alina Ulezko Antonova et al.* (Immunity, 2026).

## Overview
This codebase provides the analytical pipeline used to process, integrate, and visualize single-cell transcriptomic data for human and mouse plasmacytoid dendritic cells (pDCs). 

I have uploaded these codes hoping they will be helpful for users to reproduce the data presented here and use them for their own research and learning. I personally benefited significantly from public code when learning to analyze single-cell data, and I hope this repository proves useful to you. 

*Note: Code not available in this repository follows standard pipelines as outlined in the methods section of the publication.*

## Repository Contents

* `pDC_Fig1B_HumanSeurat.R`: Script to reproduce the main Seurat object, as presented in Figure 1B.
* `pDC_Fig1D_RadarPlot.R`: Script used to generate the radar plot shown in Figure 1D.
* `pDC_Fig3_HumanMouseCCA.R`: Script for cross-species integration. This workflow converts mouse pDC gene symbols to human orthologues, rebuilds the mouse object in the human gene space, and performs Canonical Correlation Analysis (CCA) to integrate mouse and human datasets.
*  `pDC_Fig7_BPDCN.R`: Script to reproduce the main plots in the BPDCN analysis from *Griffin et al.*, as shown in Figure 7.
*  `Link to Zenodo for Seurat object`: See bottom of this page.

## Reproducibility
To ensure results are consistent, the mapping of mouse-to-human orthologues uses the Ensembl archive (December 2021). 

## Main Seurat Object
For ease of analysis for users, I have uploaded the main seurat object (from Figure 1B) to Zenodo. 
This can be accessed in this link:
