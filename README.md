# HumanpDC

This repository contains the code associated with the publication:
**<u>Age- and tissue-dependent diversity of human plasmacytoid dendritic cells uncovers a cycling subset dominant in early life and cancer<u>**, by *Alina Ulezko Antonova et al.* (Immunity, 2026).

## Overview
This codebase provides the analytical pipeline used to process, integrate, and visualize single-cell transcriptomic data for human and mouse plasmacytoid dendritic cells (pDCs).

## Repository Contents

* `pDC_Fig1B_HumanSeurat.R`: Script to reproduce the main Seurat object, as presented in Figure 1B.
* `pDC_Fig1D_RadarPlot.R`: Script used to generate the radar plot shown in Figure 1D.
* `pDC_Fig3_HumanMouseCCA.R`: Script for cross-species integration. This workflow converts mouse pDC gene symbols to human orthologues, rebuilds the mouse object in the human gene space, and performs Canonical Correlation Analysis (CCA) to integrate mouse and human datasets.
*  `pDC_Fig7_BPDCN.R`: Script to reproduce the main plots in the BPDCN analysis from *Griffin et al.*, as shown in Figure 7.

## Reproducibility
To ensure results are consistent, the mapping of mouse-to-human orthologues uses the Ensembl archive (December 2021). 
