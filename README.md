# Spatial Transcriptomics Analysis of Mouse Brain Tissue

This repository contains the analysis pipeline and results for a spatial transcriptomics project focused on characterizing cell type distributions and gene expression patterns in mouse brain tissue using 10x Genomics Visium technology.

## Project Overview

This project involves comprehensive analysis of two mouse brain tissue sections using Visium spatial transcriptomics. The analysis includes quality control, normalization, dimensionality reduction, clustering, spatial variable gene identification, cell type deconvolution, and cell-cell communication analysis.

## Data Summary

- Technology: 10x Genomics Visium spatial transcriptomics
- Samples: 2 mouse brain tissue sections
- Feature dimensions:
  - Section 1: 32,285 genes × 3,355 spots
  - Section 2: 32,285 genes × 3,289 spots
- Spot size: 55 μm
- Spot distance: 100 μm center-to-center
- Total spots per capture area: 4,992

## Analysis Workflow

### 1. Data Structure and Organization

- Gene expression matrix: Genes × Spots (32,285 × ~3,300 per section)
- Spatial coordinates: Stored in "spatial/tissue_positions_list.csv"
- Images: High-resolution (~2000 × 2000 pixels) and low-resolution (600 × 600 pixels) images stored in the "spatial" folder

### 2. Quality Control and Preprocessing

Filtering thresholds for Visium data:
- nFeature_Spatial: 1,750-9,000 genes per spot
- nCount_Spatial: 8,000-50,000 UMIs per spot

These thresholds are higher than typical single-cell RNA-seq thresholds (nFeature_RNA: 750-2,500 genes, nCount_RNA: 1,500-6,000 UMIs) due to:
- Multiple cells captured per Visium spot (typically 1-10 cells)
- More efficient RNA recovery with direct tissue permeabilization
- Higher sequencing depth per spot

### 3. Normalization and Feature Selection

- Used SCTransform for normalization, which replaces the traditional Seurat workflow of NormalizeData, ScaleData, and FindVariableFeatures
- SCTransform better handles the technical characteristics of spatial transcriptomics data compared to standard normalization methods

### 4. Dimensionality Reduction and Clustering

- Performed PCA and identified significant principal components using elbow plots
- Generated UMAP embeddings for visualization
- Identified 16 distinct clusters across the merged dataset

### 5. Spatial Variable Gene Identification

Used Moran's I spatial autocorrelation statistic to identify spatially variable genes:

Top spatially variable genes in Section 1:
- Nrgn (p-val: 1.49E-134, avg_log2FC: 2.03, cluster: 1)
- Pcp2 (p-val: 8.97E-66, avg_log2FC: 3.17, cluster: 11)
- Pvalb (p-val: 1.84E-61, avg_log2FC: 2.73, cluster: 11)

### 6. Multi-section Integration and Batch Effect Analysis

- Merged both tissue sections and analyzed cluster compositions
- Found that clusters were evenly distributed between sections (approximately 50% from each section)
- No clusters were specific to a single section
- Tested batch correction methods but found minimal batch effects, making batch correction unnecessary

### 7. Cell Type Identification and Annotation

Identified major brain cell populations including:
- Macrophages (markers: Dock2, Apoe)
- Astrocytes (markers: Gfap, Aldh1l1)
- Oligodendrocytes

Used both reference-based annotation from scRNA-seq data integration and marker gene expression analysis.

### 8. Cell Type Deconvolution

- Applied SCDC (Single Cell Deconvolution) to estimate cell type proportions in each spot
- Compared deconvolution results with cell type prediction through scRNA-seq data integration
- Results showed consistent patterns for major cell types including astrocytes and oligodendrocytes

### 9. Cell-Cell Communication Analysis

- Used CellChat to analyze cell-cell interactions
- Focused on the glutamate signaling pathway
- Generated interaction networks between different cell populations, considering their spatial relationships

## Technical Implementation

All analyses were performed using R with the following key packages:
- Seurat (v4.x) for general spatial transcriptomics analysis
- SeuratData for reference datasets
- SCDC for cell type deconvolution
- CellChat for cell-cell communication analysis

## Key Findings

1. Consistent cellular composition across both tissue sections with no significant batch effects
2. Identification of 16 distinct spatial domains in the mouse brain tissue
3. Clear spatial organization of major brain cell types
4. Agreement between deconvolution and reference-based annotation approaches
5. Identification of key communication pathways between cell types

## Future Directions

Alternative approaches to consider for future analysis:
- Cell2Location for improved spatial deconvolution
- SpatialDE as an alternative for spatial variable gene detection
- Squidpy for spatial domain analysis
- stLearn for trajectory analysis in spatial context
- Python-based spatial analysis pipelines

## Contact

For questions or collaborations, please contact Zain at zazi00001@stud.uni-saarland.de.
