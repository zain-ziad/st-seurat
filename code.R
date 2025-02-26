# setting seed for reproducibility
set.seed(11)
options(future.globals.maxSize = 8000 * 1024^2)

# loading necessary libraries
library("Seurat")
library("CellChat")
library("SCDC")
library("patchwork")

# make plots dir
if (!dir.exists("plots")) {
  dir.create("plots")
  dir.create("plots/preprocess")
  dir.create("plots/preprocess/sec1")
  dir.create("plots/preprocess/sec2")
  dir.create("plots/dim-reduction")
  dir.create("plots/dim-reduction/sec1")
  dir.create("plots/dim-reduction/sec2")
  dir.create("plots/deg-analysis")
  dir.create("plots/deg-analysis/sec1")
  dir.create("plots/deg-analysis/sec2")
  dir.create("plots/merge-nobc")
  dir.create("plots/merge-wbc")
  dir.create("plots/integration")
  dir.create("plots/single-cell")
  dir.create("plots/celltype-pred")
  dir.create("plots/celltype-pred/sec1")
  dir.create("plots/celltype-pred/sec2")
  dir.create("plots/celltype-pred/merge")
  dir.create("plots/celltype-pred/deg-analysis")
  dir.create("plots/celltype-pred/manual")
  dir.create("plots/deconv")
  dir.create("plots/cc-com")
}

# loading 10x visium data
sec1 <- Load10X_Spatial("project_3_dataset\\Section_1", 
                        "V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.h5")
sec1$orig.ident <- "Section_1"

sec2 <- Load10X_Spatial("project_3_dataset\\Section_2", 
                        "V1_Mouse_Brain_Sagittal_Posterior_Section_2_filtered_feature_bc_matrix.h5")
sec2$orig.ident <- "Section_2"

# overview of the datasets
p1 <- SpatialDimPlot(sec1) + NoLegend() +ggtitle("Section 1")
p2 <- SpatialDimPlot(sec2) + NoLegend() +ggtitle("Section 2")
ggsave("plots/overview.pdf", p1+p2)

# for gene expression data
dim(sec1@assays$Spatial$counts)
dim(sec2@assays$Spatial$counts)

# for spot coords
print(dim(GetTissueCoordinates(sec1)))
print(dim(GetTissueCoordinates(sec2)))

# images
dim(sec1@images$slice1@image)
dim(sec2@images$slice1@image)

# visualization of a feature
p1 <- SpatialFeaturePlot(sec1, features = "Snap25", slot = "counts") +ggtitle("Section 1 Snap25")
p2 <- SpatialFeaturePlot(sec2, features = "Snap25", slot = "counts") +ggtitle("Section 2 Snap25")
ggsave("plots/visualization-of-a-feature.pdf", p1+p2)

# preprocessing
preprocess <- function(dir, data, nFeature_Spatial_up, nFeature_Spatial_down, nCount_Spatial_up, nCount_Spatial_down) {
    p1 <- SpatialFeaturePlot(data, features = "nCount_Spatial") +ggtitle("Count distribution across tissue")
    p2 <- SpatialFeaturePlot(data, features = "nFeature_Spatial") +ggtitle("Gene counts per spot")
    ggsave(paste0(dir, "before-spatial.pdf"), p1+p2)

    p_violin <- VlnPlot(data, features = c("nCount_Spatial", "nFeature_Spatial"))
    ggsave(paste0(dir, "before-vln.pdf"), p_violin)

    # create histograms
    p1_hist <- ggplot(data.frame(counts = data$nCount_Spatial), aes(x = counts)) +
    geom_histogram(bins = 100, fill = "skyblue", color = "black") +
    labs(
        x = "Counts per Spot",
        y = "Frequency") +
    theme_minimal()

    p2_hist <- ggplot(data.frame(features = data$nFeature_Spatial), aes(x = features)) +
    geom_histogram(bins = 100, fill = "lightgreen", color = "black") +
    labs(
        x = "Genes per Spot",
        y = "Frequency") +
    theme_minimal()

    ggsave(paste0(dir, "before-hist.pdf"), p1_hist+p2_hist)

    # filter on both counts and features
    data <- subset(data, 
                nFeature_Spatial > nFeature_Spatial_down & nFeature_Spatial < nFeature_Spatial_up &
                nCount_Spatial > nCount_Spatial_down & nCount_Spatial < nCount_Spatial_up)

    p1 <- SpatialFeaturePlot(data, features = "nCount_Spatial") +ggtitle("Count distribution across tissue")
    p2 <- SpatialFeaturePlot(data, features = "nFeature_Spatial") +ggtitle("Gene counts per spot")
    ggsave(paste0(dir, "after-spatial.pdf"), p1+p2)

    p_violin <- VlnPlot(data, features = c("nCount_Spatial", "nFeature_Spatial")) +ggtitle("Distribution before filtering")
    ggsave(paste0(dir, "after-vln.pdf"), p_violin)

    return(data)
}

# data preprocessing sec1
dir <- "plots/preprocess/sec1/"
summary(sec1$nFeature_Spatial)
summary(sec1$nCount_Spatial)
sec1 <- preprocess(dir, sec1, 9000, 1750, 50000, 8000)

# data preprocessing sec1
dir <- "plots/preprocess/sec2/"
summary(sec2$nFeature_Spatial)
summary(sec2$nCount_Spatial)
sec2 <- preprocess(dir, sec2, 8500, 1650, 45000, 7000)

# sctransform
sec1 <- SCTransform(sec1, assay = "Spatial", verbose = FALSE)
sec2 <- SCTransform(sec2, assay = "Spatial", verbose = FALSE)

# dimensionality reduction
dim_redution <- function(dir, data, pcs, cl_res) {
    data <- RunPCA(data, assay = "SCT", verbose = FALSE)
    elbow <- ElbowPlot(data, ndims = 50, reduction = "pca")
    ggsave(paste0(dir, "elbow-plot.pdf"), elbow)
    data <- FindNeighbors(data, reduction = "pca", dims = 1:pcs)
    data <- FindClusters(data, res = cl_res, verbose = FALSE)
    data <- RunUMAP(data, reduction = "pca", dims = 1:pcs)
    p1 <- DimPlot(data, reduction = "umap", label = TRUE)
    p2 <- SpatialDimPlot(data, label = TRUE, label.size = 3)
    ggsave(paste0(dir, "dim-plot.pdf"), p1, width = 10)
    ggsave(paste0(dir, "spatial-dim-plot.pdf"), p2)
    return(data)
}

dir <- "plots/dim-reduction/sec1/"
sec1 <- dim_redution(dir, sec1, 30, 0.8)

dir <- "plots/dim-reduction/sec2/"
sec2 <- dim_redution(dir, sec2, 30, 0.8)

# DEG (Differentially Expressed Genes) analysis based on clustering
deg_analysis <- function(dir, data, assay) {
    data <- PrepSCTFindMarkers(data)
    markers <- FindAllMarkers(data)
    data <- FindSpatiallyVariableFeatures(data, assay = assay, features = VariableFeatures(data)[1:1000],
    selection.method = "moransi")
    spatial_data <- data[[assay]]@meta.features
    top_spatial_genes <- rownames(spatial_data[
    spatial_data$MoransI_p.value < 0.05, ])[
    order(spatial_data[
        spatial_data$MoransI_p.value < 0.05, 
        "moransi.spatially.variable.rank"
    ])[1:25]]
    spatial_data <- arrange(spatial_data, moransi.spatially.variable.rank)
    write.csv(spatial_data, paste0(dir, "spatial_variable_features.csv"))
    write.csv(markers, paste0(dir, "markers.csv"))

    print("Top 3 Markers:")
    print(rownames(markers)[1:3])   
    print("Top 3 spatially variable genes:")
    print(top_spatial_genes[1:3])   
    common_genes <- intersect(rownames(markers), top_spatial_genes)
    print("Genes that are both markers and spatially variable:")
    print(common_genes)
    markers_subset <- markers[markers$gene %in% top_spatial_genes[1:3], ] %>%
        group_by(gene) %>%
        slice_max(avg_log2FC, n=1) %>%
        ungroup()
    print("Marker statistics for these genes:")
    print(markers_subset)
    write.csv(markers_subset, paste0(dir, "svg-markers.csv"))

    p1 <- SpatialFeaturePlot(sec1, features = top_spatial_genes[1:3], ncol = 3, alpha = c(0.1, 1))
    ggsave(paste0(dir, "svg.pdf"), p1)

    return(data)
}

dir <- "plots/deg-analysis/sec1/"
sec1 <- deg_analysis(dir, sec1, "SCT")

# merging without Batch-correction
merging_samples <- function(dir, data1, data2, label_celltype) {
    data <- merge(data1, y = data2, add.cell.ids = c("sec1", "sec2"), merge.data = TRUE)
    data <- SCTransform(data, assay = "Spatial", verbose = FALSE)
    DefaultAssay(data) <- "SCT"
    VariableFeatures(data) <- c(VariableFeatures(data1), VariableFeatures(data2))
    data <- dim_redution(dir, data, 20, 0.5)
    p1 <- DimPlot(data, reduction = "umap", group.by = "orig.ident")
    ggsave(paste0(dir, "samples-umap.pdf"), p1, width = 10)
    percent_data <- prop.table(table(data$seurat_clusters, data$orig.ident), margin = 1) * 100
    percent_data <- round(percent_data, 2)
    write.csv(percent_data, paste0(dir, "cluster-composition.csv"))
    if (label_celltype){
        p1 <- SpatialDimPlot(data, 
               group.by = "predicted_celltype",
               pt.size = 1.5,
               label = FALSE)
        ggsave(paste0(dir, "spatial-celltypes.pdf"), p1, width = 10)    
    }
    return(data)
}

dir <- "plots/merge-nobc/"
sec_merged <- merging_samples(dir, sec1, sec2, FALSE)

# with batch correction
integrate_datasets <- function(dir, data1, data2) {
    datalist <- list(data1, data2)
    features <- SelectIntegrationFeatures(object.list = datalist, nfeatures = 2000)
    list <- PrepSCTIntegration(object.list = datalist, anchor.features = features)
    intanchors <- FindIntegrationAnchors(
        object.list = list,
        normalization.method = "SCT",
        anchor.features = features,
        dims = 1:30
    )
    integrated <- IntegrateData(anchorset = intanchors, dims = 1:30)
    DefaultAssay(object = integrated) <- "integrated"
    integrated <- SCTransform(integrated, assay = "Spatial")
    integrated <- RunPCA(integrated, verbose = FALSE)
    integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)

    p1 <- DimPlot(integrated, group.by = "orig.ident")
    p2 <- DimPlot(integrated, group.by = "ident")
    ggsave(paste0(dir, "samples-umap.pdf"), p1, width = 10)
    ggsave(paste0(dir, "clusters-umap.pdf"), p2, width = 16)

    return(integrated)
}

dir <- "plots/integration/"
sec_integrated <- integrate_datasets(dir, sec1, sec2)

# Automatic Annotation using Data Integration
sc_data <- readRDS("project_3_dataset\\allen_cortex.rds")
sc_data <- SCTransform(sc_data, ncells = 3000)
sc_data <- RunPCA(sc_data, verbose = FALSE)
sc_data <- RunUMAP(sc_data, dims = 1:30)

dir <- "plots/single-cell/"
p1 <- DimPlot(sc_data, group.by = "class", label = TRUE)
p2 <- DimPlot(sc_data, group.by = "subclass", label = TRUE)
ggsave(paste0(dir, "scrna-class.pdf"), p1, width = 10)
ggsave(paste0(dir, "scrna-subclass.pdf"), p2, width = 10)

add_annotation <- function(dir, scdata, spdata){
    anchors <- FindTransferAnchors(reference = scdata, query = spdata, dims = 1:30, normalization.method = "SCT")
    predictions <- TransferData(anchorset = anchors, refdata = scdata$subclass, prediction.assay = TRUE, weight.reduction = spdata[["pca"]], dims = 1:30)
    spdata[["predictions"]] <- predictions
    DefaultAssay(spdata) <- "predictions"
    pred_matrix <- GetAssayData(spdata, assay = "predictions")
    non_zero_rows <- rowSums(pred_matrix) > 0
    cell_types <- rownames(pred_matrix)[non_zero_rows]
    p1 <- SpatialFeaturePlot(spdata, 
                    features = cell_types,
                    ncol = 6)
    ggsave(paste0(dir, "anno-scores.pdf"), p1, width = 12)
    pred_matrix <- GetAssayData(spdata, assay = "predictions")
    max_predictions <- apply(pred_matrix, 2, function(x) rownames(pred_matrix)[which.max(x)])
    spdata$predicted_celltype <- max_predictions

    p1 <- DimPlot(spdata, reduction = "umap", group.by = "predicted_celltype")
    ggsave(paste0(dir, "anno-umap.pdf"), p1, width = 12)

    return(spdata)
}

#sec1 <- LoadSeuratRds("sec1.rds")
#sec2 <- LoadSeuratRds("sec2.rds")

dir <- "plots/celltype-pred/sec1/"
sec1 <- add_annotation(dir, sc_data, sec1)
dir <- "plots/celltype-pred/sec2/"
sec2 <- add_annotation(dir, sc_data, sec2)

dir <- "plots/celltype-pred/merge/"
DefaultAssay(sec1) <- "SCT"
DefaultAssay(sec2) <- "SCT"
sec <- merging_samples(dir, sec1, sec2, TRUE)
p1 <- DimPlot(sec, reduction = "umap", group.by = "predicted_celltype")
ggsave(paste0(dir, "anno-umap.pdf"), p1, width = 12)

# Manual Annotation
slot(object = sec@assays$SCT@SCTModel.list[[2]], name="umi.assay")<-"Spatial"
DefaultAssay(sec) <- "SCT"

macrophage <- c("Dock2", "Apoe")
astro <- c("Gfap", "Aldh1l1")
oligo <- c("Olig2", "Mbp")

plot_cell_type_markers <- function(dir, data, gene) {
    p1 <- FeaturePlot(data, 
                    features = gene,
                    reduction = "umap") +
        ggtitle(paste0(gene, " UMAP"))

    p2 <- SpatialFeaturePlot(data, 
                        features = gene) +
        ggtitle(paste0(gene, " Spatial"))

    combined <- p1 + p2
    ggsave(filename = file.path(dir, paste0(gene, ".pdf")),
           plot = combined,
           width = 12, 
           height = 16)
}

# create plots for each cell type
dir <- "plots/celltype-pred/manual/"

for(g in macrophage) {
    plot_cell_type_markers(dir, sec, g)
}
for(g in astro) {
    plot_cell_type_markers(dir, sec, g)
}
for(g in oligo) {
    plot_cell_type_markers(dir, sec, g)
}

# deconv
# prepare the reference data
sc_data <- readRDS("project_3_dataset\\allen_cortex.rds")
cells_to_keep <- unlist(lapply(unique(sc_data$subclass), function(type) {
   cells <- names(sc_data$subclass)[sc_data$subclass == type]
   sample(cells, min(250, length(cells)))
}))

sc_data <- subset(sc_data, cells = cells_to_keep)

shared_celltypes <- intersect(unique(sc_data$subclass), unique(sec$predicted_celltype))

markers_list <- list()
for(celltype in shared_celltypes) {
   markers <- FindMarkers(sc_data, 
                         ident.1 = celltype,
                         group.by = "subclass",
                         min.pct = 0.25)
   
   markers_list[[celltype]] <- head(rownames(markers[order(markers$p_val),]), 20)
}

genes <- unique(unlist(markers_list))
sp_genes <- rownames(GetAssayData(sec, slot = "counts"))
genes_overlap <- intersect(genes, sp_genes)


sc_counts <- GetAssayData(sc_data, slot = "counts")
sc_counts <- sc_counts[genes_overlap,]

sc_pdata <- data.frame(
    celltype = sc_data$subclass,
    sample = colnames(sc_counts),
    row.names = colnames(sc_counts)
)

sc_eset <- getESET(exprs = sc_counts,
                   fdata = rownames(sc_counts),
                   pdata = sc_pdata)

sp_counts <- GetAssayData(sec, slot = "counts")
sp_counts <- sp_counts[genes_overlap,]

sp_pdata <- data.frame(
    sample = colnames(sp_counts),
    row.names = colnames(sp_counts)
)

sp_eset <- getESET(exprs = sp_counts,
                   fdata = rownames(sp_counts),
                   pdata = sp_pdata)

sc_celltypes <- unique(sc_data$subclass)
sp_celltypes <- unique(sec$predicted_celltype)
overlapping_celltypes <- intersect(sc_celltypes, sp_celltypes)

results <- SCDC_prop(bulk.eset = sp_eset, 
                    sc.eset = sc_eset,
                    ct.varname = "celltype",
                    sample = "sample",
                    ct.sub = overlapping_celltypes,
                    weight.basis = FALSE)

proportions <- results$prop.est.mvw
proportions_matrix <- t(proportions)
colnames(proportions_matrix) <- rownames(proportions)
rownames(proportions_matrix) <- colnames(proportions)

deconv_assay <- CreateAssayObject(data = proportions_matrix)
DefaultAssay(sec) <- "Spatial"
sec[["SCDC"]] <- deconv_assay

DefaultAssay(sec) <- "SCDC"
p1 <- SpatialFeaturePlot(sec, features = "Astro")
p2 <- SpatialFeaturePlot(sec, features = "Oligo")

DefaultAssay(sec) <- "predictions"
p3 <- SpatialFeaturePlot(sec, features = "Astro")
p4 <- SpatialFeaturePlot(sec, features = "Oligo")

dir <- "plots/deconv/"
ggsave(paste0(dir, "astro.pdf"), p1, width = 10)
ggsave(paste0(dir, "oligo.pdf"), p2, width = 10)
ggsave(paste0(dir, "astro_pred.pdf"), p3, width = 10)
ggsave(paste0(dir, "oligo_pred.pdf"), p4, width = 10)

# cell-cell communication
#saveRDS(sec, "sec.rds")
#sec <- LoadSeuratRds("sec.rds")

data.input <- GetAssayData(sec, slot = "counts", assay = "SCT")
data.input <- as.matrix(data.input)
deconv_data <- GetAssayData(sec, slot = "data", assay = "SCDC")
cell_labels <- apply(deconv_data, 2, function(x) {
    rownames(deconv_data)[which.max(x)]
})

meta <- data.frame(group = cell_labels, 
                  row.names = names(cell_labels))

# define image coords for both images
spatial_coords <- NULL
for(img in Images(sec)) {
    coords <- GetTissueCoordinates(sec, image = img)
    coords$image <- img
    spatial_coords <- rbind(spatial_coords, coords)
}
spatial_coords <- as.matrix(spatial_coords[, c("x", "y")])

scalefactors <- jsonlite::fromJSON(txt = file.path("project_3_dataset\\Section_1\\spatial", 'scalefactors_json.json'))
spot.size <- 65 
conversion.factor <- spot.size/scalefactors$spot_diameter_fullres
spatial.factors <- data.frame(ratio = conversion.factor, tol = spot.size/2)

cellchat <- createCellChat(object = data.input, 
                          meta = meta, 
                          group.by = "group",
                          datatype = "spatial", coordinates = spatial_coords, spatial.factors = spatial.factors)

CellChatDB <- CellChatDB.mouse
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
                              contact.dependent = TRUE, contact.range = 100)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

pathway.show <- "Glutamate"

pdf(file = "plots\\cc-com\\circle-plot.pdf",
    width = 12,
    height = 12)

netVisual_aggregate(cellchat,
                   signaling = pathway.show,
                   layout = "circle",
                   edge.weight.max = NULL,
                   edge.width.max = 8)

dev.off()

pdf(file = "plots\\cc-com\\spatial-plot.pdf",
    width = 12,
    height = 12)

netVisual_aggregate(cellchat, 
                signaling = pathway.show, 
                layout = "spatial", 
                edge.width.max = 2, 
                vertex.size.max = 1, 
                alpha.image = 0.2, 
                vertex.label.cex = 3.5)

dev.off()