# Following the tutorial from https://www.youtube.com/watch?v=5HBzgsz8qyk 

# Load the libraries
library(Seurat)
library(tidyverse)
library("hdf5r")


# load the NSCLC dataset
nsclc.sparse.m <- Read10X_h5(filename = "../data/20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5")
str(nsclc.sparse.m)
cts <- nsclc.sparse.m$`Gene Expression` # Count matrix


# Initialize the Seurat object with raw, non-normalized data
nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "NSCLC",
                                       min.cells = 3, min.features = 200)

str(nsclc.seurat.obj)
nsclc.seurat.obj

# 1. Quality Control
View(nsclc.seurat.obj@meta.data)

# % MT reads: The cells that are less viable or dead have higher percentage of 
# mitochondrial activity. We can use this information to filter out cells whose 
# reads will not be useful in our analyses. 
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
View(nsclc.seurat.obj@meta.data)

VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

# 2. Filtering
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500
                           & percent.mt < 5)

# 3. Normalize data
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)
str(nsclc.seurat.obj) # to view what processes the object has been through

# 4. Identify highly variable features
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)

# Identify 10 most highly variable genes
top10 <- head(VariableFeatures(nsclc.seurat.obj), 10)

# PLot the variable features with and without labels
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel=TRUE)

# 5. Scaling.
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)

str(nsclc.seurat.obj)

# 6. Perform Linear Dimensionality Reduction
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))

# Visualize PCA results
print(nsclc.seurat.obj[['pca']], dims = 1:5, nfeatures = 5)
DimHeatmap(nsclc.seurat.obj, dims = 1, cells=500, balanced=T)

# Determine the dimensionality of the data
ElbowPlot(nsclc.seurat.obj)

# 7. Clustering 
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15)

# Understanding resolution
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
View(nsclc.seurat.obj@meta.data)

DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.5", label = T)

# setting identity of clusters
Idents(nsclc.seurat.obj)
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1"
Idents(nsclc.seurat.obj)

# Non Linear Dimensionality reduction
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)

DimPlot(nsclc.seurat.obj, reduction = "umap")
