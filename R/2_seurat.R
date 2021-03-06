library(Seurat)
library(clustree)
library(tidyverse)

#load data
load(file.path("data", "prdata.RData"))

#run seurat
all <- CreateSeuratObject(counts = prdata,
                          min.cells = 5,
                          min.features = 500)

#regress out mitochondrial
all <- PercentageFeatureSet(all, pattern = "^mt-", col.name = "percent.mt")

# run sctransform
all <- SCTransform(all, vars.to.regress = "percent.mt", verbose = FALSE)

#pca etc
all <- RunPCA(all, verbose = FALSE)

#run elbow plot to find most relevant PCs
ElbowPlot(all)

#run UMAP and clustering on chosen PCs
all <- RunUMAP(all, dims = 1:15, verbose = FALSE)
all <- FindNeighbors(all, dims = 1:15, verbose = FALSE)
all <- FindClusters(all, resolution = seq(from=0.2,to=2,by=0.2))

#run clustree to find a cluster resolution
#url https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html#seurat-objects
clustree(all)

all <- FindClusters(all, resolution = .6)

DimPlot(all, label = TRUE) + NoLegend()

#save data
save(all, file = file.path("data", "seurat.RData"))
