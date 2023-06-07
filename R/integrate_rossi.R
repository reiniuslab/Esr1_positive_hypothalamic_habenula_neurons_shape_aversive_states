library(data.table)
library(magrittr)

library(Seurat)

gn.markers <- fread("markers.tsv")

## Patch-seq
library(tximport)
t2g <- fread("gencode.vM22.metadata.MGI.gz",header=F, stringsAsFactors = F) # gencode annotations

fls <- list.files("quant",recursive = T,pattern = "quant.sf$",full.names = T)
names(fls) <- sapply(strsplit(fls,"/"),"[[",2) # get sample names
txi <- tximport(fls,"salmon",tx2gene = t2g[,1:2]) # 2 column tx2gene

meta <- fread("metadata_annotation.txt", stringsAsFactors = F)

so.ps <- CreateSeuratObject(txi$counts, project = "patchseq")
so.ps$ephys_category <- meta$Ephys_Category

## Rossi et al., from https://github.com/stuberlab/Rossi-et-al-2021
so.rossi <- readRDS("rossi_2021/200614_MLB017.subs.dbl.rds")
mat.fluor <- so.rossi[c("eyfp", "tdtomato")]@assays$RNA@counts
so.rossi.lhb <- so.rossi[, mat.fluor[1,] > 1 & mat.fluor[2,] == 0]

## Integrate data
ls.so <- list(
  "patchseq" = so.ps,
  "rossi" = so.rossi.lhb
)

ls.so %<>% lapply(function(x){
  x %<>% NormalizeData(verbose = FALSE)
  x %<>% FindVariableFeatures(nfeatures = 1000, verbose = FALSE)
})

int.features <- SelectIntegrationFeatures(object.list = ls.so)

ls.so %<>% lapply(function(x){
  x %<>% ScaleData(features = int.features, verbose = FALSE)
  x %<>% RunPCA(features = int.features, verbose = FALSE)
})

int.anchors <- FindIntegrationAnchors(object.list = ls.so, anchor.features = int.features, reduction = "rpca")

so.merge <- IntegrateData(anchorset = int.anchors, k.weight = 20)
DefaultAssay(so.merge) <- "integrated"
so.merge$dataset <- rep(names(ls.so), sapply(ls.so, ncol))

so.merge %<>% ScaleData(verbose = FALSE)
so.merge %<>% RunPCA(npcs = 30, verbose = FALSE)
so.merge %<>% RunUMAP(reduction = "pca", dims = 1:30)
so.merge %<>% FindNeighbors(reduction = "pca", dims = 1:30)
so.merge %<>% FindClusters(resolution = 0.5)

library(cowplot)
library(Nebulosa)
ggsave2("plots/rossi_integrated_umap.pdf", width = 10, height = 9,
  plot_grid(align = "hv", axis = "trbl", ncol = 3,
    DimPlot(so.merge, reduction = "umap", group.by = "dataset", shuffle = T),
    DimPlot(so.merge, reduction = "umap", group.by = "seurat_clusters"),
    DimPlot(so.merge, reduction = "umap", group.by = "ephys_category", cells.highlight =  lapply(1:6, function(i) which(so.merge$ephys_category == i)) ) + ggtitle("ephys_category") + scale_colour_manual(labels=c("NA",6:1), values = c("grey",rev(scales::hue_pal()(6) )) ),
    plot_density(so.merge, reduction="umap", features = "Serpinb1b"),
    plot_density(so.merge, reduction="umap", features = "Nppc"),
    plot_density(so.merge, reduction="umap", features = "Esr1"),
    plot_density(so.merge, reduction="umap", features = "Plpp4"),
    plot_density(so.merge, reduction="umap", features = "Npy"),
    plot_density(so.merge, reduction="umap", features = "Kcnab1"),
    plot_density(so.merge, reduction="umap", features = "Pax6"),
    plot_density(so.merge, reduction="umap", features = "Gal"),
    plot_density(so.merge, reduction="umap", features = "Hcrt")
  )
)
