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

## Hypomap, from https://www.repository.cam.ac.uk/handle/1810/340518
so.hypomap <- readRDS("HypoMap/hypoMap.rds")
so.hypomap.lha <- subset(so.hypomap, subset= C7_named == "C7-1: GLU" & Region_summarized == "Lateral hypothalamic area")

## Integrate data
ls.so.all <- list(
  "patchseq" = so.ps,
  "rossi" = so.rossi.lhb,
  "hypomap" = so.hypomap.lha
)

ls.so.all %<>% lapply(function(x){
  x %<>% NormalizeData(verbose = FALSE)
  x %<>% FindVariableFeatures(nfeatures = 1000, verbose = FALSE)
})

int.features <- SelectIntegrationFeatures(object.list = ls.so.all)

ls.so.all %<>% lapply(function(x){
  x %<>% ScaleData(features = int.features, verbose = FALSE)
  x %<>% RunPCA(features = int.features, verbose = FALSE)
})

int.anchors <- FindIntegrationAnchors(object.list = ls.so.all, anchor.features = int.features, reduction = "rpca")

so.all <- IntegrateData(anchorset = int.anchors)
DefaultAssay(so.all) <- "integrated"
so.all$dataset <- rep(names(ls.so.all), sapply(ls.so.all, ncol))

so.all %<>% ScaleData(verbose = FALSE)
so.all %<>% RunPCA(npcs = 30, verbose = FALSE)
so.all %<>% RunUMAP(reduction = "pca", dims = 1:30)
so.all %<>% FindNeighbors(reduction = "pca", dims = 1:30)
so.all %<>% FindClusters(resolution = 0.5)

library(cowplot)
library(Nebulosa)
ggsave2("plots/hypomap_integrated_umap.pdf", width = 18, height = 18,
  plot_grid(align = "hv", axis = "trbl", ncol = 3,
    DimPlot(so.all, reduction = "umap", group.by = "dataset", cells.highlight = list("patchseq" = which(so.all$dataset == "patchseq"), "rossi" = which(so.all$dataset == "rossi") ) ) + ggtitle("dataset") + scale_colour_manual(labels=c("Hypomap",'Rossi et al.', 'Patch-seq'), values = c("grey","#E69F00", "#56B4E9") ),
    DimPlot(so.all, reduction = "umap", group.by = "seurat_clusters", label = T) + theme(legend.position = 'none'),
    DimPlot(so.all, reduction = "umap", group.by = "ephys_category", cells.highlight =  lapply(1:6, function(i) which(so.all$ephys_category == i)) ) + ggtitle("ephys_category") + scale_colour_manual(labels=c("NA",6:1), values = c("grey",rev(scales::hue_pal()(6) )) ),
    plot_density(so.all, reduction="umap", features = "Pvalb"),
    plot_density(so.all, reduction="umap", features = "Nppc"),
    plot_density(so.all, reduction="umap", features = "Esr1"),
    plot_density(so.all, reduction="umap", features = "Plpp4"),
    plot_density(so.all, reduction="umap", features = "Npy"),
    plot_density(so.all, reduction="umap", features = "Kcnab1"),
    plot_density(so.all, reduction="umap", features = "Pax6"),
    plot_density(so.all, reduction="umap", features = "Gal"),
    plot_density(so.all, reduction="umap", features = "Hcrt")
  )
)
