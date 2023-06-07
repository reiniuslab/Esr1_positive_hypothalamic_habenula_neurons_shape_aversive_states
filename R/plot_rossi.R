library(Seurat)
library(magrittr)

gn.markers <- fread("markers.tsv")

so.rossi <- readRDS("rossi_2021/200614_MLB017.subs.dbl.rds")

mat.fluor <- so.rossi[c("eyfp", "tdtomato")]@assays$RNA@counts

#so.rossi.lhb <- so.rossi[, (log2(mat.fluor[1,]+1) - log2(mat.fluor[2,]+1)) > 1]

so.rossi.lhb <- so.rossi[, mat.fluor[1,] > 1 & mat.fluor[2,] == 0]
so.rossi.lhb %<>% NormalizeData()
so.rossi.lhb %<>% FindVariableFeatures(nfeatures=1000)
so.rossi.lhb %<>% ScaleData(vars.to.regress = "percentMito")
so.rossi.lhb %<>% RunPCA()
so.rossi.lhb %<>% RunUMAP(reduction="pca", dims=1:15)
so.rossi.lhb %<>% FindNeighbors(dims = 1:15, k.param = 10)
so.rossi.lhb %<>% FindClusters(resolution = 1)

## PLOT
# dim reductions

ggsave2("plots/rossi_anno.pdf", width = 4, height = 3,
  DimPlot(so.rossi.lhb)
)

ggsave2("plots/rossi_exprs_slc.pdf", width = 8, height = 3,
 FeaturePlot(so.rossi.lhb, features = c("Slc17a6", "Slc32a1"), order=T)
)

ggsave2("plots/rossi_exprs_extended.pdf", width = 12, height = 42,
  FeaturePlot(so.rossi.lhb, features = gn.markers$gene, order=T)
)

# density
library(Nebulosa)

ggsave2("plots/rossi_density.pdf", width = 9, height = 12,
plot_grid(align = "v", axis = "rl", ncol = 1,
  plot_density(so.rossi.lhb, reduction="umap", features = c("Pvalb", "Nppc"), joint = T), # Sostdc1
  plot_density(so.rossi.lhb, reduction="umap", features = c("Esr1", "Plpp4"), joint = T),
  plot_density(so.rossi.lhb, reduction="umap", features = c("Esr1", "Kcnab1"), joint = T), # Glpr1 not detected
  plot_density(so.rossi.lhb, reduction="umap", features = c("Npy", "Pax6"), joint = T),
  plot_density(so.rossi.lhb, reduction="umap", features = c("Gal", "Hcrt"), joint = T),
  plot_density(so.rossi.lhb, reduction="umap", features = c("Avpr1a", "Hcrt"), joint = T)
))

# heatmap
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

mat.cnt.lhb <- as.matrix(so.rossi.lhb[gn.markers[type == "cell_type_specific_markers",gene],]@assays$RNA@counts)
mat.cntz.lhb <- t(apply(mat.cnt.lhb[rowSums(mat.cnt.lhb>0) > 10,], 1, scale))

pdf("plots/rossi_heatmap_extended.pdf", width = 6, height = 6)
Heatmap(
  matrix = mat.cntz.lhb, 
  name = "Z-score", 
  border = "black",
  col = colorRamp2(seq(-2,2, length.out=11), rev(brewer.pal(11, "PRGn"))), 
  column_split = so.rossi.lhb$seurat_clusters,
  row_split = gn.markers[match(row.names(mat.cntz.lhb), gene), group],
  cluster_columns = F,
  show_column_names = F, 
  show_row_dend = F,
  use_raster=T
)
dev.off()

# Z-score dot plot
library(data.table)
dat.lhb <- as.data.table(reshape2::melt(as.matrix(so.rossi.lhb[gn.markers[type == "cell_type_specific_markers", gene]]@assays$RNA@counts)))
dat.lhb[, z := scale(value), by="Var1"]
dat.lhb[, det := value > 0]
dat.lhb[, avg.det := mean(det), by="Var1"]
dat.lhb[, clusters := so.rossi.lhb$seurat_clusters[Var2]]
dat.lhb[gn.markers, group := group, on = "Var1 == gene"]

library(ggplot2)
library(cowplot)

ggsave2("plots/rossi_z.pdf", width = 4, height = 2.5,
  ggplot(dat.lhb[avg.det > 0.01, list(avg = mean(z, na.rm=T), det = mean(det, na.rm=T) ), by=c("clusters", "group")], aes(x=clusters, y=group, col=avg, size=det)) + 
  geom_point() +
  labs(x="Cluster", y=NULL, size="detection", color = "scaled expression") +
  scale_color_viridis_c(limits = c(0,0.5), oob=scales::squish) +
  theme_cowplot()
)