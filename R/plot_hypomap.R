library(Seurat)
library(magrittr)
library(data.table)

gn.markers <- fread("markers.tsv")

so.hypomap <- readRDS("HypoMap/hypoMap.rds")
so.hypomap$GLU_LHA <- ifelse(so.hypomap$C7_named == "C7-1: GLU" & so.hypomap$Region_summarized == "Lateral hypothalamic area", T, NA)

so.hypomap.lha <- so.hypomap[,so.hypomap$GLU_LHA == T]
so.hypomap.lha %<>% FindVariableFeatures()
so.hypomap.lha %<>% RunUMAP(reduction="scvi", dims = 1:30)

## PLOT
# expression levels
library(ggplot2)
library(cowplot)
library(tximport)
t2g <- fread("gencode.vM22.metadata.MGI.gz",header=F, stringsAsFactors = F) # gencode annotations

fls <- list.files("quant",recursive = T,pattern = "quant.sf$",full.names = T)
names(fls) <- sapply(strsplit(fls,"/"),"[[",2) # get sample names
txi <- tximport(fls,"salmon",tx2gene = t2g[,1:2]) # 2 column tx2gene

ggsave2("plots/hypomap_gndetected.pdf", width = 5, height = 2,
  qplot(x=c(colSums(so.hypomap.lha@assays$RNA@counts>0), colSums(txi$counts>0)), geom="density", alpha=0.3, fill=rep(c("HypoMap GLU+LHA", "Patch-seq"), c(ncol(so.hypomap.lha), ncol(txi$counts))), xlab="Genes detected", ylab="Density") + scale_fill_brewer(palette="Set1", name="Data") + guides(alpha="none") + theme_cowplot()
)

library(data.table)
dat.det.lha <- as.data.table(reshape2::melt(cbind(  
  "HypoMap GLU+LHA" = rowMeans(so.hypomap.lha[gn.markers[type == "cell_type_specific_markers", gene]]@assays$RNA@counts>0),
  "Patch-seq" = rowMeans(txi$counts[intersect(gn.markers[type == "cell_type_specific_markers", gene], row.names(txi$counts)),]>0)
)))

ggsave2("plots/hypomap_gndetratio.pdf", width = 5, height = 10,
  ggplot(dat.det.lha, aes(x=value, y=reorder(Var1, value, max), col=Var2)) +
  geom_point() +
  labs(x="Detection rate", y=NULL) +
  scale_color_brewer(palette="Set1", name="Data") +
  theme_cowplot() +
  theme(panel.border = element_rect(color="black"))
)

# Dim reductions
library(ggplot2)
library(cowplot)
ggsave2("plots/hypomap_anno.pdf", width = 18, height = 8,
plot_grid(align = "hv", axis = "trbl",
  DimPlot(so.hypomap, group.by = "C7_named"),
  DimPlot(so.hypomap, group.by = "Region_summarized"),
  DimPlot(so.hypomap, group.by = "GLU_LHA"),
  DimPlot(so.hypomap.lha, reduction = "umap", group.by = "C185_named", raster=T) + ggtitle("GLU_LHA")
))

ggsave2("plots/hypomap_exprs.pdf", width = 12, height = 8,
  FeaturePlot(so.hypomap.lha, reduction = "umap", features = c("Pvalb", "Nppc", "Esr1", "Plpp4", "Glpr1", "Npy", "Pax6", "Gal", "Hcrt", "Avpr1a", "Slc17a6"), order = T, raster = T)
)

ggsave2("plots/hypomap_exprs_extended.pdf", width = 12, height = 42,
  FeaturePlot(so.hypomap.lha, reduction = "umap", features = gn.markers$gene, order = T, raster = T)
)

# Density
library(Nebulosa)
## FA-Bk: Pvalb, Nppc
## Burst: Esr1, Plpp4
## RS-N Esr1, Glpr1
## LS-N Npy, Pax6
## LS-W: Gal, Hcrt
## RS-W Avpr1a, Hcrt

ggsave2("plots/hypomap_density.pdf", width = 9, height = 12,
plot_grid(align = "v", axis = "rl", ncol = 1,
  plot_density(so.hypomap.lha, reduction="umap", features = c("Pvalb", "Nppc"), joint = T),
  plot_density(so.hypomap.lha, reduction="umap", features = c("Esr1", "Plpp4"), joint = T),
  plot_density(so.hypomap.lha, reduction="umap", features = c("Esr1", "Kcnab1"), joint = T), # Glpr1 not detected
  plot_density(so.hypomap.lha, reduction="umap", features = c("Npy", "Pax6"), joint = T),
  plot_density(so.hypomap.lha, reduction="umap", features = c("Gal", "Hcrt"), joint = T),
  plot_density(so.hypomap.lha, reduction="umap", features = c("Avpr1a", "Hcrt"), joint = T)
))

#
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

mat.cnt.lha <- as.matrix(so.hypomap.lha[gn.markers[type == "cell_type_specific_markers",gene],]@assays$RNA@counts)
mat.cntz.lha <- t(apply(mat.cnt.lha[rowSums(mat.cnt.lha>0) > 1000,], 1, scale))

pdf("plots/hypomap_heatmap_extended.pdf", width = 12, height = 6)
Heatmap(
  matrix = mat.cntz.lha, 
  name = "Z-score", 
  border = "black",
  col = colorRamp2(seq(-2,2, length.out=11), rev(brewer.pal(11, "PRGn"))), 
  column_split = so.hypomap.lha$C185,
  row_split = gn.markers[match(row.names(mat.cntz.lha), gene), group],
  cluster_columns = F,
  show_row_dend = F,
  show_column_names = F,
  use_raster=T
)
dev.off()

# Z-score dot plot
library(data.table)
dat.lha <- as.data.table(reshape2::melt(as.matrix(so.hypomap.lha[gn.markers[type == "cell_type_specific_markers", gene]]@assays$RNA@counts)))
dat.lha[, z := scale(value), by="Var1"]
dat.lha[, det := value > 0]
dat.lha[, avg.det := mean(det), by="Var1"]
dat.lha[, c185_named := so.hypomap.lha$C185_named[Var2]]
dat.lha[gn.markers, group := group, on = "Var1 == gene"]

library(ggplot2)
library(cowplot)

ggsave2("plots/hypomap_z.pdf", width = 5, height = 4,
  ggplot(dat.lha[avg.det > 0.01, list(avg = mean(z, na.rm=T), det = mean(det, na.rm=T) ), by=c("c185_named", "group")], aes(x=c185_named, y=group, col=avg, size=det)) + 
  geom_point() +
  labs(x="Cluster", y=NULL, size="detection", color = "scaled expression") +
  scale_color_viridis_c(limits = c(0,0.5), oob=scales::squish) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
)