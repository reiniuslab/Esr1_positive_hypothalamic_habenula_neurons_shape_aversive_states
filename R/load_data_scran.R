## read salmon output
library(tximport)
library(data.table)
library(ggplot2)
library(magrittr) # enables %<>% instead of x <- fun(x)

# t10 palette for plotting
pal.t10 <- c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D", "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2", "#CDCC5D", "#6DCCDA")

t2g <- fread("gencode.vM22.metadata.MGI.gz",header=F, stringsAsFactors = F) # gencode annotations

fls <- list.files("quant",recursive = T,pattern = "quant.sf$",full.names = T)
names(fls) <- sapply(strsplit(fls,"/"),"[[",2) # get sample names
txi <- tximport(fls,"salmon",tx2gene = t2g[,1:2]) # 2 column tx2gene

## add metadata
meta <- fread("metadata_annotation.txt", stringsAsFactors = F)
meta[,Batch:=sapply(strsplit(NGI_Sample,"_"),"[[",1)]
meta[,Animal:=sapply(strsplit(Animal_Cell_ID,"A|C"),"[[",2)]
meta[,Cell:=sapply(strsplit(Animal_Cell_ID,"C"),"[[",2)]

## convert to SCE object
library(igraph)
library(scater)
library(scran)
sce <- SingleCellExperiment(assays=list(counts=txi$counts,tpm=txi$abundance),colData=meta)
sce$Ephys_Category <- factor(sce$Ephys_Category)

# identify outliers based on qc metrics
is.expr <- calcAverage(sce)>=1
is.mito <- grepl("^mt-",rownames(sce))

sce %<>% calculateQCMetrics(feature_controls = list(Mt=is.mito),use_spikes = F)
#sce %<>% runPCA(use_coldata = T, detect_outliers = T,selected_variables=c("total_features_by_counts","total_counts","pct_counts_feature_control"))
drop.counts <- isOutlier(sce$total_counts, nmads = 3, log = T, batch = sce$Batch, type = "lower")
drop.features <- isOutlier(sce$total_features_by_counts, nmads = 3, log = T, batch = sce$Batch, type = "lower")
drop.mito <- isOutlier(sce$pct_counts_feature_control, nmads = 3, batch = sce$Batch, type = "higher")
sce$outlier <- drop.counts | drop.features | drop.mito

# remove outliers and non-expressed genes, normalise data
sce.filt <- subset(sce,is.expr,outlier==F)
sce.filt %<>% computeSumFactors()
sce.filt %<>% normalize()

# identify hvgs
var.fit <- trendVar(sce.filt, use.spikes=F, block=sce.filt$Batch)
var.decomp <- decomposeVar(sce.filt, var.fit)
var.decomp$name <- rownames(var.decomp)
var.decomp <- var.decomp[with(var.decomp,order(-bio,FDR)),]
hvgs <- var.decomp[var.decomp$FDR < 0.05,'name']

# add scaled data to assays
assays(sce.filt)[["scaled_logcounts"]] <- t(apply(assays(sce.filt)[["logcounts"]],1,scale))
assays(sce.filt)[["scaled_tpm"]] <- t(apply(assays(sce.filt)[["tpm"]],1,scale))

# run dimensionality reduction
# set.seed(53) # for UMAP reproducibility
sce.filt %<>% runUMAP(feature_set = head(hvgs,500))
# set.seed(33) # for t-SNE reproducibility
sce.filt %<>% runTSNE(feature_set = head(hvgs,500))

# identify clusters
sce.filt$Cluster <- quickCluster(sce.filt,use.ranks=T,min.size=ncol(sce.filt)*0.1,subset.row=head(hvgs,500),graph.fun=cluster_louvain)

## cleanup
rm(list = c("fls","is.expr","t2g","txi"))
