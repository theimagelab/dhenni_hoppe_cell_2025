library("Seurat"); library("sctransform"); library("dplyr"); library("RColorBrewer")

# Read in and integrate datasets 
scs_obj <- readRDS(file="seuobj_normalization/scs_macs_nmfclusters.rds")
rama_obj <- readRDS(file="BCellOnlyIntron2.Donor14Nov.rds")
rama_obj@meta.data$barcode <- rownames(rama_obj@meta.data)
DefaultAssay(scs_obj) <- "RNA"
DefaultAssay(rama_obj) <- "alra"
scs_obj<-NormalizeData(scs_obj)
scs_obj<-FindVariableFeatures(scs_obj)
# rama "alra" assay has already been lognorm + variable faeatures
obj_list <- list(scs_obj, rama_obj)

features <- SelectIntegrationFeatures(object.list = obj_list)
anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features)
rama_scs <- IntegrateData(anchorset = anchors)

DefaultAssay(rama_scs) <- "integrated"

x<-rama_scs@meta.data
old_meta<-x
# Integration has messed up barcode order for SCS macs but not for rama's. Fixing here.
# This was obvious because the nmf labels on the umap were completely randomly positioned.
x %>% mutate(barcode_1 = case_when(
    is.na(TissueWithNaiveCombined) ~ paste0(barcode, "_1"),
    !is.na(TissueWithNaiveCombined) ~ paste0(barcode, "_2")
)) -> x
x$barcode_2<-rownames(x)
ind2<-match(x$barcode_1, x$barcode_2)
x_2<-x[ind2,]
rownames(x_2)<-x_2$barcode_1

rama_scs@meta.data <- x_2

rama_scs <- ScaleData(rama_scs, verbose = FALSE)
rama_scs <- RunPCA(rama_scs, npcs = 30, verbose = FALSE)
rama_scs <- RunUMAP(rama_scs, reduction = "pca", dims = 1:30)
rama_scs <- FindNeighbors(rama_scs, reduction = "pca", dims = 1:30)
rama_scs <- FindClusters(rama_scs, resolution = 0.5)
DimPlot(rama_scs, group.by = "nmf_rank")
DimPlot(rama_scs, group.by = "TissueWithNaiveCombined")
# ok x_2 is the correct one

# Cluster annotation merge into 1 column and colour palette
x <- rama_scs@meta.data
x$merged_cluster <- NA
x[!is.na(x$TissueWithNaiveCombined), "merged_cluster"]<- x[!is.na(x$TissueWithNaiveCombined), "TissueWithNaiveCombined"] 
x[!is.na(x$nmf_rank), "merged_cluster"]<- x[!is.na(x$nmf_rank), "nmf_rank"] 
rama_clusters <- x$TissueWithNaiveCombined %>% unique() %>% na.omit() 
nmf_ranks <- x$nmf_rank %>% unique() %>% na.omit() 
rama_colors <- colorRampPalette(brewer.pal(4, "Greys"))(length(rama_clusters)) 
names(rama_colors) <- rama_clusters
nmf_colors <- colorRampPalette(brewer.pal(10, "Spectral"))(length(nmf_ranks))
names(nmf_colors) <- nmf_ranks
x$color_nmf <- nmf_colors[x$nmf_rank]
x$color_rama <- rama_colors[x$TissueWithNaiveCombined]
x$color <- NA
x[!is.na(x$color_rama), "color"]<- x[!is.na(x$color_rama), "color_rama"] 
x[!is.na(x$color_nmf), "color"]<- x[!is.na(x$color_nmf), "color_nmf"] 
rama_scs@meta.data <- x
cols_vector <- c(nmf_colors, rama_colors)
###
DimPlot(rama_scs, group.by = "integrated_snn_res.0.5", label=TRUE)
DimPlot(rama_scs, group.by = "nmf_rank")
DimPlot(rama_scs, group.by = "merged_cluster", cols = cols_vector, label=TRUE)

Idents(rama_scs)<-"merged_cluster"
DefaultAssay(rama_scs) <- "RNA"
# Figuring out which cluster is the SSMs.
FeaturePlot(rama_scs, features=c("Ltb", "Lat", "Cxcr6", "Il18r1", "Ccl5", "Adgre1", "Siglec1"))
FeaturePlot(rama_scs, features=c("Mertk", "Chil3"))
# conclusion for SCS mac annotation:
# SSMs are the nub at bottom left of umap.
# NMF cluster 7 contains them. 
# But best resolved by integrated_snn_res.0.5 cluster 8.
# Use integrated cluster 8 + rama's clustering.
x %>% mutate(integrated_snn_res.0.5 = case_when(
    integrated_snn_res.0.5 == "5" ~ "5_SSM",
    TRUE ~ integrated_snn_res.0.5)
) -> x
rama_scs@meta.data <- x

# Cluster annotation merge into 1 column (again) but use seurat labels, not NMF labels.
x <- rama_scs@meta.data
x$merged_cluster_scs <- NA
x[!is.na(x$TissueWithNaiveCombined), "merged_cluster_scs"]<- x[!is.na(x$TissueWithNaiveCombined), "TissueWithNaiveCombined"] 
x[is.na(x$TissueWithNaiveCombined), "merged_cluster_scs"]<- x[is.na(x$TissueWithNaiveCombined), "integrated_snn_res.0.5"] 
rama_scs@meta.data <- x
DimPlot(rama_scs, group.by = "merged_cluster_scs", label=TRUE, pt.size=4)

rama_clusters <- x$TissueWithNaiveCombined %>% unique() %>% na.omit() 
seurat_clusters <- x$integrated_snn_res.0.5 %>% unique() %>% na.omit() 
rama_colors <- colorRampPalette(brewer.pal(4, "Set1"))(length(rama_clusters)) 
names(rama_colors) <- rama_clusters
seurat_colors <- colorRampPalette(brewer.pal(10, "Spectral"))(length(seurat_clusters))
names(seurat_colors) <- seurat_clusters
x$color_seurat <- seurat_colors[x$integrated_snn_res.0.5]
x$color_rama <- rama_colors[x$TissueWithNaiveCombined]
x$color <- NA
x[!is.na(x$color_rama), "color"]<- x[!is.na(x$color_rama), "color_rama"] 
x[is.na(x$color_rama), "color"]<- x[is.na(x$color_rama), "color_seurat"] 
rama_scs@meta.data <- x
cols_vector <- c(seurat_colors, rama_colors)
DimPlot(rama_scs, group.by = "merged_cluster_scs", label=TRUE, pt.size=4, cols = cols_vector)
Idents(rama_scs)<-"merged_cluster_scs"
save(rama_scs, file="rama_scs.RData")
# DimPlot(rama_scs, group.by = "merged_cluster_scs", cols = cols_vector, label=TRUE, repel=TRUE)
# 
# # Diff Exp
DefaultAssay(rama_scs) <- "RNA"
mk_no_latent_vars <- FindAllMarkers(rama_scs)
rama_scs@meta.data %>% mutate(expt = case_when(
    !is.na(color_rama) ~ "rama",
    !is.na(color_nmf) ~ "scs",
    TRUE ~ NA)
) -> rama_scs@meta.data
mk_latent_vars<-FindAllMarkers(rama_scs, test.use = "LR", latent.vars = "expt")

save(list(rama_scs, mk_latent_vars), file="rama_scs.RData")

load("rama_scs.RData")
