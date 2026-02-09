#setting the environment#####
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(future)
library(ggrepel) 
library(SeuratDisk)
plan("multisession", workers = 16)
plan("multicore", workers = 16)
options(future.globals.maxSize = 16000 * 1024^2)
#load hl data####
hl <- seurat_obj <- LoadH5Seurat("C1filtered.h5seurat")
hl@meta.data$sample_ID <- 'Human_lung'
#merge all our objects to hl####
#convert gene names to capital letters
rownames(E14.5_ct) <- toupper(rownames(E14.5_ct))
rownames(E14.5_mt) <- toupper(rownames(E14.5_mt))
rownames(E17.5_ct) <- toupper(rownames(E17.5_ct))
rownames(E17.5_mt) <- toupper(rownames(E17.5_mt))
#merge into one object
ours_hl <- merge(x = hl, y = c(E14.5_ct, E14.5_mt, E17.5_ct, E17.5_mt))
#run standard workflow
ours_hl <- NormalizeData(ours_hl)
ours_hl <- FindVariableFeatures(ours_hl, nfeatures = 2000)
ours_hl <- ScaleData(ours_hl)
ours_hl <- RunPCA(ours_hl)
#integrate and join the layers
ours_hl <- IntegrateLayers(object = ours_hl, method = CCAIntegration,
                                orig.reduction = "pca", 
                                new.reduction = "cca_1",
                                verbose = FALSE)
ours_hl[["RNA"]] <- JoinLayers(ours_hl[["RNA"]])
#downstream analysis
ours_hl <- FindNeighbors(ours_hl, reduction = "cca_1", dims = 1:30)
ours_hl <- FindClusters(ours_hl, resolution = .1)
ours_hl <- RunUMAP(ours_hl, reduction = "cca_1", dims = 1:30)
#save into RDS
saveRDS(ours_hl, 'ours_with_humanlung.RDS')
