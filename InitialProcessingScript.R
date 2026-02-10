'
This code has been gathered from multiple R files used in the very early stages
  of analysis and has been used to reproduce the results.
'
####----setting env----####
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(future)
library(ggrepel) 
plan("multisession", workers = 16)
plan("multicore", workers = 16)
options(future.globals.maxSize = 16000 * 1024^2)
setwd("/wynton/home/chuang/eyao/InitialProcessing/")
####----Local Functions----####
find_min_pc <- function(seurat_obj, reduction = 'pca'){
  stdv <- seurat_obj[[reduction]]@stdev
  percent_stdv <- (stdv/sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1] 
  co2 <- sort(which((percent_stdv[1:length(percent_stdv) - 1] - 
                       percent_stdv[2:length(percent_stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min_pc <- min(co1, co2)
  return(min_pc)
}

obj_create_process_fun <- function(h5.path, UMIs){ 
  counts <- Read10X_h5(h5.path)
  object <- CreateSeuratObject(counts, min.cells = 3, min.features = 200)
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^mt-|^MT-")
  
  # QC filtering
  object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 &
                     nCount_RNA > UMIs & percent.mt < 5)
  
  # Normalize and find variable features
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object, nfeatures = 2000)
  
  # Skip object if too few variable features
  if (length(VariableFeatures(object)) < 1500) {
    message("Skipping object: fewer than 1500 variable features (", 
            length(VariableFeatures(object)), ") in ", h5.path)
    return(NULL)
  }
  
  # Continue normal processing
  object <- ScaleData(object)
  object <- RunPCA(object)
  min_pc <- find_min_pc(object)
  object <- FindNeighbors(object, dims = 1:min_pc)
  object <- RunUMAP(object, dims = 1:min_pc)
  
  return(object)
}
#quick processing
q_process <- function(object, dims = 20){
  # Normalize and find variable features
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object)
  
  
  # Continue normal processing
  object <- ScaleData(object)
  object <- RunPCA(object)
  min_pc <- find_min_pc(object)
  object <- FindClusters(object, resolution = 0.1)
  object <- FindNeighbors(object, dims = 1:dims)
  object <- RunUMAP(object, dims = 1:dims)
  
  return(object)
}
#transfer labels from ref
ref_label <- function(query_obj = NULL, create.ref = FALSE){
  
  library(SeuratDisk)
  
  if(!exists("ref_obj", envir = .GlobalEnv)){
    path <- "/wynton/home/chuang/eyao/R/r-projects/sep_2025/CellRef/LungMAP_MouseLung_CellRef.v1.1.h5seurat"
    
    ref_obj <- LoadH5Seurat(path)
    DefaultAssay(ref_obj) <- "RNA"
    
    ref_obj <- NormalizeData(ref_obj)
    ref_obj <- FindVariableFeatures(ref_obj)
    ref_obj <- ScaleData(ref_obj)
    ref_obj <- RunPCA(ref_obj)
    
    if(create.ref){
      assign("ref_obj", ref_obj, envir = .GlobalEnv)
    }
  } else {
    ref_obj <- get("ref_obj", envir = .GlobalEnv)
  }
  
  if(create.ref && !is.null(query_obj)){
    stop("To create reference object, query_obj must be NULL.")
  }
  
  anchors <- FindTransferAnchors(
    reference = ref_obj,
    query = query_obj,
    dims = 1:20,
    reference.assay = "RNA"
  )
  
  predictions <- TransferData(
    anchorset = anchors,
    refdata = ref_obj$celltype_level3_fullname,
    dims = 1:20
  )
  
  query_obj <- AddMetaData(query_obj, predictions)
  
  return(query_obj)
}
####----Basic Processing----####
all_objects <- list()
inputs <- list.files("Filtered_h5/", full.names = T)
umis_metadata <- data.frame(path = inputs,
                     UMIs = c(5800, 5107, 3700, 3400), #input manually 
                     name = str_extract(inputs, '\\w+\\.\\w+(?=\\.h5$)'))
for(i in c(1:length(umis_metadata$path))){
  message(paste0("Now working on: ", umis_metadata$name[i]))
  all_objects[[umis_metadata$name[i]]] <- obj_create_process_fun(h5.path = umis_metadata$path[i], 
                                                                 UMIs = umis_metadata$UMIs[i])
}

####----Clustering E14.5_Ct----####
E14.5_ct <- all_objects[["Chuang_E14.5_Ct"]]
E14.5_ct <- FindClusters(E14.5_ct, resolution = 0.1)
FeaturePlot(E14.5_ct, 'Epcam', label = T)
E14.5_ct <- subset(E14.5_ct, seurat_clusters == '1') # pick epcam +ve cluster
E14.5_ct <- q_process(E14.5_ct)
FeaturePlot(E14.5_ct, # remove mesenchymal clusters
            c("Col1a1", "Col1a2", "Col3a1", "Col5a2", "Col23a1"), label = T)
E14.5_ct <- FindClusters(E14.5_ct, resolution = 0.2)
E14.5_ct <- subset(E14.5_ct, seurat_clusters != 2)
E14.5_ct <- q_process(E14.5_ct)
E14.5_ct <- FindClusters(E14.5_ct, resolution = 0.9)
FeaturePlot(E14.5_ct, label = T,
            c("Sox2", 'Sox9', 'Top2a', 'Cenpf', 'Gna14', 'Adam28', 'Pex5l', 'Cnpy1'))
E14.5_ct@meta.data <- E14.5_ct@meta.data %>% 
  mutate(final_clusters = case_when(
    seurat_clusters == 6 ~ 1,
    seurat_clusters == 2 ~ 2,
    seurat_clusters == 0 ~ 3,
    seurat_clusters == 10 ~ 4,
    TRUE ~ 0
  ))
Idents(E14.5_ct) <- 'final_clusters'
all_objects[["Chuang_E14.5_Ct"]] <- E14.5_ct
####----Clustering E14.5_Mt----####
e14.5_mt <- all_objects[["Chuang_E14.5_Mt"]]
e14.5_mt <- FindClusters(e14.5_mt, resolution = 0.1)
FeaturePlot(e14.5_mt, 'Epcam', label = T)
e14.5_mt <- subset(e14.5_mt, seurat_clusters == 1)
e14.5_mt <- q_process(e14.5_mt)
FeaturePlot(e14.5_mt, label = T,
            c("Sox9", 'Sox2', 'Top2a', 'Cenpf', 'Gna14', 'Adam28', 'Pex5l', 'Cnpy1'))
FeaturePlot(e14.5_mt, # remove mesenchymal clusters
            c("Col1a1", "Col1a2", "Col3a1", "Col5a2", "Col23a1"), label = T)
e14.5_mt <- FindClusters(e14.5_mt, resolution = 1)
e14.5_mt <- subset(e14.5_mt, seurat_clusters != 1)
e14.5_mt@meta.data <- e14.5_mt@meta.data %>% 
  mutate(final_clusters = case_when(
    seurat_clusters %in% c(4, 9, 11) ~ 0,
    seurat_clusters %in% c(12, 15) ~ 1,
    seurat_clusters %in% c(0, 3, 6, 13) ~ 2,
    seurat_clusters == 5 ~ 3,
    seurat_clusters == 14 ~ 4,
    TRUE ~ 5
  ))
Idents(e14.5_mt) <- "final_clusters"
all_objects[["Chuang_E14.5_Mt"]] <- e14.5_mt

####----Clustering E17.5_Ct----####
e17.5_ct <- all_objects[["Chuang_E17.5_Ct"]]
e17.5_ct <- FindClusters(e17.5_ct, resolution = 0.1)
FeaturePlot(e17.5_ct, 'Epcam', label = T)
e17.5_ct <- subset(e17.5_ct, seurat_clusters == 2)
e17.5_ct <- q_process(e17.5_ct)
FeaturePlot(e17.5_ct, # remove mesenchymal clusters
            c("Col1a1", "Col1a2", "Col3a1", "Col5a2", "Col23a1"), label = T)
e17.5_ct <- FindClusters(e17.5_ct, resolution = 0.1)
e17.5_ct <- subset(e17.5_ct, seurat_clusters != 2)
e17.5_ct <- q_process(e17.5_ct, dims = 30)
e17.5_ct <- FindClusters(e17.5_ct, resolution = 0.4)
FeaturePlot(e17.5_ct, label = T,
            c('Sox9', 'Sftpc', 'Ager', 'Hopx'))
e17.5_ct@meta.data <- e17.5_ct@meta.data %>% 
  mutate(final_clusters = case_when(
    seurat_clusters %in% c(2, 4) ~ 0,
    seurat_clusters == 3 ~ 1,
    seurat_clusters == 0 ~ 2,
    seurat_clusters == 1 ~ 3
    ))
Idents(e17.5_ct) <- "final_clusters"
all_objects[["Chuang_E17.5_Ct"]] <- e17.5_ct

####----Clustering E17.5_Mt----####
e17.5_mt <- all_objects[["Chuang_E17.5_Mt"]]
e17.5_mt <- FindClusters(e17.5_mt, resolution = 0.1)
FeaturePlot(e17.5_mt, 'Epcam', label = T)
e17.5_mt <- subset(e17.5_mt, seurat_clusters == 1)
e17.5_mt <- q_process(e17.5_mt)
ref_label(create.ref = T)
e17.5_mt <- ref_label(query_obj = e17.5_mt)
DimPlot(e17.5_mt, group.by = 'predicted.id', label = T) + NoLegend()
epithelial_cells <- c( "Alveolar type 1 cell", "AT1/AT2", "Ciliated cell",
                       "Deuterosomal cell", "Pulmonary neuroendocrine cell", 
                       "Secretory cell", "Sox9+/Id2+ epithelial cell")
e17.5_mt <- subset(e17.5_mt, predicted.id %in% epithelial_cells)
e17.5_mt <- q_process(e17.5_mt, dims = 25)
e17.5_mt <- FindClusters(e17.5_mt, resolution = 0.1)
FeaturePlot(e17.5_mt, label = T,
            c( 'Rps26', 'Pfdn5', 'Kank4', 'Chrm3'))
e17.5_mt <- FindSubCluster(e17.5_mt, cluster = 0, graph.name = 'RNA_snn',
                           resolution = 0.2)
e17.5_mt@meta.data <- e17.5_mt@meta.data %>% 
  mutate(final_clusters = case_when(
    sub.cluster == 4 ~ '0', 
    sub.cluster == 5 ~ '1', 
    sub.cluster == 9 ~ '2', 
    sub.cluster == '0_1' ~ '3', 
    sub.cluster == 1 ~ '5', 
    sub.cluster == '0_0' ~ '6', 
    sub.cluster %in% c(2, 7) ~ '4',
    TRUE ~ 'Unassigned'
  ))
e17.5_mt@meta.data$final_clusters <- factor(e17.5_mt@meta.data$final_clusters,
                                            levels = c('0', '1', '2', '3', '4', 
                                                       '5', '6', 'Unassigned'))
Idents(e17.5_mt) <- "final_clusters"
all_objects[["Chuang_e17.5_mt"]] <- e17.5_mt
