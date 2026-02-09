####-----------SET THE ENVIRONMENT-----------####
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(future)
library(ggrepel) 
library(batchelor)
plan("multicore", workers = 16) # 16 CORES
options(future.globals.maxSize = 250000 * 1024^2) # 250 GB RAM
set.seed(49)
####-----------INPUT VARIABLES-----------####
#paths
input_h5_dir <- "/wynton/home/chuang/eyao/R/r-projects/nature_oct26/h5_files"
output_dir <- "/wynton/home/chuang/eyao/R/r-projects/nature_oct26/output_merged_objects"
metadata.path <- "/wynton/home/chuang/eyao/R/r-projects/nature_oct26/geo_data_table.csv"
LungMAP_MouseLung_CellRef_h5seurat.path <- "/wynton/home/chuang/eyao/R/r-projects/sep_2025/CellRef/LungMAP_MouseLung_CellRef.v1.1.h5seurat"
#inputs
age_group_1 <- c("E12.5", "E14.5", "E15", "E15.5")
#age_group_2 <- c("E15", "E15.5", "E16.5", "E17.5", "E18.5")
ages <- c(
  "E12.5",
  "E14.5",
  "E15",
  "E15.5",
  "E16.5",
  "E17.5",
  "E18.5"
)
age_groups_list <- list(age_group_1 = ages, age_group_2 = age_group_2)
####-----------CONFIGURE REQUIRMENTS-----------####
dir.create(output_dir, recursive = TRUE)
metadata <- read_csv(metadata.path)
metadata$identifier <- paste0(metadata$Author, "_", metadata$Run)
#set reference object
ref_obj <- LoadH5Seurat(LungMAP_MouseLung_CellRef_h5seurat.path)
DefaultAssay(ref_obj) <- "RNA"
ref_obj <- NormalizeData(ref_obj)
ref_obj <- FindVariableFeatures(ref_obj)
ref_obj <- ScaleData(ref_obj)
ref_obj <- RunPCA(ref_obj)
####-----------LOCAL FUNCTIONS-----------####
#function to identify minimum number of PCs that explains > 90% variance
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
#seurat object creation and standard processing function
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
    
  # Continue normal processing
  object <- ScaleData(object)
  object <- RunPCA(object)
  min_pc <- find_min_pc(object)
  object <- FindNeighbors(object, dims = 1:min_pc)
  object <- RunUMAP(object, dims = 1:min_pc)
  return(object)
}
####-----------RUN ALL IN A FOR FUNCTION-----------####
for (age_group in names(age_groups_list)){
  group <- age_groups_list[[age_group]]
  ##start processing initial objects
  samples_to_run_g1 <- metadata[metadata$Age %in% group,,drop=FALSE]
  group_1_list <- lapply(samples_to_run_g1$Run, function(input_srr){
    message(paste0("working on", input_srr))
    md_row <- samples_to_run_g1[samples_to_run_g1$Run==input_srr,,drop=FALSE]
    loc_h5.path <- paste0(input_h5_dir,"/", input_srr, ".h5")
    object <- obj_create_process_fun(h5.path = loc_h5.path,
                                     UMIs = md_row$min_UMIs)
    object$Epcam.exp <- FetchData(object, var = "Epcam")
    object <- subset(object, Epcam.exp > 0)
    colnames(object) <- paste0(colnames(object), md_row$identifier)
    for (column in colnames(md_row)) {
      object[[column]] <- md_row[[column]]
    }
    object$n_cell <- length(Cells(object))
    return(object)
  })
  names(group_1_list) <- samples_to_run_g1$identifier
  ######cut here#####
  #merge and process
  merged_g1 <- purrr::reduce(group_1_list, merge)
  merged_g1 <- subset(merged_g1, Epcam.exp > 0)
  merged_g1 <- NormalizeData(merged_g1)
  merged_g1 <- FindVariableFeatures(merged_g1, nfeatures = 2000)
  merged_g1 <- ScaleData(merged_g1)
  merged_g1 <- RunPCA(merged_g1)
  #integrate
  merged_g1 <- IntegrateLayers(object = merged_g1,
                               method = FastMNNIntegration,
                               orig.reduction = "pca",
                               new.reduction = "mnn_1",
                               verbose = FALSE)
  merged_g1[["RNA"]] <- JoinLayers(merged_g1[["RNA"]])
  merged_g1 <- FindNeighbors(merged_g1, reduction = "mnn_1", dims = 1:30)
  merged_g1 <- RunUMAP(merged_g1, reduction = "mnn_1", dims = 1:30)
  #annotate merged object
   anchors <- FindTransferAnchors(reference = ref_obj,
                                 query = merged_g1,
                                  dims = 1:25,
                                   reference.assay = "RNA")
  
    predictions <- TransferData(anchorset = anchors,
                               refdata = ref_obj$celltype_level3_fullname,
                                dims = 1:25)
  merged_g1 <- AddMetaData(merged_g1, metadata = predictions)
  #save merged object
  saveRDS(merged_g1, paste0(output_dir, "/", paste0(group, collapse = "_"), "_merged", ".RDS"))
}

