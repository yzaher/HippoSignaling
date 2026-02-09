#set up the environment####
library(tidyverse)
library(Seurat)
library(ChIPseeker)
library(SeuratWrappers)
library(VariantAnnotation)
library(future)
library(biovizBase)
library(EnsDb.Mmusculus.v79)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicRanges)
library(GenomeInfoDb)
library(TFBSTools)
library(org.Mm.eg.db)
library(JASPAR2024)
library(Signac)
library(motifmatchr)
plan("multisession", workers = 16)
plan("multicore", workers = 16)
options(future.globals.maxSize = 60 * 1024^3)
set.seed(2568)
#
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
######
counts <- Read10X_h5('data/ours/Ctrl_E14.5/filtered_feature_bc_matrix.h5')
fragpath <- "data/ours/Ctrl_E14.5/atac_fragments.tsv.gz"
####
E14.5_Ct <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA")
#
E14.5_Ct[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath, 
  annotation = annotation)
#
#E14.5_Ct
##
DefaultAssay(E14.5_Ct) <- "ATAC"
E14.5_Ct <- NucleosomeSignal(E14.5_Ct)
E14.5_Ct <- TSSEnrichment(E14.5_Ct)
##
DensityScatter(E14.5_Ct, x = 'nCount_ATAC', y = 'TSS.enrichment',
               log_x = TRUE, quantiles = TRUE)
##
VlnPlot(
  object = E14.5_Ct,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0)
##
E14.5_Ct <- subset(
  x = E14.5_Ct,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1800 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1)
###
DefaultAssay(E14.5_Ct) <- "RNA"
E14.5_Ct <- NormalizeData(E14.5_Ct)
E14.5_Ct <- FindVariableFeatures(E14.5_Ct, nfeatures = 2000)
E14.5_Ct <- ScaleData(E14.5_Ct)
E14.5_Ct <- RunPCA(E14.5_Ct, reduction.name = 'mo_pca')
###
DefaultAssay(E14.5_Ct) <- "ATAC"
E14.5_Ct <- FindTopFeatures(E14.5_Ct, min.cutoff = 5)
E14.5_Ct <- RunTFIDF(E14.5_Ct)
E14.5_Ct <- RunSVD(E14.5_Ct)
###
E14.5_Ct <- FindMultiModalNeighbors(
  object = E14.5_Ct,
  reduction.list = list("mo_pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = c("RNA.weight", 'ATAC.weight'),
  verbose = TRUE)
##
E14.5_Ct <- RunUMAP(
  object = E14.5_Ct,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE)
##
E14.5_Ct$sample_ID <- 'E14.5_Ct (Chuang)'
#######
counts <- Read10X_h5('data/ours/Ctrl_E17.5/filtered_feature_bc_matrix.h5')
fragpath <- "data/ours/Ctrl_E17.5/atac_fragments.tsv.gz"
#
E17.5_Ct <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA")
#
E17.5_Ct[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath, 
  annotation = annotation)
#
#E17.5_Ct
##
DefaultAssay(E17.5_Ct) <- "ATAC"
E17.5_Ct <- NucleosomeSignal(E17.5_Ct)
E17.5_Ct <- TSSEnrichment(E17.5_Ct)
##
DensityScatter(E17.5_Ct, x = 'nCount_ATAC', y = 'TSS.enrichment',
               log_x = TRUE, quantiles = TRUE)
##
VlnPlot(
  object = E17.5_Ct,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0)
##
E17.5_Ct <- subset(
  x = E17.5_Ct,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1800 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1)
###
DefaultAssay(E17.5_Ct) <- "RNA"
E17.5_Ct <- NormalizeData(E17.5_Ct)
E17.5_Ct <- FindVariableFeatures(E17.5_Ct, nfeatures = 2000)
E17.5_Ct <- ScaleData(E17.5_Ct)
E17.5_Ct <- RunPCA(E17.5_Ct, reduction.name = 'mo_pca')
###
DefaultAssay(E17.5_Ct) <- "ATAC"
E17.5_Ct <- FindTopFeatures(E17.5_Ct, min.cutoff = 5)
E17.5_Ct <- RunTFIDF(E17.5_Ct)
E17.5_Ct <- RunSVD(E17.5_Ct)
###
E17.5_Ct <- FindMultiModalNeighbors(
  object = E17.5_Ct,
  reduction.list = list("mo_pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = c("RNA.weight", 'ATAC.weight'),
  verbose = TRUE)
##
E17.5_Ct <- RunUMAP(
  object = E17.5_Ct,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE)
E17.5_Ct$sample_ID <- 'E17.5_Ct (Chuang)'
########
counts <- Read10X_h5('data/ours/Mt_E14.5/filtered_feature_bc_matrix.h5')
fragpath <- "data/ours/Mt_E14.5/atac_fragments.tsv.gz"
#
E14.5_Mt <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA")
#
E14.5_Mt[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath, 
  annotation = annotation)
#
#E14.5_Mt
##
DefaultAssay(E14.5_Mt) <- "ATAC"
E14.5_Mt <- NucleosomeSignal(E14.5_Mt)
E14.5_Mt <- TSSEnrichment(E14.5_Mt)
##
DensityScatter(E14.5_Mt, x = 'nCount_ATAC', y = 'TSS.enrichment',
               log_x = TRUE, quantiles = TRUE)
##
VlnPlot(
  object = E14.5_Mt,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0)
##
E14.5_Mt <- subset(
  x = E14.5_Mt,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1800 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1)
###
DefaultAssay(E14.5_Mt) <- "RNA"
E14.5_Mt <- NormalizeData(E14.5_Mt)
E14.5_Mt <- FindVariableFeatures(E14.5_Mt, nfeatures = 2000)
E14.5_Mt <- ScaleData(E14.5_Mt)
E14.5_Mt <- RunPCA(E14.5_Mt, reduction.name = 'mo_pca')
###
DefaultAssay(E14.5_Mt) <- "ATAC"
E14.5_Mt <- FindTopFeatures(E14.5_Mt, min.cutoff = 5)
E14.5_Mt <- RunTFIDF(E14.5_Mt)
E14.5_Mt <- RunSVD(E14.5_Mt)
###
E14.5_Mt <- FindMultiModalNeighbors(
  object = E14.5_Mt,
  reduction.list = list("mo_pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = c("RNA.weight", 'ATAC.weight'),
  verbose = TRUE)
##
E14.5_Mt <- RunUMAP(
  object = E14.5_Mt,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE)
##
E14.5_Mt$sample_ID <- 'E14.5_Mt (Chuang)'

#########
counts <- Read10X_h5('data/ours/Mt_E17.5/filtered_feature_bc_matrix.h5')
fragpath <- "data/ours/Mt_E17.5/atac_fragments.tsv.gz"
#
E17.5_Mt <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA")
#
E17.5_Mt[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath, 
  annotation = annotation)
#
#E17.5_Mt
##
DefaultAssay(E17.5_Mt) <- "ATAC"
E17.5_Mt <- NucleosomeSignal(E17.5_Mt)
E17.5_Mt <- TSSEnrichment(E17.5_Mt)
##
DensityScatter(E17.5_Mt, x = 'nCount_ATAC', y = 'TSS.enrichment',
               log_x = TRUE, quantiles = TRUE)
##
VlnPlot(
  object = E17.5_Mt,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0)
##
E17.5_Mt <- subset(
  x = E17.5_Mt,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1800 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1)
###
DefaultAssay(E17.5_Mt) <- "RNA"
E17.5_Mt <- NormalizeData(E17.5_Mt)
E17.5_Mt <- FindVariableFeatures(E17.5_Mt, nfeatures = 2000)
E17.5_Mt <- ScaleData(E17.5_Mt)
E17.5_Mt <- RunPCA(E17.5_Mt, reduction.name = 'mo_pca')
###
DefaultAssay(E17.5_Mt) <- "ATAC"
E17.5_Mt <- FindTopFeatures(E17.5_Mt, min.cutoff = 5)
E17.5_Mt <- RunTFIDF(E17.5_Mt)
E17.5_Mt <- RunSVD(E17.5_Mt)
###
E17.5_Mt <- FindMultiModalNeighbors(
  object = E17.5_Mt,
  reduction.list = list("mo_pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = c("RNA.weight", 'ATAC.weight'),
  verbose = TRUE)
##
E17.5_Mt <- RunUMAP(
  object = E17.5_Mt,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE)
###
E17.5_Mt$sample_ID <- 'E17.5_Mt (Chuang)'


