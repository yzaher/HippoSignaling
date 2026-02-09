#setting the environment#####
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(future)
library(ggrepel) 
setwd("/media/youssef/My_drive/my_R_packages/February/")
plan("multisession", workers = 16)
plan("multicore", workers = 16)
options(future.globals.maxSize = 16000 * 1024^2)

#loading all seurat objects####
E14.5_ct <- read_rds('compressed/our_clustered_objects/cl_1_E14.5_Ct.RDS')
E14.5_mt <- read_rds('epcam_Mt_E14.5.Rds')
E17.5_ct <- read_rds('our_updated_objects/epcam_Ct_E17.5.Rds')
E17.5_mt <- read_rds('our_updated_objects/epcam_Mt_E17.5.Rds')
#finding markers between samples####
all <- read_rds('our_at1s_merged/ours_at1s.RDS')
Idents(all) <- 'sample_ID'
at1_14.5mt_17.5ct_marks <- FindMarkers(all, ident.1 = 'E14.5_Mt (Chuang)', ident.2 = 'E17.5_Ct (Chuang)')%>% filter(p_val_adj < 0.05)
at1_17.5ct_17.5mt_marks <-  FindMarkers(all, ident.1 = 'E17.5_Ct (Chuang)', ident.2 = 'E17.5_Mt (Chuang)')%>% filter(p_val_adj < 0.05)
#finding E14.5_ct cluster markers####
E14.5_ct_0_1_marks <- FindMarkers(E14.5_ct, ident.1 = 0, ident.2 = 1)    %>% filter(p_val_adj < 0.05)
E14.5_ct_0_2_marks <- FindMarkers(E14.5_ct, ident.1 = 0, ident.2 = 2)    %>% filter(p_val_adj < 0.05)
E14.5_ct_1_3_marks <- FindMarkers(E14.5_ct, ident.1 = 1, ident.2 = 3)    %>% filter(p_val_adj < 0.05)
E14.5_ct_2_3_marks <- FindMarkers(E14.5_ct, ident.1 = 2, ident.2 = 3)    %>% filter(p_val_adj < 0.05)
E14.5_ct_1_2_marks <- FindMarkers(E14.5_ct, ident.1 = 1, ident.2 = 2)    %>% filter(p_val_adj < 0.05)
#finding E14.5_mt cluster markers####
E14.5_mt_01_marks <- FindMarkers(E14.5_mt, ident.1 = 0 , ident.2 = 1)    %>% filter(p_val_adj < 0.05)
E14.5_mt_02_marks <- FindMarkers(E14.5_mt, ident.1 = 0 , ident.2 = 2)    %>% filter(p_val_adj < 0.05)
E14.5_mt_13_marks <- FindMarkers(E14.5_mt, ident.1 = 1 , ident.2 = 3)    %>% filter(p_val_adj < 0.05)
E14.5_mt_25_marks <- FindMarkers(E14.5_mt, ident.1 = 2 , ident.2 = 5)    %>% filter(p_val_adj < 0.05)
#finding E17.5_ct cluster markers####
E17.5_ct_01_marks <- FindMarkers(E17.5_ct, ident.1 = 0 , ident.2 = 1)    %>% filter(p_val_adj < 0.05)
E17.5_ct_12_marks <- FindMarkers(E17.5_ct, ident.1 = 1 , ident.2 = 2)    %>% filter(p_val_adj < 0.05)
E17.5_ct_13_marks <- FindMarkers(E17.5_ct, ident.1 = 1 , ident.2 = 3)    %>% filter(p_val_adj < 0.05)
E17.5_ct_23_marks <- FindMarkers(E17.5_ct, ident.1 = 2 , ident.2 = 3)    %>% filter(p_val_adj < 0.05)
#finding E17.5_mt cluster markers####
E17.5_mt_0_2 <- FindMarkers(E17.5_mt, ident.1 = 0, ident.2 = 2, group.by = 'final')  %>% filter(p_val_adj < 0.05)  
E17.5_mt_0_1 <- FindMarkers(E17.5_mt, ident.1 = 0, ident.2 = 1, group.by = 'final')  %>% filter(p_val_adj < 0.05)
E17.5_mt_1_2 <- FindMarkers(E17.5_mt, ident.1 = 1, ident.2 = 2, group.by = 'final')  %>% filter(p_val_adj < 0.05)
E17.5_mt_1_3 <- FindMarkers(E17.5_mt, ident.1 = 1, ident.2 = 3, group.by = 'final')  %>% filter(p_val_adj < 0.05)
E17.5_mt_3_4 <- FindMarkers(E17.5_mt, ident.1 = 3, ident.2 = 4, group.by = 'final')  %>% filter(p_val_adj < 0.05)
E17.5_mt_3_5 <- FindMarkers(E17.5_mt, ident.1 = 3, ident.2 = 5, group.by = 'final')  %>% filter(p_val_adj < 0.05)
E17.5_mt_3_6 <- FindMarkers(E17.5_mt, ident.1 = 3, ident.2 = 6, group.by = 'final')  %>% filter(p_val_adj < 0.05)
#find E14.5_ct cluster colors####
plot <- DimPlot(E14.5_ct, label = TRUE)
colors_used_e14.5ct <- unique(ggplot_build(plot)$data[[1]]$colour) # 
#find E14.5_mt cluster colors####
i <- DimPlot(E14.5_mt)
colors_E14.5_mt <- rownames(table(ggplot_build(i)$data[[1]]$colour)) # cluster 2 3 4 1 5 0 
#find E17.5_ct clusters colors####
p <- DimPlot(E17.5_ct)
colors_e17.5_ct <- rownames(table(ggplot_build(p)$data[[1]]$colour)) # cluster 2 1 3 0
#find E17.5_mt cluster colors####
o <- DimPlot(E17.5_mt)
colors_E17.5_mt <- rownames(table(ggplot_build(o)$data[[1]]$colour)) # cluster 4 3 5 2 6 1 0 7 
#plotting samples markers####
cr_plot_degs(at1_14.5mt_17.5ct_marks, color_to_the_right = "#4477aa", color_to_the_left = "#ee6677")
cr_plot_degs(at1_17.5ct_17.5mt_marks, color_to_the_right = "#2ca02c", color_to_the_left = "#4477aa")
write.csv(at1_14.5mt_17.5ct_marks, "pub_output/at1_14.5mt_17.5ct_marks.csv")
write.csv(at1_17.5ct_17.5mt_marks, "pub_output/at1_17.5ct_17.5mt_marks.csv")
#plotting E14.5_ct markers####
#reverse ifelse in the function(cr_plot_degs) to -1####
cr_plot_degs(E14.5_ct_0_1_marks, color_to_the_right = colors_used_e14.5ct[5], color_to_the_left = colors_used_e14.5ct[1])
cr_plot_degs(E14.5_ct_0_2_marks, color_to_the_right = colors_used_e14.5ct[4], color_to_the_left = colors_used_e14.5ct[1])
cr_plot_degs(E14.5_ct_1_3_marks, color_to_the_right = colors_used_e14.5ct[2], color_to_the_left = colors_used_e14.5ct[5])
cr_plot_degs(E14.5_ct_2_3_marks, color_to_the_right = colors_used_e14.5ct[2], color_to_the_left = colors_used_e14.5ct[4])
cr_plot_degs(E14.5_ct_1_2_marks, color_to_the_right = colors_used_e14.5ct[4], color_to_the_left = colors_used_e14.5ct[5])
#plotting E14.5_mt markers####
#reverse ifelse in the function(cr_plot_degs) to 1
cr_plot_degs(E14.5_mt_01_marks, color_to_the_right = colors_E14.5_mt[6], color_to_the_left = colors_E14.5_mt[4])
cr_plot_degs(E14.5_mt_02_marks, color_to_the_right = colors_E14.5_mt[6], color_to_the_left = colors_E14.5_mt[1])
cr_plot_degs(E14.5_mt_13_marks, color_to_the_right = colors_E14.5_mt[4], color_to_the_left = colors_E14.5_mt[2])
cr_plot_degs(E14.5_mt_25_marks, color_to_the_right = colors_E14.5_mt[1], color_to_the_left = colors_E14.5_mt[5])
#plotting E17.5_ct markers####
cr_plot_degs(E17.5_ct_01_marks, color_to_the_left = "#7CAE00", color_to_the_right = "#F8766D")
cr_plot_degs(E17.5_ct_12_marks, color_to_the_left = "#00BFC4", color_to_the_right = "#7CAE00")
cr_plot_degs(E17.5_ct_13_marks, color_to_the_left = "#C77CFF", color_to_the_right = "#7CAE00")
cr_plot_degs(E17.5_ct_23_marks, color_to_the_left = "#C77CFF", color_to_the_right = "#00BFC4")
#plotting E17.5_mt markers####
#plotting
cr_plot_degs(E17.5_mt_0_1, color_to_the_left = "#E38900", color_to_the_right =  "#F8766D")
cr_plot_degs(E17.5_mt_1_2, color_to_the_left =  "#C49A00", color_to_the_right = "#E38900")
cr_plot_degs(E17.5_mt_1_3, color_to_the_left = "#99A800", color_to_the_right = "#E38900")
cr_plot_degs(E17.5_mt_3_4, color_to_the_left = "#53B400", color_to_the_right = "#99A800")
cr_plot_degs(E17.5_mt_3_5, color_to_the_left = "#00BC56", color_to_the_right = "#99A800")
cr_plot_degs(E17.5_mt_3_6, color_to_the_left = "#00C094", color_to_the_right = "#99A800")
#
#
"#00B6EB"
"#00BC56" (5)
"#00BFC4"
"#00C094"  (6)
"#06A4FF"
"#53B400" (4)
"#99A800" (3)
"#A58AFF"
"#C49A00" (2)
"#DF70F8"
"#E38900" (1)
"#F8766D" (0)
"#FB61D7"
"#FF66A8"
######
E17.5_mt_0_23456 <- FindMarkers(E17.5_mt, ident.1 = "0", ident.2 = c("2","3","4","5","6"), group.by = "final") %>% filter(p_val_adj < 0.05)
cr_plot_degs(E17.5_mt_0_23456, color_to_the_left = "#4477aa", color_to_the_right = "#F8766D")
#
E17.5_mt_1_23456 <- FindMarkers(E17.5_mt, ident.1 = "1", ident.2 = c("2","3","4","5","6"), group.by = "final")%>% filter(p_val_adj < 0.05)
cr_plot_degs(E17.5_mt_1_23456, color_to_the_left = "#4477aa", color_to_the_right = "#E38900")
