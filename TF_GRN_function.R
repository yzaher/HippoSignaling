make_tf_network_orig <- function(degs, object, adj_table, tf_list, 
                            n_top = 50, avg_fc_cutoff = 1, pct_cutoff = 0) {  #added here pct_cutoff 0.2 -> 0 #oct_25
  
  library(dplyr)
  library(igraph)
  library(ggraph)
  library(purrr)

#  degs$avg_log2FC <- (degs$avg_log2FC * -1)
  
  # --- Filter DEGs ---
  degs_filtered <- degs %>% 
  filter(abs(avg_log2FC) > avg_fc_cutoff) #added here abs



  
  # --- Filter adjacency table ---
  top_adj <- adj_table %>% 
    filter(target %in% degs_filtered$gene,
           TF %in% tf_list)
  
  # --- Select top DEGs ---
  selected_genes <- degs_filtered %>% 
    arrange(desc(avg_log2FC)) %>% 
    filter(abs(pct.1) > pct_cutoff) %>% 
    filter(abs(pct.2) > pct_cutoff) %>% #added here whole line
    slice_head(n = n_top) %>% 
    slice_tail(n = n_top) #added here whole line
  
  selected_gene_names <- selected_genes$gene[!(selected_genes$gene %in% adj_table$TF)]
  
  # --- Mark top DEGs ---
  top_adj <- top_adj %>%
    mutate(Top100 = target %in% selected_gene_names)
  
  # --- Correlation calculation ---
  working_mat <- GetAssayData(object, assay = "RNA", layer = "data")
  
  top_adj$cor <- map_dbl(seq_len(nrow(top_adj)), function(i) {
    df <- top_adj[i, ]
    if (df$TF %in% rownames(working_mat) &&
        df$target %in% rownames(working_mat)) {
      x <- working_mat[df$TF, ]
      y <- working_mat[df$target, ]
      if (sd(x) > 0 && sd(y) > 0) {
        ct <- suppressWarnings(cor.test(x, y, exact = FALSE))
        if (!is.na(ct$p.value) && ct$p.value < 0.05) return(ct$estimate)
      }
    }
    return(NA_real_)
  })
  
  # --- Edge color based on correlation ---
  top_adj <- top_adj %>%
    mutate(
      edge_color = case_when(
        is.na(cor) ~ "gray70",
        cor > 0 ~ "skyblue",
        cor < 0 ~ "red3"
      )
    )
  
  # --- Node info ---
  all_nodes <- unique(c(top_adj$TF, top_adj$target))
  node_info <- data.frame(
    name = all_nodes,
    type = ifelse(all_nodes %in% tf_list, "TF", "DEGs"),
    stringsAsFactors = FALSE
  )
  
  node_info$Colors <- case_when(
    node_info$name %in% selected_gene_names ~ "Within top 100 DEGs", #added here 50 to 100
    node_info$type == "TF" ~ "TF",
    TRUE ~ "DEGs"
  )
  
  # --- Graph ---
  adj_net <- graph_from_data_frame(
    d = top_adj[, c("TF", "target", "edge_color")],
    vertices = node_info,
    directed = FALSE
  )
  
  # --- Plot ---
  p <- ggraph(adj_net, layout = "stress") +   
    geom_edge_link(
      aes(color = edge_color),
      alpha = 0.8, width = 1.4, show.legend = TRUE  #added here width 1.2 to 1.4 #oct_25
    ) +
    geom_node_point(aes(color = Colors), size = 15) + #added here size 6 to 10 #oct_25
    geom_node_text(aes(label = name), repel = TRUE, size = 20, color = "black", face ="bold") + #added here size 3.8 to 5 and face ="bold"  #oct_25
    scale_color_manual(
      name = "Node type",
      values = c(
        "TF" = "tomato",
        "DEGs" = "skyblue",
        "Within top 100 DEGs" = "purple" #added here 50 to 100
      )
    ) +
    scale_edge_color_identity(
      name = "Correlation direction",
      breaks = c("skyblue", "red3", "gray70"),
      labels = c("Positive", "Negative", "Undetermined"),
      guide = "legend"
    ) +
    theme_void(base_size = 14) +
    theme(
      legend.position = "top",
      legend.title = element_text(size = 50), #added here size 13 to 50 #oct_25
      legend.text = element_text(size = 50),                #added here size 12 to 50 #oct_25
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}
