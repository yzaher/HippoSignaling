dir <- -1
#####
cr_plot_degs <- function(degs_list, color_to_the_left,
                         color_to_the_right,
                         axes_size = 30,
                         gene_name_size = 11,
                         legend_size = 20) {
  # Prep the data
  temp <- degs_list %>%
    mutate(
      pct_diff = pct.1 - pct.2,
      gene = rownames(.),
      expression_direction = case_when(
        avg_log2FC > 0 ~ 1,
        avg_log2FC < 0 ~ -1,
        avg_log2FC == 0 ~ 0
      )
    ) %>%
    filter(p_val < 0.05)
  
  # Compute thresholds for labeling extreme genes
  top_fc <- quantile(temp$avg_log2FC, 0.999)
  lower_fc <- quantile(temp$avg_log2FC, 0.001)
  top_pct_diff <- quantile(temp$pct_diff, 0.999)
  lower_pct_diff <- quantile(temp$pct_diff, 0.001)
  
  # Add gene labels for extreme points
  temp_6 <- temp %>%
    mutate(
      names = case_when(
        avg_log2FC > top_fc ~ gene,
        avg_log2FC < lower_fc ~ gene,
        pct_diff > top_pct_diff ~ gene,
        pct_diff < lower_pct_diff ~ gene,
        TRUE ~ NA_character_
      ),
      colors = ifelse(expression_direction == dir, color_to_the_right, color_to_the_left)
    )
  
  # Plot
  plot <- ggplot(temp_6, aes(x = -pct_diff, y = abs(avg_log2FC), color = colors, size = abs(avg_log2FC))) +
    geom_point(alpha = 0.7) +
    geom_text_repel(aes(label = names), size = gene_name_size, box.padding = 0.8, point.padding = 0.8, force = 1.4, max.overlaps = Inf) +
    geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "black") +
    theme_minimal() +
    labs(
      x = "Difference in % of cells expressing the gene between both states",
      y = "Average log2FC",
      color = "Direction of Expression Change",
      size = "Difference in avg_log2FC"
    ) +
    scale_color_identity() +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(5, "lines"),
      legend.text = element_text(size = 20),# face = "bold"),
      axis.title.x = element_text(size = axes_size),#, face = "bold"),  # bigger x-axis label
      axis.title.y = element_text(size = axes_size),# face = "bold"),  # bigger y-axis label
      axis.text.x = element_text(size = legend_size),  # bigger tick labels
      axis.text.y = element_text(size = legend_size)
    )
  plot <- plot +
    guides(
      size = guide_legend(
        title = "Difference in avg_log2FC",
        title.theme = element_text(size = legend_size),  # legend title size
        label.theme = element_text(size = 20)   # legend labels size
      ),
      color = guide_legend(
        title = "Direction of Expression Change",
        title.theme = element_text(size = legend_size),
        label.theme = element_text(size = 20)
      )
    )
  
  # Save the plot
  plotname <- paste0("pub_output/", deparse(substitute(degs_list)), ".png")
  png(filename = plotname, width = 13, height = 8, units = "in", res = 400)
  print(plot)
  dev.off()
}
cr_plot_degs(cl_0_1_marks, color_to_the_right = colors_used_e14.5ct[5], color_to_the_left = colors_used_e14.5ct[1])
