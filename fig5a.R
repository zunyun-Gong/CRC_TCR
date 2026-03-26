#!/usr/bin/env Rscript
#  Fig5a - MSS CDR3 motifs 


library(readxl)
library(ggplot2)
library(dplyr)

# 数据读取
file_path <- "D:/TCR项目/result paper(1)/result paper/supplyment.table12-13.xlsx"
data <- read_excel(file_path, sheet = "Supplyment Table12", skip = 2)

# 计算 logFDR
data$logFDR <- -log10(pmax(data$FDR_fisher, 1e-300))

# 合并 KRAS 和 NRAS
data$Gene_Group <- ifelse(
  grepl("KRAS|NRAS", data$Uniq_Gene),
  "KRAS & NRAS",
  as.character(data$Uniq_Gene)
)

# 只保留突变数 > 5 的基因
gene_counts <- data %>% group_by(Gene_Group) %>% summarise(count = n())
genes_to_keep <- gene_counts$Gene_Group[gene_counts$count > 5]
data_filtered <- data %>% filter(Gene_Group %in% genes_to_keep)

# 按突变数降序排序基因
genes_sorted <- gene_counts %>% 
  filter(count > 5) %>% 
  arrange(desc(count)) %>% 
  pull(Gene_Group)

# 基因配色方案
gene_colors <- c(
  "KRAS & NRAS" = "#e66642",  # 橙红
  "ZFP36L2"     = "#386cb0" ,  # 绿色
  "APC"         = "#7fc97f",  # 紫色
  "TP53"        = "#beaed4",  # 橙色
  "PIK3CA"      = "#fdc086",  # 蓝色
  "PCBP1"       = "#F7F013",  # 黄色
  "SMAD4"       = "#A39063"   # 浅蓝
)

# 创建绘图数据 - X轴为排序index，基因间空20个
current_index <- 1
plot_data <- data.frame()
gene_boundaries <- list()
gene_labels_pos <- list()

for (i in seq_along(genes_sorted)) {
  gene <- genes_sorted[i]
  gene_data <- data_filtered %>%
    filter(Gene_Group == gene) %>%
    arrange(desc(overlap_sp_num))  # 按共出现次数排序
  
  n_mutations <- nrow(gene_data)
  
  # X轴直接为排序的index
  gene_data$x <- seq(current_index, current_index + n_mutations - 1)
  plot_data <- bind_rows(plot_data, gene_data)
  
  # 记录基因边界（用于背景色带）
  gene_boundaries[[length(gene_boundaries) + 1]] <- data.frame(
    gene = gene,
    start = current_index - 0.5,
    end = current_index + n_mutations - 0.5,
    idx = i
  )
  
  # 记录基因标签位置（中间）
  gene_labels_pos[[gene]] <- current_index + (n_mutations - 1) / 2
  
  # 当前基因结束后，加20个空位（最后一个基因不加）
  if (i < length(genes_sorted)) {
    current_index <- current_index + n_mutations + 20
  } else {
    current_index <- current_index + n_mutations
  }
}

gene_boundaries_df <- bind_rows(gene_boundaries)

# 将基因颜色添加到plot_data
plot_data$gene_color <- gene_colors[plot_data$Gene_Group]

# 将 total_positve 转换为因子（离散值）用于形状映射
plot_data$total_positve_factor <- as.factor(plot_data$total_positve)

# 定义形状映射
shape_mapping <- c(
  "1" = 1,   # 空心圈
  "2" = 10,  # 空心圈加十字
  "3" = 19   # 实心圈
)

# 计算基因分隔线位置（每个基因结束后10个位置）
gene_dividers <- data.frame(
  xintercept = gene_boundaries_df$end[-nrow(gene_boundaries_df)] + 10
)

# 绘图
p <- ggplot(plot_data, aes(x = x, y = logFDR)) +
  # 添加基因间黑色虚线分隔
  geom_vline(
    data = gene_dividers,
    aes(xintercept = xintercept),
    linetype = "dashed",
    color = "black",
    linewidth = 0.8
  ) +
  geom_point(
    aes(color = Gene_Group, shape = total_positve_factor, size = overlap_sp_num),
    stroke = 0.6
  ) +
  scale_size_continuous(
    name = "Co-occurence count",
    range = c(3.5, 8)
  ) +
  scale_color_manual(
    values = gene_colors,
    name = "Gene",
    breaks = genes_sorted  # 按突变数排序显示图注
  ) +
  scale_shape_manual(
    values = shape_mapping,
    name = "Tools support TCR-pMHC binding",
    labels = c("1", "2", "3")
  ) +
  guides(
    shape = guide_legend(order = 1, override.aes = list(size = 5)),
    size = guide_legend(order = 2, override.aes = list(size = 4)),
    color = guide_legend(order = 3, override.aes = list(size = 4))
  ) +
  # X轴为数字，每50显示一个刻度
  scale_x_continuous(
    breaks = seq(50, 400, by = 50)
  ) +
  labs(
    x = "Mutation Index",
    y = "-log10(FDR)",
    title = "MSS"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    legend.position = "right"
  ) +
  coord_cartesian(clip = "off", ylim = c(1.2, 3.0))
p
# 保存图片
output_file <- "/root/paper_code/Fig5a.pdf"
ggsave(output_file, p, width = 16, height = 7, dpi = 300, bg = "white")
