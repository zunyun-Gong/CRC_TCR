
library(ggplot2)
library(reshape2)
library(dplyr)
library(readxl)

setwd("D:/TCR项目/result5/")

cat("=== 步骤1: 读取motif丰度数据 ===\n")
motif_data <- read.table("MSSall.stat.matrix.sigmotif.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
cat("原始数据维度:", nrow(motif_data), "motifs x", ncol(motif_data), "samples\n")

# 转换为01矩阵
binary_matrix <- as.matrix((motif_data > 0) * 1)


cat("\n=== 步骤2: 读取样本分类信息 ===\n")
sample_info <- read_excel("supplyment.table1_202411298.xlsx", sheet = 1, skip = 2)
mss_samples <- sample_info %>% 
  filter(grepl("^MSS_[123]$", cluster)) %>%
  select(Sample_TCR, cluster)

binary_matrix_new = binary_matrix[,mss_samples$Sample_TCR]

motif_freq <- rowSums(binary_matrix_new) / ncol(binary_matrix_new)

# 筛选在至少5%样本中出现的motif
selected_motifs <- read.csv("sel_motif.txt",header = F)$V1
binary_filtered <- binary_matrix_new[selected_motifs, ]


# 背景数量
background_counts <- table(mss_samples$cluster)


# 匹配样本
common_samples <- intersect(colnames(binary_filtered), mss_samples$Sample_TCR)
binary_mss <- binary_filtered[, common_samples]
sample_clusters <- mss_samples$cluster[match(common_samples, mss_samples$Sample_TCR)]
names(sample_clusters) <- common_samples

chi_square_results <- data.frame(
  motif = selected_motifs,
  chisq_stat = NA,
  p_value = NA,
  stringsAsFactors = FALSE
)

for (i in seq_along(selected_motifs)) {
  motif <- selected_motifs[i]
  motif_presence <- binary_mss[motif, ]
  
  # 构建列联表
  contingency_table <- table(motif_presence, sample_clusters)
  print(contingency_table[2,])
  # chisq_test <- chisq.test(contingency_table[2,],p= background_counts/sum(background_counts))
  chisq_test <- chisq.test(contingency_table)
  chi_square_results$chisq_stat[i] <- chisq_test$statistic
  chi_square_results$p_value[i] <- chisq_test$p.value
}

# 使用原始p值标注显著性
chi_square_results$significance <- ""
chi_square_results$significance[chi_square_results$p_value < 0.05] <- "*"
chi_square_results$significance[chi_square_results$p_value < 0.01] <- "**"
chi_square_results$significance[chi_square_results$p_value < 0.001] <- "***"


# 显示显著motif
sig_motifs <- chi_square_results[chi_square_results$p_value < 0.05, c("motif", "p_value", "significance")]
sig_motifs <- sig_motifs[order(sig_motifs$p_value), ]




plot_data <- data.frame()

for (motif in selected_motifs) {
  motif_presence <- binary_mss[motif, ]
  
  for (cluster in c("MSS_1", "MSS_2", "MSS_3")) {
    cluster_samples <- names(sample_clusters)[sample_clusters == cluster]
    freq <- sum(motif_presence[cluster_samples]) / length(cluster_samples)
    
    plot_data <- rbind(plot_data, data.frame(
      motif = motif,
      cluster = cluster,
      frequency = freq,
      stringsAsFactors = FALSE
    ))
  }
}

# 添加显著性标记
plot_data$significance <- chi_square_results$significance[match(plot_data$motif, chi_square_results$motif)]
plot_data$cluster <- factor(plot_data$cluster, levels = c("MSS_1", "MSS_2", "MSS_3"))

# # 按MSS_1频率排序
motif_order <- plot_data %>%
  filter(cluster == "MSS_1") %>%
  arrange(desc(frequency)) %>%
  pull(motif)

plot_data$motif <- factor(plot_data$motif, levels = motif_order)

cat("\n=== 步骤5: 绘制堆积柱状图（使用原始p值）===\n")
options(bitmapType = "cairo")

# 下面的热图
cat("=== 读取sel_motif.txt ===\n")
sel_motifs <- readLines("sel_motif.txt")
sel_motifs <- trimws(sel_motifs)
sel_motifs <- sel_motifs[sel_motifs != ""]
cat("选择的motif数量:", length(sel_motifs), "\n")

cat("\n=== 读取数据 ===\n")
motif_data <- read.table("MSSall.stat.matrix.sigmotif.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
existing_motifs <- intersect(sel_motifs, rownames(motif_data))
cat("数据中存在:", length(existing_motifs), "个motif\n")
motif_subset <- motif_data[existing_motifs, ]

sample_info <- read_excel("supplyment.table1_202411298.xlsx", sheet = 1, skip = 2)
mss_samples <- sample_info %>% 
  filter(grepl("^MSS_[123]$", cluster)) %>%
  select(Sample_TCR, cluster)

compare_data <- read.table("total.compareMotifabundence.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

common_samples <- intersect(colnames(motif_subset), mss_samples$Sample_TCR)
motif_mss <- motif_subset[, common_samples]
sample_clusters <- mss_samples$cluster[match(common_samples, mss_samples$Sample_TCR)]
names(sample_clusters) <- common_samples

cat("\n=== 计算各MSS亚型的均值 ===\n")
mean_matrix <- matrix(NA, nrow = length(existing_motifs), ncol = 3)
rownames(mean_matrix) <- existing_motifs  # 保持原始顺序
colnames(mean_matrix) <- c("MSS_1", "MSS_2", "MSS_3")

for (motif in existing_motifs) {
  for (cluster in c("MSS_1", "MSS_2", "MSS_3")) {
    cluster_samples <- names(sample_clusters)[sample_clusters == cluster]
    mean_matrix[motif, cluster] <- mean(as.numeric(motif_mss[motif, cluster_samples]), na.rm = TRUE)
  }
}

cat("\n=== 按行标准化 (Z-score) ===\n")
# 按行标准化: (x - mean) / sd
mean_matrix_scaled <- t(apply(mean_matrix, 1, function(x) {
  if (sd(x) == 0) return(c(0, 0, 0))  # 避免除0
  (x - mean(x)) / sd(x)
}))

cat("标准化后范围:", min(mean_matrix_scaled), "to", max(mean_matrix_scaled), "\n")

cat("\n=== 准备显著性标记 ===\n")
sig_matrix <- matrix("", nrow = length(existing_motifs), ncol = 3)
rownames(sig_matrix) <- existing_motifs
colnames(sig_matrix) <- c("MSS_1", "MSS_2", "MSS_3")

for (motif in existing_motifs) {
  if (motif %in% rownames(compare_data)) {
    # MSS_1 vs others
    p1 <- compare_data[motif, "MSS_1_vs_others.pvalue"]
    if (!is.na(p1)) {
      if (p1 < 0.001) sig_matrix[motif, "MSS_1"] <- "***"
      else if (p1 < 0.01) sig_matrix[motif, "MSS_1"] <- "**"
      else if (p1 < 0.05) sig_matrix[motif, "MSS_1"] <- "*"
    }
    
    # MSS_2 vs others
    p2 <- compare_data[motif, "others_vs_MSS_2.pvalue"]
    if (!is.na(p2)) {
      if (p2 < 0.001) sig_matrix[motif, "MSS_2"] <- "***"
      else if (p2 < 0.01) sig_matrix[motif, "MSS_2"] <- "**"
      else if (p2 < 0.05) sig_matrix[motif, "MSS_2"] <- "*"
    }
    
    # MSS_3 vs others
    p3 <- compare_data[motif, "others_vs_MSS_3.pvalue"]
    if (!is.na(p3)) {
      if (p3 < 0.001) sig_matrix[motif, "MSS_3"] <- "***"
      else if (p3 < 0.01) sig_matrix[motif, "MSS_3"] <- "**"
      else if (p3 < 0.05) sig_matrix[motif, "MSS_3"] <- "*"
    }
  }
}


heat_df <- melt(mean_matrix_scaled)
colnames(heat_df) <- c("Motif", "Group", "Zscore")


sig_df <- melt(sig_matrix)
colnames(sig_df) <- c("Motif", "Group", "Significance")

heat_df$Significance <- sig_df$Significance


heat_df$Motif <- factor(heat_df$Motif, levels = rownames(mean_matrix_scaled))
heat_df$Group <- factor(heat_df$Group, levels = colnames(mean_matrix_scaled))

######## 绘图 ##########

##################### 拼接图片 ##################


colors <- c(MSS_1 = "#52B450", MSS_2 = "#EFC922",MSS_3="#508CCA")

# 准备显著性标记数据
sig_data <- plot_data %>%
  filter(significance != "") %>%
  group_by(motif) %>%
  summarise(
    max_freq = sum(frequency),
    sig = first(significance),
    .groups = "drop"
  )

plot_data <- plot_data[order(plot_data[[1]]), ]

# 再把第一列转成按当前顺序排列的因子
plot_data[[1]] <- factor(plot_data[[1]], levels = unique(plot_data[[1]]))

p1 <- ggplot(plot_data, aes(x = motif, y = frequency, fill = cluster)) +
  # --- 修改 1：增加 width 到 0.9 (箱体变更宽) ---
  geom_bar(stat = "identity", position = "stack", width = 0.9) +
  
  scale_fill_manual(values = colors, name = "MSS Cluster") +
  
  # 保持左侧留白（控制离Y轴的距离）
  scale_x_discrete(expand = expansion(add = c(1.2, 0.6))) +
  
  # 保持Y轴贴底设置
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
  
  labs(x = "Motif", y = "Frequency") + # 稍微简化了标题以便展示
  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.y = element_text(size = 10), 
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

# 添加显著性标记
if (nrow(sig_data) > 0) {
  p1 <- p1 + geom_text(
    data = sig_data,
    # --- 修改 2：极度减小偏移量，让星号下降 ---
    # 将原来的 + 0.01 改为 + 0.002，几乎贴在柱子上
    aes(x = motif, y = max_freq + 0.002, label = sig),
    inherit.aes = FALSE,
    size = 4,
    color = "black",
    fontface = "bold",
    vjust = 0 
  )
}
p1


heat_df[[1]] = factor(heat_df[[1]], levels = unique(plot_data[[1]]))

p2 <- ggplot(heat_df, aes(y = Group, x = Motif, fill = Zscore)) +
  geom_tile(color = "white", linewidth = 0.3)+
  geom_text(aes(label = Significance), 
            color = "black", 
            size = 3) +
  scale_fill_gradient2(
    low = "#15ACFD",    # 浅蓝色
    mid = "white",        # 中间是白色
    high = "#FD2715",      # 浅红色
    midpoint = 0,
    limits = c(-1, 1),    # 将色阶压缩到[-1, 1]
    oob = scales::squish,
    name = "Z-score"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),  # x轴字体倾斜
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    axis.title = element_blank()
  )
p2

p1 <- p1 +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

library(patchwork)
combined_plot <- p1 / p2 + plot_layout(heights = c(1.2, 0.3))  # 控制上下比例
ggsave("fig5b.pdf", plot = combined_plot, width = 8, height = 6)
