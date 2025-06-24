library(ggplot2)
library(pheatmap)

setwd("/home/song/wangs/pwas_pwp")
PRM_file <- "QC/QC.rds"

PRM <- readRDS("PWP/PRM_esm3b.rds")
#
PRM_list <- readRDS("PWP/PRM.rds")

cor(c(PRM[["usalignMatrix_TMScore"]]), c(PRM[["esm3b_SimilarityMatrix"]]))

cor(c(PRM[["usalignMatrix_TMScore"]]), c(PRM[["mafftSequenceSimilarityMatrix"]]))
cor(c(PRM[["usalignMatrix_TMScore"]]), c(PRM[["muscleSequenceSimilarityMatrix"]]))
cor(c(PRM[["usalignMatrix_TMScore"]]), c(PRM[["tcoffeeSequenceSimilarityMatrix"]]))
cor(c(PRM[["mafftSequenceSimilarityMatrix"]]), c(PRM[["tcoffeeSequenceSimilarityMatrix"]]))
cor(c(PRM[["mafftSequenceSimilarityMatrix"]]), c(PRM[["muscleSequenceSimilarityMatrix"]]))
cor(c(PRM[["usalignMatrix_TMScore"]]), c(PRM[["haplotype"]]))

matrix_list <- list(
  Structure = PRM[["usalignMatrix_TMScore"]],
  "Sequence (MAFFT)" = PRM[["mafftSequenceSimilarityMatrix"]],
  "Sequence (Muscle)" = PRM[["muscleSequenceSimilarityMatrix"]],
  "Sequence (T-Coffee)" = PRM[["tcoffeeSequenceSimilarityMatrix"]],
  Haplotype = PRM[["haplotype"]]
)

# 计算相关性矩阵
cor_results <- matrix(nrow = length(matrix_list), ncol = length(matrix_list))
rownames(cor_results) <- colnames(cor_results) <- names(matrix_list)

# 填充相关性矩阵
for (i in 1:length(matrix_list)) {
  for (j in 1:length(matrix_list)) {
    cor_results[i, j] <- cor(c(matrix_list[[i]]), c(matrix_list[[j]]), use = "complete.obs")
  }
}
cor_results <- round(cor_results, 3)

#cor_results[upper.tri(cor_results)] <- NA

write.csv(cor_results,"correlation_heatmap.csv")
# 使用 pheatmap 绘制热图

pheatmap(cor_results,
         display_numbers = cor_results,                # 显示相关系数
         fontsize_number = 15,                  # 数字字体大小
         number_color = "black",
         color = colorRampPalette(c("#D6E6F2","#1868B2"))(500), # 定义颜色
         cluster_rows = FALSE,                    # 显示行谱系图
         cluster_cols = FALSE,                     # 显示列谱系图
         main = "Correlation of Different Similarity Matrices",
         filename = "correlation_heatmap.pdf")  # 标题

#用这个代码

PRM_list <- readRDS("PWP/PRM.rds")

sim_metric <- "haplotype"
PRM_list <- PRM_list[names(PRM_list) != sim_metric]

n <- sapply(PRM_list, nrow) %>% unique()
H <- diag(n) - matrix(1, nrow=n, ncol=n)/n

X1 <- sapply(PRM_list, function(mat) {
  return(mat[upper.tri(mat, diag=TRUE)])
})

X2 <- sapply(PRM_list, function(mat) {
  out <- H %*% tcrossprod(mat, H)
  return(out[upper.tri(out, diag=TRUE)])
})

#R1 <- cor(X1)
R2 <- cor(X2)

R2 <- round(R2, 3)

pheatmap(R2,
         display_numbers = R2,                
         fontsize_number = 15,# 数字字体大小
         angle_col = 45,
         number_color = "black",
         color = colorRampPalette(c("#D6E6F2","#1868B2"))(500), # 定义颜色
         cluster_rows = FALSE,                    # 显示行谱系图
         cluster_cols = FALSE,                     # 显示列谱系图
         #main = "Correlation of Different Similarity Matrices",
         filename = "correlation_centering.pdf")



sim_metrics <- c(
  "usalignMatrix_TMScore",
  "mafftSequenceSimilarityMatrix",
  "tcoffeeSequenceSimilarityMatrix",
  "muscleSequenceSimilarityMatrix",
  "haplotype"
) 
PPC_list <- readRDS("QC/PPC_list.rds")

# 5133 lg1 not su1
for (i in sim_metrics) {
  su1 <- as.data.frame(PPC_list[[i]][["5133"]])
  su1_sorted <- su1[order(su1[, "PC1"]), ]
  cor_matrix_melt <- melt(as.matrix(su1_sorted))
  
  ggplot(cor_matrix_melt, aes(Var2, Var1, fill = value)) +
    geom_tile() +
    geom_text(aes(label = round(value, 3)), color = "black") + # 显示相关系数
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) + # 颜色渐变
    theme_minimal() +
    theme(axis.text.x = element_text( hjust = 1)) +
    labs(title = "", x = "", y = "")
  ggsave(paste0("su1/",i,"5133.pdf"),width = 14,height = 10)
  
}

#su1 13166
i <- "mafftSequenceSimilarityMatrix"
su1 <- as.data.frame(PPC_list[[i]][["13166"]])
su1_sorted <- su1[order(su1[, "PC1"]), ]

cor_matrix_melt <- melt(as.matrix(su1_sorted))

ggplot(cor_matrix_melt, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 3)), color = "black") + # 显示相关系数
  scale_fill_gradient2(low = "#41B6C4", high = "#075930", mid = "white", midpoint = 0) + # 颜色渐变
  theme_minimal() +
  theme(axis.text.x = element_text( hjust = 1)) +
  labs(title = "", x = "", y = "")
ggsave(paste0("su1/su1",i,"13166.pdf"),width = 10,height = 15)


"#41B6C4","#7FCDBB","#C7E9B4","#EDF8B1"


#data <- read.table("nam_information.txt",sep='\t')
#su1_sorted$match_key <- substr(rownames(su1_sorted), 1, 8)
#data$match_key <- data[[2]]
#merged <- merge(su1_sorted, data, by = "match_key",, all.x = TRUE, suffixes = c(".f1", ".f2"))


