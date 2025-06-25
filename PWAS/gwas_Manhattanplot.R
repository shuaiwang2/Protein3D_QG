library(dplyr)

setwd("/home/song/wangs/snpgwas")
# Manhattan plot
dat <-read.table("./output/GE_GWAS13.assoc.txt",header = T)
data <-read.table("./pwas_gwas_merge.txt",header = T)

#files <- c("./output/GE_GWAS10.assoc.txt", "./output/GE_GWAS13.assoc.txt", "GE_GWAS20.assoc.txt")
# # 读取所有文件，并添加文件名列
# df_list <- lapply(files, function(f) {
#   read.table(f, header = TRUE, sep = "\t") %>%
#     mutate(trait = basename(f))  # 加文件名列（只取文件名，不带路径）
# })
# 
# 
# data <- bind_rows(df_list)
# data$trait <- gsub("GE_GWAS13.assoc.txt", "Ear row number", data$trait)
# data$trait <- gsub("GE_GWAS20.assoc.txt", "Node number above ear", data$trait)
# data$trait <- gsub("GE_GWAS10.assoc.txt", "Upper leaf angle", data$trait)
data <- data %>%
  pivot_longer(
    cols = c(Structure,min_p),      # 要转换的列
    names_to = "variable",         # 新列：原列名变成变量名
    values_to = "min_p"            # 新列：原来的值
  )

colnames(data) <- c("trait","gene","chr","ps","variable","min_p")




threshold <- (0.05/1090989) #max number of number of analyzed SNPs/var for all traits

pwas <-data
pwas$chr <- factor(pwas$chr, levels = paste0("chr", 1:10))

# 1)计算chr长度
chr_len <- pwas %>% 
  group_by(chr) %>% 
  summarise(chr_len=max(ps))
# 2） 计算每条chr的初始位置
chr_pos <- chr_len  %>% 
  mutate(total = cumsum(chr_len) - chr_len) %>%
  select(-chr_len)
#3)计算累计SNP的位置
pwas <- chr_pos  %>%
  left_join(pwas, ., by="chr") %>%
  arrange(chr, ps) %>%
  mutate( BPcum = ps + total)

X_axis <- pwas %>%
  mutate(BPcum = as.numeric(BPcum)) %>%
  group_by(chr) %>%
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

X_axis<-X_axis[c(1:10),]
#pwas$logp <- -log10(pwas$p_wald)
#max_y <- max(pwas$logp, na.rm = TRUE)

sim_colors <- c(
  "Structure"             = "#238b45",  # 高亮深绿
  #"Sequence (Muscle)"     = "#e31a1c",  # 高亮红
  "min_p"   = "#6a3d9a"  # 紫
  #"Sequence (MAFFT)"      = "#1f78b4"   # 蓝
)


pwas <- pwas[(pwas$trait =="upperleafangle" | pwas$trait =="earrowno"| pwas$trait =="nodenumberaboveear"),]

ggplot(pwas, aes(x=BPcum, y=-log10(min_p))) +
  #设置点的大小，透明度
  geom_point( aes(color=variable), size=0.5) +
  #设置颜色
  scale_color_manual(values=sim_colors)+
  #设定X,Y轴
  scale_x_continuous(expand = c(0, 0),label = X_axis$chr, breaks= X_axis$center ) +
  scale_y_continuous(expand = c(0, 0)) +
  #添加阈值线
  geom_hline(yintercept = -log10(BC_threshold), color ='red',linetype="dotted",size = 0.8) + 
  geom_hline(yintercept = 5.32, color ='#e31a1c',linetype="dotted",size = 0.8) +
  #设置主题
  labs(x = "Chromosome", y = expression(-log[10]~"(p)")) +
  facet_wrap(~ trait, ncol = 1, scales = "free_y") + 
  theme_bw() +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18),
    panel.grid=element_blank(),
    panel.border = element_blank(),
    axis.line= element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text = element_text(size = 13))
ggsave("gwas_pwas.pdf",width = 9, height = 8)
ggsave("gwas_pwas.jpeg",width = 9, height = 8)
