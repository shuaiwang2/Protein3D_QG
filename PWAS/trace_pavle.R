#

QC_file <- "QC/QC.rds"
PWAS_file <- "PWAS/PWAS.rds"
fig_dir <- "figs/PWAS/HQK/"

QC <- readRDS(QC_file)
pwas <- readRDS(PWAS_file)

hit_max <- read.csv("PWAS/QC_hits_max_B73",sep=' ')
BC_threshold <- 4.744733e-06


cor_results <- data.frame(trait = character(0), correlation = numeric(0))


# 循环遍历每个 trait
for (i in unique(pwas$trait)) {
  
  pwas_filtered <- pwas %>%
    filter(trait == i)
  
  colnames(pwas_filtered)[colnames(pwas_filtered) == "sim_metric"] <- "metric"
  
  merged_data <- merge(pwas_filtered, QC[, c("pangene", "trace")], by = "pangene")
  
  correlation_result <- merged_data %>%
    group_by(metric) %>%
    summarise(correlation = cor(pval, trace, use = "complete.obs"))
  
  cor_results <- rbind(cor_results, data.frame(trait = i,metric = correlation_result$metric
                                               , correlation = correlation_result$correlation))
  
}




# pwas-interpretation.R 过滤一下
hits_max2 <- hits_max[hits_max$mafftSequence > -log10(BC_threshold) | hits_max$Structure > -log10(BC_threshold), ]

hits_max2 <- merge(hits_max, QC[, c("pangene", "trace")], by = "pangene")

hits_max2 <- hits_max2[c(1,2,3,6,7)]

# 第一二三四象限
QC_hits_max_B73_1 <- hits_max2[hits_max2$mafftSequence  > -log10(BC_threshold) & hits_max2$Structure > -log10(BC_threshold), ]
QC_hits_max_B73_2 <- hits_max2[hits_max2$mafftSequence  < -log10(BC_threshold) & hits_max2$Structure > -log10(BC_threshold), ]
QC_hits_max_B73_3 <- hits_max2[hits_max2$mafftSequence  < -log10(BC_threshold) & hits_max2$Structure < -log10(BC_threshold), ]
QC_hits_max_B73_4 <- hits_max2[hits_max2$mafftSequence  > -log10(BC_threshold) & hits_max2$Structure < -log10(BC_threshold), ]



QC_hits_max_B73_1 <- QC_hits_max_B73_1 %>%
    pivot_longer(cols = matches("^(mafftSequence|muscleSequence|tcoffeeSequence|Structure)"),  # 正则表达式匹配多个列名
                 names_to = "Metric",values_to = "pvaule") %>%
                  mutate(file_name = "QC_hits_max_B73_1")

QC_hits_max_B73_2 <- QC_hits_max_B73_2 %>%
  pivot_longer(cols = matches("^(mafftSequence|muscleSequence|tcoffeeSequence|Structure)"),  # 正则表达式匹配多个列名
               names_to = "Metric",values_to = "pvaule")%>%
                mutate(file_name = "QC_hits_max_B73_2")

QC_hits_max_B73_3 <- QC_hits_max_B73_3 %>%
  pivot_longer(cols = matches("^(mafftSequence|muscleSequence|tcoffeeSequence|Structure)"),  # 正则表达式匹配多个列名
               names_to = "Metric",values_to = "pvaule")%>%
                mutate(file_name = "QC_hits_max_B73_3")

QC_hits_max_B73_4 <- QC_hits_max_B73_4 %>%
  pivot_longer(cols = matches("^(mafftSequence|muscleSequence|tcoffeeSequence|Structure)"),  # 正则表达式匹配多个列名
               names_to = "Metric",values_to = "pvaule")%>%
              mutate(file_name = "QC_hits_max_B73_4")


QC_hits_max_B73 <- rbind(QC_hits_max_B73_1,QC_hits_max_B73_2,QC_hits_max_B73_3,QC_hits_max_B73_4)



QC_hits_max_B73$Metric[QC_hits_max_B73$Metric == "mafftSequence"] <- "Sequence (MAFFT)"
QC_hits_max_B73$Metric[QC_hits_max_B73$Metric == "muscleSequence"] <- "Sequence (Muscle)"
QC_hits_max_B73$Metric[QC_hits_max_B73$Metric == "tcoffeeSequence"] <- "Sequence (T-Coffee)"

QC_hits_max_B73$file_name[QC_hits_max_B73$file_name == "QC_hits_max_B73_1"] <- "First quadrant"
QC_hits_max_B73$file_name[QC_hits_max_B73$file_name == "QC_hits_max_B73_2"] <- "Second quadrant"
QC_hits_max_B73$file_name[QC_hits_max_B73$file_name == "QC_hits_max_B73_4"] <- "Fourth quadrant"

QC_hits_max_B73_cor <- QC_hits_max_B73 %>%
  group_by(file_name, Metric) %>%
  mutate(Correlation = cor(pvaule, trace)) %>%
  ungroup() %>%
  select(file_name, Metric, Correlation) %>%
  distinct()

QC_hits_max_B73_cor$file_name <- factor(QC_hits_max_B73_cor$file_name, 
                                        levels = c("First quadrant", 
                                                   "Second quadrant", 
                                                   "Fourth quadrant"))

#不是一个数据
ggplot(QC_hits_max_B73_cor,aes(x=Metric, y=Correlation)) +
  geom_point() +
  facet_wrap(~ file_name, scales = "free_y") +
  theme_bw(base_size=12) +
  labs(x="Metric", y="Trace")+
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) 
ggsave(filename = paste0(fig_dir, "QC_hits_max_B73_cor.pdf"),width = 12, height = 5)

ggplot(QC_hits_max_B73,aes(x=Metric, y=trace,fill=Metric)) +
  geom_boxplot() +
  facet_wrap(~ file_name, scales = "free_y") +
    theme_bw(base_size=15) +
  scale_fill_manual(values = c("#3B549D","#33A2FF"))+
  labs(x="Metric", y="Trace")+
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) 
ggsave(filename = paste0(fig_dir, "QC_hits_max_B73_trace.pdf"),width = 10, height = 5)



# 第一次 整体画
cor_results$trait=sub("_.+", "", cor_results$trait)

cor_results$trait[cor_results$trait == "DTA"] <- "Days to anthesis"
cor_results$trait[cor_results$trait == "DTS"] <- "Days to silking"
cor_results$trait[cor_results$trait == "aminoacids"] <- "Amino acids"
cor_results$trait[cor_results$trait == "asi"] <- "Anthesis-silking interval"
cor_results$trait[cor_results$trait == "chlorophylla"] <- "Chlorophyll a"
cor_results$trait[cor_results$trait == "chlorophyllb"] <- "Chlorophyll b"
cor_results$trait[cor_results$trait == "cobdiam"] <- "Cob diameter"
cor_results$trait[cor_results$trait == "coblength"] <- "Cob length"
cor_results$trait[cor_results$trait == "earmass"] <- "Ear mass"
cor_results$trait[cor_results$trait == "earrowno"] <- "Ear row number"
cor_results$trait[cor_results$trait == "fructose"] <- "Fructose"
cor_results$trait[cor_results$trait == "fumarate"] <- "Fumarate"
cor_results$trait[cor_results$trait == "glucose"] <- "Glucose"
cor_results$trait[cor_results$trait == "glutamate"] <- "Glutamate"
cor_results$trait[cor_results$trait == "kernelnoperrow"] <- "kernel number per row"
cor_results$trait[cor_results$trait == "leaflength"] <- "Leaf length"
cor_results$trait[cor_results$trait == "leafwidth"] <- "Leaf width"
cor_results$trait[cor_results$trait == "malate"] <- "Malate"
cor_results$trait[cor_results$trait == "nitrate"] <- "Nitrate"
cor_results$trait[cor_results$trait == "nodenumberaboveear"] <- "Node number above ear"
cor_results$trait[cor_results$trait == "nodenobelowear"] <- "Node number below ear"
cor_results$trait[cor_results$trait == "numbraceroots"] <- "Number brace roots"
cor_results$trait[cor_results$trait == "plantheight"] <- "Plant height"
cor_results$trait[cor_results$trait == "protein"] <- "Protein"
cor_results$trait[cor_results$trait == "starch"] <- "Starch"
cor_results$trait[cor_results$trait == "sucrose"] <- "Sucrose"
cor_results$trait[cor_results$trait == "tasslength"] <- "Tassel length"
cor_results$trait[cor_results$trait == "tassprimbranchno"] <- "Tassel primary branch number"
cor_results$trait[cor_results$trait == "totalkernelno"] <- "Total kernel number"
cor_results$trait[cor_results$trait == "totalkernelweight"] <- "Total kernel weight"
cor_results$trait[cor_results$trait == "upperleafangle"] <- "Upper leaf angle"
cor_results$trait[cor_results$trait == "weight20kernels"] <- "Weight 20 kernels"


ggplot(cor_results,aes(x=trait,y=correlation))+geom_point()+
  theme_bw(base_size=15) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
  labs(x="Trait", y="Correlation") 
ggsave(filename = paste0(fig_dir, "cor_trace_PC", ".pdf"),width = 10, height = 5)


#不同pangene之间的PC
PPC_file <- "QC/PPC_list.rds"
PPC_list <- readRDS(PPC_file)
sim_metric <- "haplotype"
PPC_list <- PPC_list[names(PPC_list) != sim_metric]



result <- data.frame(
  BigList = character(0),   # 第一列：大list的名称
  SmallList = character(0), # 第二列：小list的名称
  ColumnCount = integer(0)  # 第三列：列数
)


for (big_list_name in names(PPC_list)) {

  small_lists <- PPC_list[[big_list_name]]

  for (small_list_name in names(small_lists)) {
    matrix_data <- small_lists[[small_list_name]]
    
    if (is.null(matrix_data) || ncol(matrix_data) == 0) {
      #column_count <- 0
    } else {
      column_count <- ncol(matrix_data)
    }
    
    result <- rbind(result, data.frame(BigList = big_list_name, 
                                       SmallList = small_list_name, 
                                       ColumnCount = column_count))
  }
}

max(result[result$BigList=="usalignMatrix_TMScore",]$ColumnCount) #25 
min(result[result$BigList=="usalignMatrix_TMScore",]$ColumnCount) #1
mean(result[result$BigList=="usalignMatrix_TMScore",]$ColumnCount) #7.8

max(result[result$BigList!="usalignMatrix_TMScore",]$ColumnCount) #25 
min(result[result$BigList!="usalignMatrix_TMScore",]$ColumnCount) #1
mean(result[result$BigList!="usalignMatrix_TMScore",]$ColumnCount) #7.8

max(result[result$BigList=="mafftSequenceSimilarityMatrix",]$ColumnCount) #25 
min(result[result$BigList=="mafftSequenceSimilarityMatrix",]$ColumnCount) #1
mean(result[result$BigList=="mafftSequenceSimilarityMatrix",]$ColumnCount) #6.27


result$BigList[result$BigList == "mafftSequenceSimilarityMatrix"] <- "Sequence (MAFFT)"
result$BigList[result$BigList == "muscleSequenceSimilarityMatrix"] <- "Sequence (Muscle)"
result$BigList[result$BigList == "tcoffeeSequenceSimilarityMatrix"] <- "Sequence (T-Coffee)"
result$BigList[result$BigList == "usalignMatrix_TMScore"] <- "Structure"

ggplot(data=result, aes(x=BigList,y=ColumnCount,fill=BigList))+
  geom_boxplot()+
  scale_fill_manual(values = c("#3B549D", "#D6AEDD","#73ABCF","#33A2FF")) +
  theme_bw(base_size=15)+
  #scale_fill_manual(values=sim_colors)+
  labs(x="Metric", y="PC count") 
ggsave(filename = paste0(fig_dir, "PC_count", ".pdf"),width = 10, height = 4.5)

