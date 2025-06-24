install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)
library(reshape2)


pheno <- read.table("C:\\Users\\a\\Desktop\\Protein_QG\\Table\\Supplemental_Table_S3.csv",sep = ',',header = T,row.names = 1)


pheno <- read.table("C:\\Users\\a\\Desktop\\Protein_QG\\Table\\all_NAM_phenos.csv",sep = ',',header = T)
pheno1 <- pheno[,c(24:59)]

# colnames(pheno) <-c("Days to anthesis","Days to silking","Anthesis-silking interval",
#                     "Plant height", "Leaf length","Leaf width","Tassel primary branch number",
#                     "Tassel length","Upper leaf angle","Cob diameter","Cob length","Ear row number",
#                     "kernel number per row","Ear mass","Total kernel weight","Weight 20 kernels",
#                     "Total kernel number","Node number below ear","Node number above ear",
#                     "Number brace roots","Chlorophyll a","Chlorophyll b","Malate","Fumarate",
#                     "Glutamate","Amino acids","Protein","Nitrate","Starch","Sucrose",
#                     "Glucose","Fructose"
#                    )

pheno_cor <- cor(pheno1,use = "pairwise.complete.obs", method = "pearson")
#cor_long <- melt(pheno_cor)
write.csv(pheno_cor,"pheno_cor.csv")
#chart.Correlation(pheno_cor, histogram=TRUE, pch="+")
