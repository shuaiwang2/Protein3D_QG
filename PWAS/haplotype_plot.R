library(ggplot2)
library(tidyr)
library(tibble)  # for rownames_to_column()

# str(PPC.obs)
# str(y.obs)

dir.create("haplotype", showWarnings = FALSE)



# PWAS-HQK.R
PPC <- PPC_list[[sim_metric]]
PPC <- PPC[lengths(PPC) > 0]

y <- pheno[, trait]
obs <- is.finite(y)
N <- sum(obs)
y.obs <- y[obs]
Z.obs <- proteome[, obs]

k <-0
PPC_pangene <- permute(PPC[[pangene]], k)

PPC.obs <- match_to_matrix(Z.obs[pangene, ], PPC_pangene)


# 
y.df <- data.frame(ID = names(y.obs), Trait = as.numeric(y.obs))

#
PPC.df <- as.data.frame(PPC.obs)
PPC.df$ID <- rownames(PPC.df)

df_combin2 <- merge(PPC.df, y.df, by = "ID")




df_combin3 <- df_combin2[,c(PC,"Trait")]

#df_combin3$PC16_factor <- factor(sprintf("%.8f", df_combin3$PC16),ordered = TRUE)

#df_combin3$PC16_factor <- factor(df_combin3$PC16,ordered = TRUE)

            
ggplot(df_combin3, aes(x = PC3, y = Trait)) +
  #geom_boxplot()+
  geom_point()+
  geom_smooth(method = lm)+
  xlab("PC values") + ylab("Phenotypic values") +
  theme_bw(base_size=20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(filename = paste0("haplotype/","point","/",pangene,"_",PC, ".pdf"), width = 6, height = 6)


PPC_pangene1 <- as.data.frame(PPC_pangene)
PPC_pangene1$ID <- rownames(PPC_pangene1)

#nam id
nam_information<- read.table("nam_information.txt",sep='\t')
PPC_pangene1$V2 <- substr(rownames(PPC_pangene1), 1, 8)
pcc_id <-merge(nam_information, PPC_pangene1, by ="V2")
mat <- pcc_id[,c("V1",PC)]

df <- melt(mat)

ggplot(df, aes(x = V1, y = variable, fill = value)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(value,6)), color = "black", angle = 90) +
  scale_fill_gradient(low = "white", high = "white") +
  theme_void(base_size=14) +
  theme(axis.title = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))

ggsave(filename = paste0("haplotype/","point","/",pangene,"heatmap_",PC, ".pdf"), width = 12, height = 3)




