setwd("C:\\Users\\a\\Desktop\\Protein_QG\\282")

library(tidyverse)

df <- read.table("combin_pi.txt", header = FALSE, stringsAsFactors = FALSE)

df <- df %>%
  mutate(V1 = gsub("[/_]", "-", V1)) %>%  #  `/` or  `_` to `-`
  separate(V1, into = c("dir", "bin", "plddt","rsa"), sep = "-")

df_summary <- df %>%
  group_by(plddt,rsa) %>%
  summarise(weighted_pi = sum(V2 * V3) / sum(V3), .groups = "drop")

df_summary$plddt <- gsub("plddt1", "[0,50)", df_summary$plddt)
df_summary$plddt <- gsub("plddt2", "[50,70)", df_summary$plddt)
df_summary$plddt <- gsub("plddt3", "[70,90)", df_summary$plddt)
df_summary$plddt <- gsub("plddt4", "[90,100]", df_summary$plddt)


df_summary$rsa_num <- sapply(df_summary$rsa, function(x) {
  n <- as.numeric(sub("rsa", "", x))
  mid_point <- (n - 1) * 0.1 + 0.05
  return(mid_point)
})

ggplot(df_summary,aes(rsa_num,weighted_pi))+
  geom_point(size = 0.8,)+
  geom_line(col="black")+
  theme_bw(base_size = 10)+
  scale_fill_manual(values=c("#66A61E", "#7570B3"))+
  facet_grid(~ plddt, scales = "free_x")+

  scale_x_continuous(
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1),     
    expand = c(0, 0)       
  )+
  labs(x = "RSA", y =  expression(Pi[N]/Pi[S]))+
  ylim(0,1)+
  theme(legend.title = element_blank(),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
ggsave("pi2.pdf",width=6,height=2.5,units="in",dpi=1000)

