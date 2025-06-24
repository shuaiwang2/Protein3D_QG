# Title     : SNP input
# Objective : Weight SNPs by AlphaFold annotation
# Created by: Guillaume Ramstein (ramstein@qgg.au.dk)

#--------------------------------------------------------
# Script parameters
#--------------------------------------------------------
# Working directory
wd <- ifelse(Sys.info()["user"] == "au686885",
             "C:/Users/au686885/OneDrive - Aarhus Universitet/QGG_plants/alphafold_QG/",
             "~/QGG_plants/alphafold_QG/")
dir.create(wd, showWarnings = FALSE)
setwd("C:\\Users\\a\\Desktop\\Protein_QG\\Fig\\Fig2")

# Input
annot_file <- "C:\\Users\\a\\Desktop\\Protein_QG\\Fig\\Fig2\\DSSP.overlap_282.txt"
annot_header <- c("SNP_type", "chr", "position", "allele_count", "transcript",
                  "protein_length", "protein_seq_id", "refAle", "altAle",
                  "outRefAle", "outAle1", "outAle2","outAle3", "outAle4",
                  "codon_index", "refFrq", "altFrq", "refCodon", "refAA",
                  "altCodon", "altAA", "blosum62_score", "cds_index", 
                  "GERP_score", "pLDDT", "ASA", "maxASA")

gds_file <- "data/coding_gene.gds"

# Graphics parameters
base_size <- 20

q2.5 <- qnorm(0.025)
q97.5 <- qnorm(0.975)

# Output
# output_folder <- "data/SNP_effects/"
# dir.create(output_folder, showWarnings=FALSE)
# 
# fig_dir <- "figs/SNP_effects/"
# dir.create(fig_dir, showWarnings=FALSE)

#--------------------------------------------------------
# Functions
#--------------------------------------------------------
library(foreach)
library(data.table)
library(tidyverse)
library(dplyr)
library(patchwork)
library(gridGraphics)

library(mgcv)
library(tidymv)

#--------------------------------------------------------
# Data
#--------------------------------------------------------
# SNP annotations
annot_full <- read.table(annot_file) %>%
  as.data.frame() %>%
  setNames(annot_header) %>%
  mutate(rASA=pmin(1, ASA/maxASA),
         transcript=factor(transcript),
         pLDDT_class=cut(pLDDT, c(0, 50, 70, 90, 100),right = FALSE)) %>%
  dplyr::rename(CHROM=chr,
                POS=position)

annot <- dplyr::filter(annot_full,
                       is.finite(GERP_score) & 
                         GERP_score >= 0 & 
                         SNP_type == "nonsynonymous") %>%
  mutate(transcript=factor(transcript))


dat <- annot %>% 
  mutate(RSA_class = ifelse(rASA >= 0.5, "exposed", "buried"))  #%>%
#mutate(length_class=cut(protein_length,c(0,median(annot$protein_length),max(annot$protein_length))))


t.test(dat$pLDDT[dat$RSA_class == "exposed"], dat$pLDDT[dat$RSA_class == "buried"])

m0 <- lm(GERP_score ~ RSA_class,data=dat)
m1 <- lm(GERP_score ~ RSA_class *pLDDT_class,data=dat)
m2 <- lm(GERP_score ~ RSA_class * pLDDT_class*protein_length,data=dat)
anova(m0,m1)
anova(m1,m2)

m1 <- lm(GERP_score ~ pLDDT_class*protein_length,data=dat)
m2 <- lm(GERP_score ~ pLDDT_class*protein_length+RSA_class,data=dat)
m3 <- lm(GERP_score ~ pLDDT_class*protein_length+RSA_class+ RSA_class:pLDDT_class,data=dat)
anova(m1,m2)
anova(m2,m3)

ggplot(dat,aes(RSA_class,GERP_score, fill = RSA_class))+
  geom_boxplot(alpha=0.9,width = 0.4,size = 0.2, color = "black") +
  scale_fill_manual(values = c("#BEB8DC", "#82B0D2")) +
  #theme_bw(base_size = 15)+
  labs(x = "Residues", y = "GERP")+
  guides(fill=FALSE)+
  theme(panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_blank())

ggsave("exposed-buried.pdf",width = 3, height = 4)



##GERP and RSA by plddt class

ggplot(annot,aes(GERP_score,rASA))+
  #geom_point(size=0.01,alpha=0.1)+
  geom_smooth(linetype="solid", alpha=0.4,color ="black")+
  theme_bw(base_size = 15)+
  facet_grid(~ pLDDT_class)+
  labs(x = "GERP", y = "RSA")+
  geom_label(data=data.frame(x=max(annot$GERP_score),
                             y=0,
                             pLDDT_class=tail(levels(annot$pLDDT_class), 1),
                             label="p-value < 2.2e-16"),
             aes(x=x, y=y, label=label),
             hjust="right",
             size=4)
ggsave("GERP_GERPRSA_plddtclass.pdf",width=8,height=4,units="in",dpi=1000)

m0 <- lm(rASA ~ GERP_score,data=annot)
m1 <- lm(rASA ~ GERP_score *pLDDT_class,data=annot)
m2 <- lm(rASA ~ GERP_score * pLDDT_class*protein_length,data=annot)
anova(m0,m1)
anova(m1,m2)

# ggplot(annot,aes(rASA,GERP_score))+
#   #geom_point(size=0.01,alpha=0.1)+
#   geom_smooth(size=1, linetype="solid", alpha=0.4,color ="black")+
#   theme_bw(base_size = 12)+
#   facet_grid(~ pLDDT_class)+
#   labs(x = "RSA", y = "GERP")
# ggsave("GERP_RSAGERP_plddtclass.pdf",width=10,height=5,units="in",dpi=1000)

dat33 <- annot %>%
  group_by(pLDDT_class) %>%  
  arrange(rASA) %>%  
  mutate(CDF = ecdf(rASA)(rASA)) %>%  
  ungroup()

ggplot(dat33,aes(x=rASA,y=CDF))+ geom_path()+
  expand_limits(x = 0, y = 0) +
  xlim(0,1)+
  facet_grid(. ~ pLDDT_class)+
  theme_bw(base_size = 12)+
  ylab("Cumulative distribution") + 
  xlab("RSA") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave("class_ExtendFig.2-1.pdf", width = 6, height = 3)



# ggplot(dat,aes(RSA_class,GERP_score, fill = RSA_class))+geom_boxplot(alpha=0.9,width = 0.4) +
#   scale_fill_manual(values = c("#BEB8DC", "#82B0D2")) +
#   theme_bw(base_size = 16)+
#   facet_grid(~ pLDDT_class)+
#   stat_compare_means(method = "t.test") +
#   labs(x = "Residues", y = "GERP")+
#   guides(fill=FALSE)
# ggsave("plDDT_exposed-buried.pdf",width = 8, height = 4)
# 
# ggplot(dat,aes(RSA_class,GERP_score, fill = RSA_class))+geom_boxplot(alpha=0.9,width = 0.4) +
#   scale_fill_manual(values = c("#BEB8DC", "#82B0D2")) +
#   theme_bw(base_size = 20)+
#   facet_grid(length_class ~ pLDDT_class)+
#   stat_compare_means(method = "t.test",label.y = min(dat$GERP_score)+0.5, vjust = 1.5) +
#   labs(x = "Residues", y = "GERP")+
#   guides(fill=FALSE)
# ggsave("Length_plDDT_exposed-buried.pdf",width = 10, height = 7)



#--------------------------------------------------------
# Analysis
#--------------------------------------------------------
# GAM linear regression
fit_rASA <- gam(GERP_score ~ s(rASA, bs="cr") + 
                   as.factor(CHROM),
                 method="REML",
                 data=annot)

fit_pLDDT <- gam(GERP_score ~ s(pLDDT, bs="cr") + 
                   as.factor(CHROM),
                 method="REML",
                 data=annot)

fit_pLDDT.rASA <- gam(GERP_score ~ s(rASA, bs="cr") + 
                        s(pLDDT, bs="cr") + 
                        as.factor(CHROM),
                      method="REML",
                      data=annot)

fit_pLDDTxrASA <- gam(GERP_score ~ s(pLDDT, bs="cr") + 
                        s(rASA, by=pLDDT_class, bs="cr") +
                        as.factor(CHROM),
                      method="REML",
                      data=annot)

# Model comparison
aov_list <- list()
aov_list[["pLDDT"]] <- anova(fit_rASA, fit_pLDDT.rASA, test="F")
aov_list[["rASA"]] <- anova(fit_pLDDT, fit_pLDDT.rASA, test="F")
aov_list[["pLDDTxrASA"]] <- anova(fit_pLDDT.rASA, fit_pLDDTxrASA, test="F")

F_label <- sapply(aov_list, function(aov){
  
  p <- aov$`Pr(>F)`[2]
  
  paste0("F-value = ", round(aov$F[2], 0),
         ifelse(p < 0.001,
                " ***",
                sub("e", "E", paste0(" (P = ", signif(p, 2), ")")))
  )

})

p_label <- sapply(aov_list, function(aov){
  
  res <- aov[2, ]
  
  p <- pf(res$F, res$Df, res$`Resid. Df`, lower.tail=FALSE)
  
  ifelse(p < 2.2e-16, "p-value < 2.2e-16", paste0("p-value = ", signif(p, 2)))
  
})

# Non-linear relationship between pLDDT and rASA
rho_test <- with(annot, cor.test(rASA, pLDDT, method="spearman"))
rho <- signif(rho_test$estimate, 2)
pval <- rho_test$p.value
rho_label <- ifelse(pval < 0.001, " ***", sub("e", "E", paste0("(P = ", signif(pval, 2), ")")))

F1 <- ggplot(annot, aes(x=rASA, y=pLDDT)) +
  geom_point(size=1, alpha=0.01, color="#000066") +
  geom_smooth(method='gam', formula=y ~ s(x, bs = "cr")) +
  annotate("text",
           x=min(annot$rASA),
           y=min(annot$pLDDT),
           hjust="left",
           label=bquote(rho~"="~.(rho)~.(rho_label)),
           size=6) +
  theme_bw(base_size=base_size) +
  labs(x="RSA", y="plDDT")

ggsave(filename = "SupplementaryFig.2-1.pdf",width = 5,height = 5,units = "in",dpi = 1000)

#png(paste0(fig_dir, "F1.png"), res=300, width=1800, height=1800)
#print(F1)
#dev.off()

# Non-linear associations among pLDDT, rASA, and GERP score
pLDDT_DF <- predict_gam(fit_pLDDT.rASA,
                        values=list(rASA=median(annot$rASA),
                                    pLDDT=seq(0, 100, 0.1),
                                    CHROM="1")) %>%
  dplyr::mutate(fit=fit-coef(fit_pLDDT.rASA)["(Intercept)"]) %>%
  dplyr::select(pLDDT, fit, se.fit) %>%
  as.data.frame()

rASA_DF <- predict_gam(fit_pLDDT.rASA,
                       values=list(rASA=seq(0, 1, 0.001),
                                   pLDDT=median(annot$pLDDT),
                                   CHROM="1")) %>%
  dplyr::mutate(fit=fit-coef(fit_pLDDT.rASA)["(Intercept)"]) %>%
  dplyr::select(rASA, fit, se.fit) %>%
  as.data.frame()

F2 <- list(
  ggplot(rASA_DF, aes(x=rASA, y=fit)) +
    geom_smooth_ci(size=1, linetype="solid", ci_z=q97.5, ci_alpha=0.2) +
    annotate("label",
             x=min(rASA_DF$rASA),
             y=min(rASA_DF$fit-rASA_DF$se.fit*q97.5),
             hjust="left",
             label=p_label["rASA"],
             size=6) +
    theme_bw(base_size=base_size) +
    labs(x="RSA", y="Average GERP score deviation", tag="a"),
  ggplot(pLDDT_DF, aes(x=pLDDT, y=fit)) +
    geom_smooth_ci(size=1, linetype="solid", ci_z=q97.5, ci_alpha=0.2) +
    annotate("label",
             x=min(pLDDT_DF$pLDDT),
             y=min(pLDDT_DF$fit-pLDDT_DF$se.fit*q97.5),
             hjust="left",
             label=p_label["pLDDT"],
             size=6) +
    theme_bw(base_size=base_size) +
    labs(x="plDDT", y="Average GERP score deviation", tag="b")

)

ggsave(filename ="SupplementaryFig.2-2A,B.pdf",plot = (F2[[1]] | F2[[2]]),width = 10,height = 5,units = "in",dpi = 1000)
#png(paste0(fig_dir, "F2.png"), res=300, width=3600, height=1800)
#print(
#  (F2[[1]] | F2[[2]])
#)
#dev.off()

# Differential relationship between rASA and GERP score, depending on pLDDT
rASAxpLDDT_DF <- foreach(pLDDT_class=unique(annot$pLDDT_class), .combine=rbind) %do% {
  
  tmp <- annot[annot$pLDDT_class == pLDDT_class, ]
  
  predict_gam(
    fit_pLDDTxrASA,
    values=list(
      rASA=seq(min(tmp$rASA), max(tmp$rASA), 0.001),
      pLDDT=median(tmp$pLDDT),
      pLDDT_class=pLDDT_class,
      CHROM="1"
    )
  )
  
} %>%
  dplyr::mutate(fit=fit-coef(fit_pLDDTxrASA)["(Intercept)"]) %>%
  dplyr::select(rASA, pLDDT_class, fit, se.fit) %>%
  as.data.frame()

F3 <- ggplot(rASAxpLDDT_DF, aes(x=rASA, y=fit)) +
  facet_grid(~ pLDDT_class) +
  #geom_point()+
  geom_smooth_ci(size=1, linetype="solid", ci_z=q97.5, ci_alpha=0.2) +
  theme_bw(base_size=base_size) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        plot.title=element_text(size=base_size)) +
  labs(x=paste0("RSA"),
       y="Average GERP score deviation",
       title="") +
  geom_label(data=data.frame(x=max(rASAxpLDDT_DF$rASA),
                             y=min(rASAxpLDDT_DF$fit-rASAxpLDDT_DF$se.fit*q97.5),
                             pLDDT_class=tail(levels(rASAxpLDDT_DF$pLDDT_class), 1),
                             label=p_label["pLDDTxrASA"]),
             aes(x=x, y=y, label=label),
             hjust="right",
             size=6)
ggsave(filename="F3.pdf",width=10,height=5,units="in",dpi=1000)
#png(paste0(fig_dir, "F3.png"), res=300, width=3600, height=1800)
#print(F3)
#dev.off()
