# PWAS: interpretation of results
# Objective : Compare significance from different descriptors of genetic variability
# Created by: Guillaume Ramstein (ramstein@qgg.au.dk)
# Created on: 4/8/2024

#--------------------------------------------------------
# Script parameters
#--------------------------------------------------------
# Working directory
wd <- ifelse(Sys.info()["user"] == "au686885",
             "C:/Users/au686885/OneDrive - Aarhus Universitet/QGG_plants/alphafold_QG/",
             "~/QGG_plants/alphafold_QG/")
dir.create(wd, showWarnings = FALSE)
setwd("/home/song/wangs/pwas_pwp")

# Protein annotation
pangene_file <- "mergedPanGeneTables"

gff_file <- "Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz"

# cyc_file <- "data/corncyc_pathways.20230103"
# 
# uniprot_file <- "data/UniProt.tsv"
# 
# GOslim_file <- "data/goslim_plant.obo"

# PWAS
PWAS_file <- "PWAS/HQK/PWAS.rds"


# pwas1 <- readRDS("PWAS/HQK/PWAS.rds")
# pwas2 <- readRDS("PWAS/HQK2/PWAS.rds")
# pwas <- rbind(pwas1,pwas2)
# saveRDS(pwas, paste0("PWAS/PWAS.rds"))

#sim_metrics <- c("Structure", "Sequence")
sim_metrics <- c("mafftSequence","muscleSequence","tcoffeeSequence","Structure")

#sim_colors <- c("Structure"="#338800", "Sequence"="#6699FF", "Haplotype"="#FF5500")[sim_metrics]

sim_colors <- c("Structure"="#3B549D","muscleSequence"="#D6AEDD","tcoffeeSequence"="#73ABCF","mafftSequence"="#33A2FF")[sim_metrics]

#permutation_file <- "data/PWAS/Z_K/permutation_tests/PWAS_perm=1-20.rds"

# Size of peak window
w <- 5e6

# Output directory
fig_dir <- "figs/PWAS/HQK/"
dir.create(fig_dir, showWarnings=FALSE)
dir.create(paste0(fig_dir, "PWAS_peaks/"), showWarnings=FALSE)
dir.create(paste0(fig_dir, "QQ_plots/"), showWarnings=FALSE)

# out_dir <- "data/PWAS/HQK/output/"
# dir.create(out_dir, showWarnings=FALSE)
# info_file <- paste0(out_dir, "PlantCyc_info.csv")
# GOE_file <- paste0(out_dir, "GO_enrichment.csv")

#--------------------------------------------------------
# Functions
#--------------------------------------------------------
library(tidyverse)
library(data.table)
library(reshape2)
library(foreach)

library(rtracklayer)

#library(ontologyIndex)

#--------------------------------------------------------
# Gene annotation
#--------------------------------------------------------
# B73 gene IDs
gff <- readGFF(gff_file,
               filter=list(type=c("mRNA"))) %>% 
  as.data.frame()

gff$gene <- unlist(gff$Parent)

# B73 gene coordinates
gene_info <- fread(pangene_file, select=1:2) %>%
  dplyr::mutate(pan_gene_id=sub("pan_gene_", "", pan_gene_id)) %>%
  dplyr::rename(pangene=pan_gene_id, ID=B73) %>%
  merge(gff[, c("seqid", "start", "end", "strand", "ID", "gene")], by="ID") %>%
  dplyr::mutate(seqid=as.character(seqid))

#--------------------------------------------------------
# PWAS results
#--------------------------------------------------------
pwas <- readRDS(PWAS_file) %>%
  reshape2::dcast(trait + pangene ~ sim_metric, value.var="pval") %>%
  dplyr::mutate(pangene=as.integer(pangene)) %>%
  dplyr::rename(mafftSequence=mafftSequenceSimilarityMatrix,
                muscleSequence=muscleSequenceSimilarityMatrix,
                tcoffeeSequence=tcoffeeSequenceSimilarityMatrix,
                Structure=usalignMatrix_TMScore) %>%
  merge(gene_info, by="pangene") %>%
  dplyr::arrange(trait, pangene)

# Peak information
BC_threshold <- 0.05/length(unique(pwas$pangene))
min_p.value <- apply(pwas[, sim_metrics], 1, min, na.rm=TRUE)


#significant <- min_p.value < BC_threshold
significant <- ifelse(min_p.value > 0, min_p.value < BC_threshold, FALSE)

idx <- order(min_p.value)
idx <- idx[idx %in% which(significant)]

excluded <- integer(0)

hits <- data.frame()
peaks <- data.frame()

for (i in idx) {
  
  if (! i %in% excluded) {
    
    trait <- pwas$trait[i]
    seqid <- pwas$seqid[i]
    pos <- pwas$start[i]
    
    i_seq <- which(pwas$trait == trait &
                     pwas$seqid == seqid & 
                     pwas$start > pos - w &
                     pwas$start < pos + w) %>%
      sort()
    
    excluded <- c(excluded, i_seq)
    
    # Hits
    hits <- rbind(hits,
                  data.frame(trait=sub("_.+", "", pwas[i, "trait"]),
                             pangene=as.character(pwas[i, "pangene"]),
                             Distance=pwas$start[i_seq]-pos,
                             Structure=-log10(pwas[i_seq, "Structure"]),
                             mafftSequence=-log10(pwas[i_seq, "mafftSequence"]),
                             muscleSequence=-log10(pwas[i_seq, "muscleSequence"]),
                             tcoffeeSequence=-log10(pwas[i_seq, "tcoffeeSequence"])))
    
    # Peaks by metric
    tmp <- pwas[i_seq, ]
    
    out <- foreach(sim_metric=sim_metrics, .combine=rbind) %do% {
      
      j <- which.min(tmp[, sim_metric])
      
      data.frame(trait=sub("_.+", "", pwas[i, "trait"]),
                 hit=as.character(pwas[i, "pangene"]),
                 metric=sim_metric,
                 peak=as.character(tmp[j, "pangene"]),
                 score=-log10(tmp[j, sim_metric]))
      
    }
    
    peaks <- rbind(peaks, out)
    
  }
  
}

hits_molten <- reshape2::melt(hits, measure.vars=sim_metrics, value.name = "Score")

hits_max <- hits_molten %>%
  group_by(trait, pangene, variable) %>%
  summarize(max_score=max(Score, na.rm=TRUE), .groups="drop_last") %>%
  reshape2::dcast(trait + pangene ~ variable, value.var="max_score") %>%
  as.data.frame()

# PWAS peaks
for (trait in unique(hits_molten$trait)) {
  
  print(trait)
  
  p <- ggplot(hits_molten[hits_molten$trait == trait, ],
              aes(x=Distance, group=variable, color=variable, y=Score)) +
    geom_point(alpha=0.6) +
    geom_smooth(method="gam", formula=y ~ s(x, bs = "cs")) +
    facet_wrap(~ pangene, scales="free") +
    labs(color="Protein\nvariable", title=trait) +
    geom_hline(yintercept=-log10(BC_threshold), linetype="dashed") +
    theme_bw(base_size=16) +
    scale_x_continuous(labels = scales::scientific) +
    scale_color_manual(values=sim_colors)
  
  png(paste0(fig_dir, "PWAS_peaks/", trait, ".png"),
      res=300, width=3000, height=2400)
  print(p)
  dev.off()
  
}

# QQ-plots
pwas_molten <- reshape2::melt(pwas, measure.vars=sim_metrics,
                              value.name = "pvalue") %>%
  mutate(trait=sub("_.+", "", trait))

for (trait in unique(pwas_molten$trait)) {
  
  print(trait)
  
  q <- ggplot(pwas_molten[pwas_molten$trait == trait, ] %>%
                group_by(variable) %>%
                mutate(expected=(rank(pvalue, ties.method="first")-.5)/n()),
              aes(color=variable,
                  x=-log10(expected),
                  y=-log10(pvalue))) +
    geom_point(alpha=0.6) +
    geom_abline(intercept=0, slope=1, linetype="longdash") +
    labs(color="Protein\nvariable", title=trait) +
    theme_bw(base_size=10) +
    scale_color_manual(values=sim_colors, 
                       labels = c("Sequence (MAFFT)", "Sequence (MUSCLE)", "Sequence (T-Coffee)","Structure"))+  
    labs(
      x = expression(paste("Expected ", -log[10]~"(p)")),
      y = expression(paste("Observed ", -log[10]~"(p)"))) +
        theme(
          axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12))
  ggsave(filename = paste0(fig_dir, "QQ_plots/", trait, ".png"),width = 5, height = 4)
  
  # png(paste0(fig_dir, "QQ_plots/", trait, ".png"),
  #     res=300, width=3000, height=2400)
  # print(q)
  # dev.off()
  # 
}

pwas_molten$trait[pwas_molten$trait == "DTA"] <- "Days to anthesis"
pwas_molten$trait[pwas_molten$trait == "DTS"] <- "Days to silking"
pwas_molten$trait[pwas_molten$trait == "aminoacids"] <- "Amino acids"
pwas_molten$trait[pwas_molten$trait == "asi"] <- "Anthesis-silking interval"
pwas_molten$trait[pwas_molten$trait == "chlorophylla"] <- "Chlorophyll a"
pwas_molten$trait[pwas_molten$trait == "chlorophyllb"] <- "Chlorophyll b"
pwas_molten$trait[pwas_molten$trait == "cobdiam"] <- "Cob diameter"
pwas_molten$trait[pwas_molten$trait == "coblength"] <- "Cob length"
pwas_molten$trait[pwas_molten$trait == "earmass"] <- "Ear mass"
pwas_molten$trait[pwas_molten$trait == "earrowno"] <- "Ear row number"
pwas_molten$trait[pwas_molten$trait == "fructose"] <- "Fructose"
pwas_molten$trait[pwas_molten$trait == "fumarate"] <- "Fumarate"
pwas_molten$trait[pwas_molten$trait == "glucose"] <- "Glucose"
pwas_molten$trait[pwas_molten$trait == "glutamate"] <- "Glutamate"
pwas_molten$trait[pwas_molten$trait == "kernelnoperrow"] <- "kernel number per row"
pwas_molten$trait[pwas_molten$trait == "leaflength"] <- "Leaf length"
pwas_molten$trait[pwas_molten$trait == "leafwidth"] <- "Leaf width"
pwas_molten$trait[pwas_molten$trait == "malate"] <- "Malate"
pwas_molten$trait[pwas_molten$trait == "nitrate"] <- "Nitrate"
pwas_molten$trait[pwas_molten$trait == "nodenumberaboveear"] <- "Node number above ear"
pwas_molten$trait[pwas_molten$trait == "nodenumberbelowear"] <- "Node number below ear"
pwas_molten$trait[pwas_molten$trait == "numbraceroots"] <- "Number brace roots"
pwas_molten$trait[pwas_molten$trait == "plantheight"] <- "Plant height"
pwas_molten$trait[pwas_molten$trait == "protein"] <- "Protein"
pwas_molten$trait[pwas_molten$trait == "starch"] <- "Starch"
pwas_molten$trait[pwas_molten$trait == "sucrose"] <- "Sucrose"
pwas_molten$trait[pwas_molten$trait == "tasslength"] <- "Tassel length"
pwas_molten$trait[pwas_molten$trait == "tassprimbranchno"] <- "Tassel primary branch number"
pwas_molten$trait[pwas_molten$trait == "totalkernelno"] <- "Total kernel number"
pwas_molten$trait[pwas_molten$trait == "totalkernelweight"] <- "Total kernel weight"
pwas_molten$trait[pwas_molten$trait == "upperleafangle"] <- "Upper leaf angle"
pwas_molten$trait[pwas_molten$trait == "weight20kernels"] <- "Weight 20 kernels"

# hist p-value
ggplot(pwas_molten, aes(x = -log10(pvalue), fill = variable)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  facet_wrap(~ trait, scales = "free") +
  #scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +  
  labs(x = expression(-log[10]~"(p)"), y = "Count") +
  theme_bw(base_size=10) +
  theme(legend.position = "top")+
  xlim(3, NA)

ggsave(filename = paste0(fig_dir, "hist/", "logtrait3", ".pdf"),width = 15, height = 14)

#all qqplot 
# 
# pwas_molten_exp <- pwas_molten %>%
#   group_by(trait, variable) %>%
#   mutate(expected = (rank(pvalue, ties.method = "first") - 0.5) / n())
# 
# 
# qq_plot <- ggplot(pwas_molten_exp,
#                   aes(x = -log10(expected), y = -log10(pvalue), color = variable)) +
#   geom_point(alpha = 0.6, size = 0.8) +
#   geom_abline(intercept = 0, slope = 1, linetype = "longdash") +
#   scale_color_manual(values = sim_colors,
#                      labels = c("Sequence (MAFFT)", "Sequence (MUSCLE)",
#                                 "Sequence (T-Coffee)", "Structure")) +
#   facet_wrap(~ trait, scales = "free") +
#   labs(
#     color = "Protein\nvariable",
#     x = expression(paste("Expected ", -log[10]~"(p)")),
#     y = expression(paste("Observed ", -log[10]~"(p)"))
#   ) +
#   theme_bw(base_size = 10) +
#   theme(
#     legend.position = "top",
#     strip.text = element_text(size = 10),
#     axis.text.x = element_text(size = 8),
#     axis.text.y = element_text(size = 8)
#   )
# 
# ggsave(filename = paste0(fig_dir, "QQ_plots/all_traits_QQ.pdf"),
#        plot = qq_plot, width = 10, height = 8)

#--------------------------------------------------------
# LD analysis
#--------------------------------------------------------
# Significance by distance
p <- ggplot(hits_molten,
            aes(x=Distance*10^(-6), group=variable, color=variable, y=Score)) +
  geom_smooth(method="gam", formula=y ~ s(x, bs = "cs")) +
  labs(color="Protein\nvariable",
       x="Distance from peak (Mb)",
       y=bquote(-log[10]~"(p)"),
       title="Decay of significance over distance") +
  # geom_hline(yintercept=-log10(BC_threshold), linetype="dashed") +
  theme_bw(base_size=24) +
  scale_color_manual(values=sim_colors, 
                     labels = c("Sequence (MAFFT)", "Sequence (MUSCLE)", "Sequence (T-Coffee)","Structure"))

ggsave(filename = paste0(fig_dir,"Average_score_by_distance.pdf"),width = 10, height = 8)

# png(paste0(fig_dir,"Average_score_by_distance.png"),
#     res=300, width=3000, height=2400)
# print(p)
# dev.off()

p <- ggplot(hits_molten,
            aes(x=cut_width(abs(Distance*10^(-6)), w*10^(-6)/10),
                fill=variable,
                y=Score)) +
  geom_boxplot(outlier.size=0.1) +
  labs(color="Protein\nvariable",
       fill="Protein\nvariable",
       x="Distance from peak (Mb)",
       title="Decay of significance over distance") +
  geom_hline(yintercept=-log10(BC_threshold), linetype="dashed") +
  theme_bw(base_size=24) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
  scale_fill_manual(values=sim_colors)

ggsave(filename = paste0(fig_dir,"Score_by_absolute_distance.pdf"),width = 10, height = 8)

# png(paste0(fig_dir,"Score_by_absolute_distance.png"),
#     res=300, width=3000, height=2400)
# print(p)
# dev.off()

# Comparison with other tests
if (all(c("mafftSequence", "Structure") %in% sim_metrics) ) {
  
  hits_max <- hits_max[hits_max$mafftSequence > -log10(BC_threshold) | hits_max$Structure > -log10(BC_threshold), ]
  p <- ggplot(hits_max,
              aes(x=mafftSequence,
                  y=Structure)) +
    geom_point(alpha=0.5) +
    geom_abline(intercept=0, slope=1) +
    labs(x=bquote("Sequence-based (MAFFT) tests :"~-log[10]~"(p)"), 
         y=bquote("Structure-based tests:"~-log[10]~"(p)"),
         title="Significance of associations by protein variable") +
    geom_hline(yintercept=-log10(BC_threshold), linetype="dashed") +
    geom_vline(xintercept=-log10(BC_threshold), linetype="dashed") +
    theme_bw(base_size=24)
  
  ggsave(filename = paste0(fig_dir,"Comparison_with_mafftsequence-based_tests.pdf"),width = 10, height = 8)
  
  # png(paste0(fig_dir,"Comparison_with_mafftsequence-based_tests.png"),
  #     res=300, width=3000, height=2400)
  # print(p)
  # dev.off()
  
}

if (all(c("muscleSequence", "Structure") %in% sim_metrics)) {
  
  hits_max <- hits_max[hits_max$muscleSequence > -log10(BC_threshold) | hits_max$Structure > -log10(BC_threshold), ]
  p <- ggplot(hits_max,
              aes(x=muscleSequence,
                  y=Structure)) +
    geom_point(alpha=0.5) +
    geom_abline(intercept=0, slope=1) +
    labs(x=bquote("Sequence-based (MUSCLE) tests:"~-log[10]~"(p)"), 
         y=bquote("Structure-based tests:"~-log[10]~"(p)"),
         title="Significance of associations by protein variable") +
    geom_hline(yintercept=-log10(BC_threshold), linetype="dashed") +
    geom_vline(xintercept=-log10(BC_threshold), linetype="dashed") +
    theme_bw(base_size=24)
  
  ggsave(filename = paste0(fig_dir,"Comparison_with_musclesequence-based_tests.pdf"),width = 10, height = 8)
  
  # png(paste0(fig_dir,"Comparison_with_musclesequence-based_tests.png"),
  #     res=300, width=3000, height=2400)
  # print(p)
  # dev.off()
  
}

if (all(c("tcoffeeSequence", "Structure") %in% sim_metrics)) {
  
  hits_max <- hits_max[hits_max$tcoffeeSequence > -log10(BC_threshold) | hits_max$Structure > -log10(BC_threshold), ]
  p <- ggplot(hits_max,
              aes(x=tcoffeeSequence,
                  y=Structure)) +
    geom_point(alpha=0.5) +
    geom_abline(intercept=0, slope=1) +
    labs(x=bquote("Sequence-based (T-Coffee) tests:"~-log[10]~"(p)"), 
         y=bquote("Structure-based tests:"~-log[10]~"(p)"),
         title="Significance of associations by protein variable") +
    geom_hline(yintercept=-log10(BC_threshold), linetype="dashed") +
    geom_vline(xintercept=-log10(BC_threshold), linetype="dashed") +
    theme_bw(base_size=24)
  
  ggsave(filename = paste0(fig_dir,"Comparison_with_tcoffeesequence-based_tests.pdf"),width = 10, height = 8)
  
  # png(paste0(fig_dir,"Comparison_with_tcoffeesequence-based_tests.png"),
  #     res=300, width=3000, height=2400)
  # print(p)
  # dev.off()
  
}

if (all(c("mafftSequence", "muscleSequence") %in% sim_metrics)) {
  
  hits_max <- hits_max[hits_max$mafftSequence > -log10(BC_threshold) | hits_max$muscleSequence > -log10(BC_threshold), ]
  p <- ggplot(hits_max,
              aes(x=mafftSequence,
                  y=muscleSequence)) +
    geom_point(alpha=0.5) +
    geom_abline(intercept=0, slope=1) +
    labs(x=bquote("Sequence-based (MAFFT) tests :"~-log[10]~"(p)"), 
         y=bquote("Sequence-based (MUSCLE) tests:"~-log[10]~"(p)"),
         title="Significance of associations by protein variable") +
    geom_hline(yintercept=-log10(BC_threshold), linetype="dashed") +
    geom_vline(xintercept=-log10(BC_threshold), linetype="dashed") +
    theme_bw(base_size=24)
  
  ggsave(filename = paste0(fig_dir,"Comparison_with_mafft_muscle.pdf"),width = 10, height = 8)
  
  # png(paste0(fig_dir,"Comparison_with_mafft_muscle.png"),
  #     res=300, width=3000, height=2400)
  # print(p)
  # dev.off()
  
}

if (all(c("mafftSequence", "tcoffeeSequence") %in% sim_metrics)) {
  
  hits_max <- hits_max[hits_max$mafftSequence > -log10(BC_threshold) | hits_max$tcoffeeSequence > -log10(BC_threshold), ]
  p <- ggplot(hits_max,
              aes(x=mafftSequence,
                  y=tcoffeeSequence)) +
    geom_point(alpha=0.5) +
    geom_abline(intercept=0, slope=1) +
    labs(x=bquote("Sequence-based (MAFFT) tests :"~-log[10]~"(p)"), 
         y=bquote("Sequence-based (T-coffee) tests:"~-log[10]~"(p)"),
         title="Significance of associations by protein variable") +
    geom_hline(yintercept=-log10(BC_threshold), linetype="dashed") +
    geom_vline(xintercept=-log10(BC_threshold), linetype="dashed") +
    theme_bw(base_size=24)
  
  ggsave(filename = paste0(fig_dir,"Comparison_with_mafftTcoffee.pdf"),width = 10, height = 8)
  
  # png(paste0(fig_dir,"Comparison_with_mafftTcoffee.png"),
  #     res=300, width=3000, height=2400)
  # print(p)
  # dev.off()
  
}




if (all(c("Haplotype", "Structure") %in% sim_metrics)) {
  
  hits_max <- hits_max[hits_max$Haplotype > -log10(BC_threshold) | hits_max$Structure > -log10(BC_threshold), ]
  p <- ggplot(hits_max,
              aes(x=Haplotype,
                  y=Structure)) +
    geom_point(alpha=0.5) +
    geom_abline(intercept=0, slope=1) +
    labs(x="Haplotype-based tests: -log10(p)", y="Structure-based tests: -log10(p)") +
    geom_hline(yintercept=-log10(BC_threshold), linetype="dashed") +
    geom_vline(xintercept=-log10(BC_threshold), linetype="dashed") +
    theme_bw(base_size=24)
  
  ggsave(filename = paste0(fig_dir,"Comparison_with_haplotype-based_tests.pdf"),width = 10, height = 8)
  
  # png(paste0(fig_dir,"Comparison_with_haplotype-based_tests.png"),
  #     res=300, width=3000, height=2400)
  # print(p)
  # dev.off()
  
}


pangene_all <- read.table("pan_gene_matrix_v3_cyverse.csv",sep=",",header = T) %>%
  dplyr::mutate(pangene=sub("pan_gene_", "", Pan_gene_ID))

QC_hits_max_B73 <-left_join(hits_max, pangene_all[, c("pangene", "B73")], by="pangene")

QC_hits_max_B73$B73 <- sub("_T", "_P", QC_hits_max_B73$B73)

write.table(#QC_hits_max_B73_mafft,file="PWAS/QC_hits_max_B73_mafft",
  QC_hits_max_B73,file="PWAS/QC_hits_max_B73",sep=" ",quote=FALSE,col.names=T)

hits_max_location <-left_join(QC_hits_max_B73_mafft, gff[, c("seqid","start","end", "transcript_id")], by="transcript_id")

colnames(QC_hits_max_B73_mafft)[7] <- "transcript_id"

transcript_id 

QC_hits_max_B73_mafft <- hits_max[hits_max$mafftSequence  < -log10(BC_threshold) & hits_max$Structure >-log10(BC_threshold), ]
QC_hits_max_B73_mafft1 <- hits_max[hits_max$mafftSequence  > -log10(BC_threshold) & hits_max$Structure >-log10(BC_threshold), ]
QC_hits_max_B73_mafft2 <- hits_max[hits_max$mafftSequence  > -log10(BC_threshold) & hits_max$Structure < -log10(BC_threshold), ]

QC_hits_max_B73_mafft <- QC_hits_max_B73[QC_hits_max_B73$mafftSequence  < -log10(BC_threshold) & QC_hits_max_B73$Structure >-log10(BC_threshold), ]

QC_hits_max_B73_mafft2 <- QC_hits_max_B73[QC_hits_max_B73$mafftSequence  < QC_hits_max_B73$Structure, ]

QC_hits_max_B73_mafft3 <- QC_hits_max_B73[QC_hits_max_B73$Structure >-log10(BC_threshold), ]

QC_hits_max_B73_mafft_candid <- QC_hits_max_B73_mafft2[c(QC_hits_max_B73_mafft2$B73 =="Zm00001eb067740_T001" |
                                                         QC_hits_max_B73_mafft2$B73 =="Zm00001eb117820_T001" | 
                                                         QC_hits_max_B73_mafft2$B73 =="Zm00001eb188370_T001" | 
                                                         QC_hits_max_B73_mafft2$B73 =="Zm00001eb287100_T001"),]

#trait classification 
#1-1 and 1-2
QC_hits_max_B73_mafft <- QC_hits_max_B73[QC_hits_max_B73$mafftSequence < -log10(BC_threshold) & QC_hits_max_B73$Structure >-log10(BC_threshold), ]
#Class1 <- QC_hits_max_B73_mafft[QC_hits_max_B73_mafft$trait %in% c("DTA","DTS","asi","chlorophyllb"),]
Class1 <- QC_hits_max_B73_mafft[QC_hits_max_B73_mafft$trait %in% c("DTA","DTS","asi","cobdiam", "coblength","earmass","earrowno","kernelnoperrow","tasslength","tassprimbranchno","totalkernelweight","weight20kernels"),]
Class2 <- QC_hits_max_B73_mafft[QC_hits_max_B73_mafft$trait %in% c("leaflength","upperleafangle","nodenumberaboveear","nodenumberbelowear","numbraceroots"),]
# repeat Class11 <- QC_hits_max_B73_mafft[QC_hits_max_B73_mafft$trait %in% c("cobdiam", "coblength","earmass","earrowno","kernelnoperrow","tasslength","tassprimbranchno","totalkernelweight","weight20kernels"),]

#2-1 and 2-2
QC_hits_max_B73_mafft2 <- QC_hits_max_B73[QC_hits_max_B73$mafftSequence +1e-10 < QC_hits_max_B73$Structure, ]
Class1 <- QC_hits_max_B73_mafft2[QC_hits_max_B73_mafft2$trait %in% c("DTA","DTS","asi","cobdiam", "coblength","earmass","earrowno","kernelnoperrow","tasslength","tassprimbranchno","totalkernelweight","weight20kernels"),]
Class2 <- QC_hits_max_B73_mafft2[QC_hits_max_B73_mafft2$trait %in% c("leaflength","upperleafangle","nodenumberaboveear","nodenumberbelowear","numbraceroots"),]
#Class11 <- QC_hits_max_B73_mafft2[QC_hits_max_B73_mafft2$trait %in% c("cobdiam", "coblength","earmass","earrowno","kernelnoperrow","tasslength","tassprimbranchno","totalkernelweight","weight20kernels"),]

#3-1 and 3-2
QC_hits_max_B73_mafft3 <- QC_hits_max_B73[QC_hits_max_B73$Structure >-log10(BC_threshold), ]
Class1 <- QC_hits_max_B73_mafft3[QC_hits_max_B73_mafft3$trait %in% c("DTA","DTS","asi","cobdiam", "coblength","earmass","earrowno","kernelnoperrow","tasslength","tassprimbranchno","totalkernelweight","weight20kernels"),]
Class2 <- QC_hits_max_B73_mafft3[QC_hits_max_B73_mafft3$trait %in% c("leaflength","upperleafangle","nodenumberaboveear","nodenumberbelowear","numbraceroots"),]
#Class11 <- QC_hits_max_B73_mafft3[QC_hits_max_B73_mafft3$trait %in% c("cobdiam", "coblength","earmass","earrowno","kernelnoperrow","tasslength","tassprimbranchno","totalkernelweight","weight20kernels"),]






# # Comparison with permutation tests
# null <- readRDS(permutation_file) %>%
#   group_by(trait, pangene) %>%
#   summarize(permutation=min(pval)) %>%
#   merge(pwas[, c("trait", "pangene", "usalignMatrix_TMScore", "haplotype")], by=c("trait", "pangene"))
# 
# plot(-log10(null$permutation), -log10(null$usalignMatrix_TMScore))
# abline(0, 1, col="blue")

#--------------------------------------------------------
# Functional enrichment
#--------------------------------------------------------
# Functional annotations
cyc <- fread(cyc_file)

uniprot <- fread(uniprot_file) %>%
  mutate(`Gene-name`=sub("_.+", "", EnsemblPlants))

# Peak information
peak_info <- left_join(peaks, gene_info[, c("pangene", "gene")], by=c("peak"="pangene")) %>%
  left_join(cyc, by=c("gene"="Gene-name")) %>%
  arrange(trait, hit, metric) %>%
  as.data.frame()

fwrite(peak_info, info_file)

# PWAS information
pwas_info <- readRDS(PWAS_file) %>%
  left_join(gene_info[, c("pangene", "gene")], by="pangene") %>%
  left_join(uniprot, by=c("gene"="Gene-name")) %>%
  arrange(trait, pangene, sim_metric) %>%
  dplyr::select(-c(n, n_missing, permutation)) %>%
  as.data.frame()

# GO information
slim_IDs <- sub("GO:", "", get_ontology(GOslim_file)$id)

GO_terms <- pwas_info$`Gene Ontology (GO)` %>%
  na.omit() %>%
  strsplit("; ") %>% 
  unlist() %>%
  unique()

GO_index <- GO_terms %>%
  sub("]$", "", .) %>%
  strsplit(" [[]GO:") %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::rename(desc=V1, ID=V2) %>%
  cbind(term=.GlobalEnv$GO_terms) %>%
  dplyr::filter(ID %in% .GlobalEnv$slim_IDs)

GO_list <- pwas_info$`Gene Ontology (GO)` %>%
  strsplit("; ")

GO_match <- matrix(0, nrow=nrow(GO_index), ncol=length(GO_list))

obs <- length(GO_list) > 0 & !is.na(GO_list)
GO_match[, obs] <- sapply(GO_list[obs], function(GO) {
  as.numeric(GO_index$term %in% GO)
})

rownames(GO_match) <- GO_index$term

# GO enrichment test
GOE <- foreach(trait=unique(pwas_info$trait), .combine=rbind) %:%
  foreach(sim_metric=unique(pwas_info$sim_metric), .combine=rbind) %do% {
    
    idx <- which(pwas_info$trait == trait & pwas_info$sim_metric == sim_metric)
    fdr <- p.adjust(pwas_info$pval[idx], "fdr")
    
    foreach(GO_term=rownames(GO_match), .combine=rbind) %do% {
      
      tab <- table(GO_match[GO_term, idx], fdr < 0.05)
      
      test <- try(fisher.test(tab), silent=TRUE)
      
      if (any(class(test) == "try-error")) {
        
        out <- data.frame()
        
      } else {
        
        test <- fisher.test(tab)
        
        out <- data.frame(
          trait=trait,
          metric=sim_metric,
          GO_term=GO_term,
          OR=test$estimate,
          pval=test$p.value
        )
        
      }
      
      return(out)
      
    }
    
  }

GOE$selected <- GOE$pval < 0.05/nrow(GO_match)

fwrite(GOE, GOE_file)
