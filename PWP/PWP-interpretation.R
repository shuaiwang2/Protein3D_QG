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
#setwd("/media/shuaiwang/3d/3dres/fig4")

# PWP
PWP_file <- "PWP/PWP_analysis.csv"
#GS_files <- c("LOFO"="/media/shuaiwang/3d/guil/Data/Proteomewide_prediction/LOFO_analysis.csv",
#              "CV"="/media/shuaiwang/3d/guil/Data/Proteomewide_prediction/CV_analysis.csv"
#              )

GS_files <- c("LOFO"="PWP/LOFO_analysis.csv")


VC_names <- c("GRM"="GRM",
              "mafftSequenceSimilarityMatrix"="mafftSequence",
              "muscleSequenceSimilarityMatrix"="muscleSequence",
              "tcoffeeSequenceSimilarityMatrix"="tcoffeeSequence",
              "haplotype"="Haplotype",
              "usalignMatrix_TMScore"="Structure",
              "E"="Error")

# Output directory
fig_dir <- "figs/PWP/"
dir.create(fig_dir, showWarnings=FALSE)

#--------------------------------------------------------
# Functions
#--------------------------------------------------------
library(tidyverse)
library(ggrepel)
library(data.table)
library(reshape2)
library(foreach)

#--------------------------------------------------------
# Variance partition
#--------------------------------------------------------
# PWP results
pwp <- fread(PWP_file) %>%
  mutate(input=gsub("=", "", input),
         trait=sub("_.+", "", trait)) %>%
  as.data.frame()

colnames(pwp)[-(1:3)] <- VC_names[colnames(pwp)[-(1:3)]]

llik <- pwp %>%
  reshape2::dcast(trait ~ input, value.var="llik") %>%
  arrange(trait) %>%
  column_to_rownames(var="trait")

# Describe the likelihood value for different comparison.  First round comments.
write.csv(llik, file = paste0(out_dir, "llik.csv"))

# pval <- data.frame(baseline=pchisq(2*(llik$structure-llik$baseline), df=1, lower.tail=FALSE),
#                    sequence=pchisq(2*(llik$sequence_structure-llik$sequence), df=1, lower.tail=FALSE),
#                    haplotype=pchisq(2*(llik$haplotype_structure-llik$haplotype), df=1, lower.tail=FALSE),
#                    row.names=rownames(llik)
# )

# Variance partition

dat <- pwp %>%
  rowwise() %>%
  mutate(across(4:ncol(pwp), ~ . / sum(c_across(4:ncol(pwp)), na.rm = TRUE) * 100)) %>%
  ungroup()


vc <- reshape2::melt(dat[, -which(colnames(pwp) == "llik")],
                     measure.vars=colnames(pwp)[-(1:3)], value.name="sigma2") %>%
  arrange(trait, input) %>%
  mutate(input=sub("_", "+", input),
         trait=sub("number", "no", trait))
  

vc$trait[vc$trait == "DTA"] <- "Days to anthesis"
vc$trait[vc$trait == "DTS"] <- "Days to silking"
vc$trait[vc$trait == "aminoacids"] <- "Amino acids"
vc$trait[vc$trait == "asi"] <- "Anthesis-silking \n interval"
vc$trait[vc$trait == "chlorophylla"] <- "Chlorophyll a"
vc$trait[vc$trait == "chlorophyllb"] <- "Chlorophyll b"
vc$trait[vc$trait == "cobdiam"] <- "Cob diameter"
vc$trait[vc$trait == "coblength"] <- "Cob length"
vc$trait[vc$trait == "earmass"] <- "Ear mass"
vc$trait[vc$trait == "earrowno"] <- "Ear row number"
vc$trait[vc$trait == "fructose"] <- "Fructose"
vc$trait[vc$trait == "fumarate"] <- "Fumarate"
vc$trait[vc$trait == "glucose"] <- "Glucose"
vc$trait[vc$trait == "glutamate"] <- "Glutamate"
vc$trait[vc$trait == "kernelnoperrow"] <- "kernel number\n per row"
vc$trait[vc$trait == "leaflength"] <- "Leaf length"
vc$trait[vc$trait == "leafwidth"] <- "Leaf width"
vc$trait[vc$trait == "malate"] <- "Malate"
vc$trait[vc$trait == "nitrate"] <- "Nitrate"
vc$trait[vc$trait == "nodenoaboveear"] <- "Node number\n above ear"
vc$trait[vc$trait == "nodenobelowear"] <- "Node number\n below ear"
vc$trait[vc$trait == "numbraceroots"] <- "Number brace \n roots"
vc$trait[vc$trait == "plantheight"] <- "Plant height"
vc$trait[vc$trait == "protein"] <- "Protein"
vc$trait[vc$trait == "starch"] <- "Starch"
vc$trait[vc$trait == "sucrose"] <- "Sucrose"
vc$trait[vc$trait == "tasslength"] <- "Tassel length"
vc$trait[vc$trait == "tassprimbranchno"] <- "Tassel primary \n branch number"
vc$trait[vc$trait == "totalkernelno"] <- "Total kernel \n number"
vc$trait[vc$trait == "totalkernelweight"] <- "Total kernel \n weight"
vc$trait[vc$trait == "upperleafangle"] <- "Upper leaf angle"
vc$trait[vc$trait == "weight20kernels"] <- "Weight 20 kernels"



vc$input[vc$input == "baseline"] <- "Baseline"
vc$input[vc$input == "sequencemafft"] <- "Sequence"
vc$input[vc$input == "structure"] <- "Structure"
vc$input[vc$input == "sequencemafft+structure"] <- "Sequence+Structure"
selected_inputs <- c("Baseline", "Sequence", "Structure", "Sequence+Structure")

vc <- vc[vc$variable != "Error",]
vc$variable = factor(vc$variable, levels=c("Structure","mafftSequence","GRM"))


#"GRM","mafftSequence","Structure","Error"

p <- ggplot(filter(vc, input %in% selected_inputs) %>%
         mutate(input=factor(input, levels=.GlobalEnv$selected_inputs)),
       aes(x=input, y=sigma2, fill=variable)) +
  geom_bar(position="stack", stat="identity") +
  #facet_wrap(~ trait, scales="free_y") +
  facet_wrap(~ trait, nrow=4, labeller = labeller(trait = function(x) gsub("\n", "\n", x))) +
  scale_fill_manual(values=c("#4A73D1","#7AC9E6","#E1E4D6","#C6DBEF","yellow","green"))+
  #scale_fill_manual(values=c("#6BAED6","#C6DBEF","#EDB120","#2171B5"))+
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
  labs(x="Genomic model", y="Variance (%)", fill="Variance\ncomponent")
ggsave(filename = paste0(fig_dir, "variance_partition.pdf"),width = 14, height = 8)


  #png(paste0(fig_dir, "variance_partition.png"), res=300, width=3600, height=3000)
#print(p)
#dev.off()
#--------------------------------------------------------
# Genomic prediction
#--------------------------------------------------------
GS_name <- "LOFO"
for (GS_name in names(GS_files)) {
  
  print(GS_name)
  
  GS_file <- GS_files[GS_name]
  
  out_dir <- paste0(fig_dir, GS_name, "/")
  dir.create(out_dir, showWarnings=FALSE)
  
  # Results
  GS_molten <- fread(GS_file) %>%
    mutate(input=gsub("=", "", input),
           trait=sub("_BLUP.+", "", trait)) %>%
    group_by(trait, input) %>%
    summarize(Corr.mean=mean(Corr), Corr.SE=sd(Corr)/sqrt(n())) %>%
    as.data.frame()
  
  r.se <- GS_molten %>%
    reshape2::dcast(trait ~ input, value.var="Corr.SE")
  
  r <- GS_molten %>%
    reshape2::dcast(trait ~ input, value.var="Corr.mean") %>%
    merge(r.se, by="trait", suffixes=c("", ".SE"))
  
  #ggplot(GS_molten,aes(x=baseline, y=Corr.mean, fill=trait))+geom_boxplot(alpha = 0.35,size=0.1)
  
  r$trait[r$trait == "DTA"] <- "Days to anthesis"
  r$trait[r$trait == "DTS"] <- "Days to silking"
  r$trait[r$trait == "aminoacids"] <- "Amino acids"
  r$trait[r$trait == "asi"] <- "Anthesis-silking interval"
  r$trait[r$trait == "chlorophylla"] <- "Chlorophyll a"
  r$trait[r$trait == "chlorophyllb"] <- "Chlorophyll b"
  r$trait[r$trait == "cobdiam"] <- "Cob diameter"
  r$trait[r$trait == "coblength"] <- "Cob length"
  r$trait[r$trait == "earmass"] <- "Ear mass"
  r$trait[r$trait == "earrowno"] <- "Ear row number"
  r$trait[r$trait == "fructose"] <- "Fructose"
  r$trait[r$trait == "fumarate"] <- "Fumarate"
  r$trait[r$trait == "glucose"] <- "Glucose"
  r$trait[r$trait == "glutamate"] <- "Glutamate"
  r$trait[r$trait == "kernelnoperrow"] <- "kernel number per row"
  r$trait[r$trait == "leaflength"] <- "Leaf length"
  r$trait[r$trait == "leafwidth"] <- "Leaf width"
  r$trait[r$trait == "malate"] <- "Malate"
  r$trait[r$trait == "nitrate"] <- "Nitrate"
  r$trait[r$trait == "nodenumberaboveear"] <- "Node number above ear"
  r$trait[r$trait == "nodenumberbelowear"] <- "Node number below ear"
  r$trait[r$trait == "numbraceroots"] <- "Number brace roots"
  r$trait[r$trait == "plantheight"] <- "Plant height"
  r$trait[r$trait == "protein"] <- "Protein"
  r$trait[r$trait == "starch"] <- "Starch"
  r$trait[r$trait == "sucrose"] <- "Sucrose"
  r$trait[r$trait == "tasslength"] <- "Tassel length"
  r$trait[r$trait == "tassprimbranchno"] <- "Tassel primary branch number"
  r$trait[r$trait == "totalkernelno"] <- "Total kernel number"
  r$trait[r$trait == "totalkernelweight"] <- "Total kernel weight"
  r$trait[r$trait == "upperleafangle"] <- "Upper leaf angle"
  r$trait[r$trait == "weight20kernels"] <- "Weight 20 kernels"
  
  
  t.test(x=r$structure, y=r$sequencemafft, paired = TRUE, alternative = "greater") # p-value = 1.1e-07   
  t.test(x=r$structure, y=r$sequencemuscle, paired = TRUE, alternative = "greater") # p-value = 1.4e-07
  t.test(x=r$structure, y=r$sequencetcoffee, paired = TRUE, alternative = "greater") # p-value = 1.0e-09

  t.test(x=r$sequencemafft_structure, y=r$sequencemafft, paired = TRUE, alternative = "greater") # p-value = 2.7e-07
  t.test(x=r$sequencemuscle_structure, y=r$sequencemuscle, paired = TRUE, alternative = "greater") # p-value = 3.1e-07
  t.test(x=r$sequencemuscle_structure, y=r$sequencetcoffee, paired = TRUE, alternative = "greater") #p-value = 2.4e-07
  
  t.test(x=r$sequencemafft, y=r$sequencemuscle, paired = TRUE, alternative = "greater") #p-value = 0.72
  t.test(x=r$sequencemafft, y=r$sequencetcoffee, paired = TRUE, alternative = "greater") #p-value = 0.077
  
  t.test(x=r$structure, y=r$haplotype, paired = TRUE, alternative = "greater") #p-value = 0.0037
  t.test(x=r$haplotype_structure, y=r$haplotype, paired = TRUE, alternative = "greater") #p-value = 0.00019
  
  # Comparison: structure vs. baseline
  p <- ggplot(r, aes(x=baseline, y=structure, label=trait)) +
    geom_point(size=3, alpha=0.9) +
    geom_errorbar(aes(xmin=baseline-baseline.SE,
                      xmax=baseline+baseline.SE),
                  width=0.005,
                  alpha=0.9) +
    geom_errorbar(aes(ymin=structure-structure.SE,
                      ymax=structure+structure.SE),
                  width=0.005,
                  alpha=0.9) +
    geom_text_repel() +
    geom_abline(slope=1, intercept=0, linetype="longdash") +
    theme_bw(base_size=24) +
    labs(x="Genomics",
         y="Genomics + Protein structures")
         #title="Prediction accuracy (+/- S.E.)"
  ggsave(filename = paste0(out_dir, "structure-vs-baseline.pdf"),width = 10, height = 9)
  #png(paste0(out_dir, "structure-vs-baseline.png"), res=300, width=3000, height=2800)
  #print(p)
  #dev.off()
  
  # Comparison: structure vs. haplotype
  p <- ggplot(r, aes(x=haplotype, y=haplotype_structure, label=trait)) +
    geom_point(size=3, alpha=0.9) +
    geom_errorbar(aes(xmin=haplotype-haplotype.SE,
                      xmax=haplotype+haplotype.SE),
                  width=0.005,
                  alpha=0.9) +
    geom_errorbar(aes(ymin=haplotype_structure-haplotype_structure.SE,
                      ymax=haplotype_structure+haplotype_structure.SE),
                  width=0.005,
                  alpha=0.9) +
    geom_text_repel() +
    geom_abline(slope=1, intercept=0, linetype="longdash") +
    theme_bw(base_size=24) +
    labs(x="Genomics + Protein IBD",
         y="Genomics + Protein IBD + Protein structures")
    #title="Prediction accuracy (+/- S.E.)"
  ggsave(filename = paste0(out_dir, "haplotype_structure-vs-haplotype.pdf"),width = 10, height = 9)
  #png(paste0(out_dir, "haplotype_structure-vs-haplotype.png"), res=300, width=3000, height=2800)
  #print(p)
  #dev.off()
  
  p <- ggplot(r, aes(x=haplotype, y=structure, label=trait)) +
    geom_point(size=3, alpha=0.9) +
    geom_errorbar(aes(xmin=haplotype-haplotype.SE,
                      xmax=haplotype+haplotype.SE),
                  width=0.005,
                  alpha=0.9) +
    geom_errorbar(aes(ymin=structure-structure.SE,
                      ymax=structure+structure.SE),
                  width=0.005,
                  alpha=0.9) +
    geom_text_repel() +
    geom_abline(slope=1, intercept=0, linetype="longdash") +
    theme_bw(base_size=24) +
    labs(x="Genomics + Protein IBD",
         y="Genomics + Protein structures")
  #title="Prediction accuracy (+/- S.E.)"
  ggsave(filename = paste0(out_dir, "structure-vs-haplotype.pdf"),width = 10, height = 9)
  
  #png(paste0(out_dir, "structure-vs-haplotype.png"), res=300, width=3000, height=2800)
  #print(p)
  #dev.off()
  
  # Comparison: structure vs. sequence
  p <- ggplot(r, aes(x=sequencemafft, y=sequencemafft_structure, label=trait)) +
    geom_point(size=3, alpha=0.9) +
    geom_errorbar(aes(xmin=sequencemafft-sequencemafft.SE,
                      xmax=sequencemafft+sequencemafft.SE),
                  width=0.005,
                  alpha=0.9) +
    geom_errorbar(aes(ymin=sequencemafft_structure-sequencemafft_structure.SE,
                      ymax=sequencemafft_structure+sequencemafft_structure.SE),
                  width=0.005,
                  alpha=0.9) +
    geom_text_repel() +
    geom_abline(slope=1, intercept=0, linetype="longdash") +
    theme_bw(base_size=24) +
    labs(x="Genomics + Protein sequences (MAFFT)",
         y="Genomics + Protein sequences + Protein structures")
  #title="Prediction accuracy (+/- S.E.)"
  ggsave(filename = paste0(out_dir, "sequencemafft_structure-vs-sequence.pdf"),width = 10, height = 9)
  
  #png(paste0(out_dir, "sequence_structure-vs-sequence.png"), res=300, width=3000, height=2800)
  #print(p)
  #dev.off()
  
  p <- ggplot(r, aes(x=sequencemuscle, y=sequencemuscle_structure, label=trait)) +
    geom_point(size=3, alpha=0.9) +
    geom_errorbar(aes(xmin=sequencemuscle-sequencemuscle.SE,
                      xmax=sequencemuscle+sequencemuscle.SE),
                  width=0.005,
                  alpha=0.9) +
    geom_errorbar(aes(ymin=sequencemuscle_structure-sequencemuscle_structure.SE,
                      ymax=sequencemuscle_structure+sequencemuscle_structure.SE),
                  width=0.005,
                  alpha=0.9) +
    geom_text_repel() +
    geom_abline(slope=1, intercept=0, linetype="longdash") +
    theme_bw(base_size=24) +
    labs(x="Genomics + Protein sequences (MUSCLE)",
         y="Genomics + Protein sequences + Protein structures")
  #title="Prediction accuracy (+/- S.E.)"
  ggsave(filename = paste0(out_dir, "sequencemuscle_structure-vs-sequence.pdf"),width = 10, height = 9)
  
  p <- ggplot(r, aes(x=sequencetcoffee, y=sequencetcoffee_structure, label=trait)) +
    geom_point(size=3, alpha=0.9) +
    geom_errorbar(aes(xmin=sequencetcoffee-sequencetcoffee.SE,
                      xmax=sequencetcoffee+sequencetcoffee.SE),
                  width=0.005,
                  alpha=0.9) +
    geom_errorbar(aes(ymin=sequencetcoffee_structure-sequencetcoffee_structure.SE,
                      ymax=sequencetcoffee_structure+sequencetcoffee_structure.SE),
                  width=0.005,
                  alpha=0.9) +
    geom_text_repel() +
    geom_abline(slope=1, intercept=0, linetype="longdash") +
    theme_bw(base_size=24) +
    labs(x="Genomics + Protein sequences (T-Coffee)",
         y="Genomics + Protein sequences + Protein structures")
  #title="Prediction accuracy (+/- S.E.)"
  ggsave(filename = paste0(out_dir, "sequencetcoffee_structure-vs-sequence.pdf"),width = 10, height = 9)
  
  p <- ggplot(r, aes(x=sequencemafft, y=structure, label=trait)) +
    geom_point(size=3, alpha=0.9) +
    geom_errorbar(aes(xmin=sequencemafft-sequencemafft.SE,
                      xmax=sequencemafft+sequencemafft.SE),
                  width=0.005,
                  alpha=0.9) +
    geom_errorbar(aes(ymin=structure-structure.SE,
                      ymax=structure+structure.SE),
                  width=0.005,
                  alpha=0.9) +
    geom_text_repel() +
    geom_abline(slope=1, intercept=0, linetype="longdash") +
    theme_bw(base_size=24) +
    labs(x="Genomics + Protein sequences (MAFFT)",
         y="Genomics + Protein structures")
  #title="Prediction accuracy (+/- S.E.)"
  ggsave(filename = paste0(out_dir, "structure-vs-sequencemafft.pdf"),width = 10, height = 9)
  
  p <- ggplot(r, aes(x=sequencemuscle, y=structure, label=trait)) +
    geom_point(size=3, alpha=0.9) +
    geom_errorbar(aes(xmin=sequencemuscle-sequencemuscle.SE,
                      xmax=sequencemuscle+sequencemuscle.SE),
                  width=0.005,
                  alpha=0.9) +
    geom_errorbar(aes(ymin=structure-structure.SE,
                      ymax=structure+structure.SE),
                  width=0.005,
                  alpha=0.9) +
    geom_text_repel() +
    geom_abline(slope=1, intercept=0, linetype="longdash") +
    theme_bw(base_size=24) +
    labs(x="Genomics + Protein sequences (MUSCLE)",
         y="Genomics + Protein structures")
  #title="Prediction accuracy (+/- S.E.)"
  ggsave(filename = paste0(out_dir, "structure-vs-sequencemuscle.pdf"),width = 10, height = 9)
  
  p <- ggplot(r, aes(x=sequencetcoffee, y=structure, label=trait)) +
    geom_point(size=3, alpha=0.9) +
    geom_errorbar(aes(xmin=sequencetcoffee-sequencetcoffee.SE,
                      xmax=sequencetcoffee+sequencetcoffee.SE),
                  width=0.005,
                  alpha=0.9) +
    geom_errorbar(aes(ymin=structure-structure.SE,
                      ymax=structure+structure.SE),
                  width=0.005,
                  alpha=0.9) +
    geom_text_repel() +
    geom_abline(slope=1, intercept=0, linetype="longdash") +
    theme_bw(base_size=24) +
    labs(x="Genomics + Protein sequences (T-Coffee)",
         y="Genomics + Protein structures")
  #title="Prediction accuracy (+/- S.E.)"
  ggsave(filename = paste0(out_dir, "structure-vs-sequencetcoffee.pdf"),width = 10, height = 9)
  
  #png(paste0(out_dir, "structure-vs-sequence.png"), res=300, width=3000, height=2800)
  #print(p)
  #dev.off()
  
  p <- ggplot(r, aes(x=sequencemuscle, y=sequencemafft, label=trait)) +
    geom_point(size=3, alpha=0.9) +
    geom_errorbar(aes(xmin=sequencemuscle-sequencemuscle.SE,
                      xmax=sequencemuscle+sequencemuscle.SE),
                  width=0.005,
                  alpha=0.9) +
    geom_errorbar(aes(ymin=sequencemafft-sequencemafft.SE,
                      ymax=sequencemafft+sequencemafft.SE),
                  width=0.005,
                  alpha=0.9) +
    geom_text_repel() +
    geom_abline(slope=1, intercept=0, linetype="longdash") +
    theme_bw(base_size=24) +
    labs(x="Genomics + Protein sequences (MUSCLE)",
         y="Genomics + Protein sequences (MAFFT)")
  #title="Prediction accuracy (+/- S.E.)"
  ggsave(filename = paste0(out_dir, "sequencemafft-vs-sequencemuscle.pdf"),width = 10, height = 9)
  
  p <- ggplot(r, aes(x=sequencetcoffee, y=sequencemafft, label=trait)) +
    geom_point(size=3, alpha=0.9) +
    geom_errorbar(aes(xmin=sequencetcoffee-sequencetcoffee.SE,
                      xmax=sequencetcoffee+sequencetcoffee.SE),
                  width=0.005,
                  alpha=0.9) +
    geom_errorbar(aes(ymin=sequencemafft-sequencemafft.SE,
                      ymax=sequencemafft+sequencemafft.SE),
                  width=0.005,
                  alpha=0.9) +
    geom_text_repel() +
    geom_abline(slope=1, intercept=0, linetype="longdash") +
    theme_bw(base_size=24) +
    labs(x="Genomics + Protein sequences (T-Coffee)",
         y="Genomics + Protein sequences (MAFFT)")
  #title="Prediction accuracy (+/- S.E.)"
  ggsave(filename = paste0(out_dir, "sequencemafft-vs-sequencetcoffee.pdf"),width = 10, height = 9)
  
  mean((r$sequencemafft_structure / r$sequencemafft - 1))
  mean((r$sequencemuscle_structure / r$sequencemuscle - 1))
  mean((r$sequencetcoffee_structure / r$sequencetcoffee - 1))
  
  mean((r$structure / r$sequencemafft - 1))
  mean((r$structure / r$sequencemuscle - 1))
  mean((r$structure / r$sequencetcoffee - 1))
  
  mean((r$structure / r$haplotype - 1))
  mean((r$haplotype_structure / r$haplotype - 1))
  mean((r$sequencemafft / r$sequencetcoffee - 1))
  mean((r$sequencemafft / r$sequencemuscle - 1))
  
  # p <- ggplot(r, aes(x=100*(structure/sequencemafft - 1))) +
  #   geom_histogram(fill="#94B5D8", color="black", binwidth=1) +
  #   geom_vline(xintercept=0, linetype="longdash") +
  #   theme_bw(base_size=24) +
  #   labs(x="Relative improvement (%)",
  #        y="Count")
  #        #title="Differences in prediction accuracy")
  # ggsave(filename = paste0(out_dir, "relative_improvement_structure-vs-sequencemafft.pdf"),width = 10, height = 9)
  # 
  # #png(paste0(out_dir, "relative_improvement_structure-vs-sequence.png"), res=300, width=3000, height=2800)
  # #print(p)
  # #dev.off()
  # p <- ggplot(r, aes(x=100*(sequencemafft_structure/structure - 1))) +
  #   geom_histogram(fill="#94B5D8", color="black", binwidth=1) +
  #   geom_vline(xintercept=0, linetype="longdash") +
  #   theme_bw(base_size=24) +
  #   labs(x="Relative improvement (%)",
  #        y="Count")
  # #title="Differences in prediction accuracy")
  # ggsave(filename = paste0(out_dir, "relative_improvement_sequencemafft_structure_vs_structure.pdf"),width = 10, height = 9)
  # 
  # p <- ggplot(r, aes(x=100*(haplotype_structure/structure - 1))) +
  #   geom_histogram(fill="#94B5D8", color="black", binwidth=1) +
  #   geom_vline(xintercept=0, linetype="longdash") +
  #   theme_bw(base_size=24) +
  #   labs(x="Relative improvement (%)",
  #        y="Count")
  # #title="Differences in prediction accuracy")
  # ggsave(filename = paste0(out_dir, "relative_improvement_haplotype_structure_vs_structure.pdf"),width = 10, height = 9)
  # 
}




