library(ggplot2)
library(dplyr)
library(ggpubr)
library(agricolae)
library(car)

# Input
setwd("C:\\Users\\a\\Desktop\\Protein_QG\\Fig\\Fig2")
dat = read.table("C:\\Users\\a\\Desktop\\Protein_QG\\Fig\\Fig2\\DSSP.overlap_282.txt")
colnames(dat) = c("SNP_type", "chr", "position", "allele_count", "transcript", "protein_length", "protein_seq_id", "refAle", "altAle", "outRefAle", "outAle1", "outAle2","outAle3", "outAle4", "codon_index", "refFrq", "altFrq", "refCodon", "refAA", "altCodon", "altAA", "blosum62_score", "cds_index", "GERP_score", "plDDT", "ASA", "maxASA")

# data processing
#dat = dat[-which(dat$SNP_type == "nonsynonymous" & (dat$blosum62_score <= 0)),]
#dat = dat[which(dat$allele_count>50),]
#remove the start codon and the last codon, since they might have biase for rASA calculation
#dat = dat[which(dat$codon_index>1),]
#dat = dat[which(dat$codon_index< dat$protein_length),]

dat$GERP_score = as.numeric(dat$GERP_score)
dat$rasa <- pmin(dat$ASA / dat$maxASA, 1)
dat <- dat %>%mutate(plddt_class=cut(plDDT, c(0, 50, 70, 90, 100),right = FALSE))

#Fisher.test to synonymous and nonsynonymous SNP.
nrow(dat[dat$SNP_type == "nonsynonymous" & dat$plDDT >=70, ])
nrow(dat[dat$SNP_type == "nonsynonymous" & dat$plDDT <70, ])
nrow(dat[dat$SNP_type == "synonymous" & dat$plDDT >=70, ])
nrow(dat[dat$SNP_type == "synonymous" & dat$plDDT <70, ])
matrix_plddt <- matrix(
  c(220977,274380,311283, 179532),
  nrow = 2,
  dimnames = list(
    plDDT = c("higherplDDT", "smallerplDDT"),
    variation = c("nonsynonymous", "synonymous")
  )
)
fisher.test(matrix_plddt)
#  p-value < 2.2e-16

nrow(dat[dat$SNP_type == "synonymous" & dat$rasa >=0.5 & dat$rasa <=1, ])
nrow(dat[dat$SNP_type == "synonymous" & dat$rasa < 0.5, ])
nrow(dat[dat$SNP_type == "nonsynonymous" & dat$rasa >=0.5 & dat$rasa <=1 , ])
nrow(dat[dat$SNP_type == "nonsynonymous" & dat$rasa < 0.5, ])
matrix_rsa <- matrix(
  c(203863, 286022,293625, 199481),
  nrow = 2,
  dimnames = list(
    RSA = c("higherrsa", "smallerrsa"),
    variation = c("synonymous", "nonsynonymous")
  )
)
fisher.test(matrix_rsa)
#  p-value < 2.2e-16

nrow(dat[dat$SNP_type == "nonsynonymous", ])
nrow(dat[dat$SNP_type == "synonymous", ])

ggplot(data=dat, aes(x=rasa, group=SNP_type, fill=SNP_type))+
  geom_histogram(breaks=seq(0, 1, 0.025), alpha = 0.35, position = "identity", aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]), ..count..[..group..==2]/sum(..count..[..group..==2]))) ) + 
  scale_fill_manual(values=c("#66A61E", "#7570B3"))+
  labs(x="RSA", y="Proportion of SNPs", title="", fill="") +
  theme(panel.border = element_rect(color = "black", fill = NA),
        #axis.line = element_line(color ="black"), 
        panel.background = element_blank(), 
        legend.position = c(0.85, 0.85) )

ggsave(filename = "rasa_distri.pdf",width = 7,height = 5,units = "in",dpi = 1000)

ggplot(data=dat, aes(x=plDDT, group=SNP_type, fill=SNP_type))+
  geom_histogram(breaks=seq(20, 100, 2), alpha = 0.35, position = "identity", aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]), ..count..[..group..==2]/sum(..count..[..group..==2]))) ) + 
  scale_x_continuous(breaks = c(20,50,70,90,100),limits = c(20, 100))+
  scale_fill_manual(values=c("#66A61E", "#7570B3")) +
  labs(x="pLDDT", y="Proportion of SNPs", title="", fill="") + 
  theme(panel.border = element_rect(color = "black", fill = NA),
        #axis.line = element_line(color ="black"), 
        panel.background = element_blank(), 
        legend.position = c(0.20, 0.85) )
ggsave(filename = "plddt_distri.pdf",width = 7,height = 5,units = "in",dpi = 1000)


# Cumulative distribution
dat3 <- dat %>%
  arrange(rasa) %>%
  mutate(CDF = ecdf(rasa)(rasa))
ggplot(dat3,aes(x=rasa,y=CDF))+ geom_path()+
  expand_limits(x = 0, y = 0) +
  xlim(0,1)+
  theme_bw(12)+
  ylab("Cumulative distribution") + 
  xlab("RSA") 
ggsave("ExtendFig.2-1.pdf", width = 6, height = 3)


## Distribution of DAF.
dat$maf = dat$refFrq
dat[which(dat$maf > dat$altFrq),]$maf = dat[which(dat$maf > dat$altFrq),]$altFrq

dat$daf =0
dat[which((dat$outAle1 == dat$refAle | dat$outAle1 == "R" ) | ( dat$outAle1 == "N" & (dat$outAle2 == dat$refAle | dat$outAle2 == "R") )),]$daf = dat[which((dat$outAle1 == dat$refAle | dat$outAle1 == "R" ) | ( dat$outAle1 == "N" & (dat$outAle2 == dat$refAle | dat$outAle2 == "R") )),]$altFrq
dat[which((dat$outAle1 == dat$altAle) | (dat$outAle1 == "N"  & dat$outAle2 == dat$altAle )),]$daf = dat[which((dat$outAle1 == dat$altAle) | (dat$outAle1 == "N"  & dat$outAle2 == dat$altAle )),]$refFrq
dat=dat[-which(dat$daf == 0),]

dat$daf_bin=cut(dat$daf,c(seq(from = 0, to = 1, by = 0.025)), include.lowest = FALSE)
# Fitting linear model to compare the difference of DAF
m0 <- lm(rasa ~ SNP_type,data = dat)
m1 <- lm(rasa ~ SNP_type + daf_bin,data = dat)
m2 <- lm(rasa ~ SNP_type * daf_bin,data = dat)
anova(m0, m1) # significance of DAF bins
anova(m1, m2) # significance of DAF bins depending on NS


dat2 <- dat[dat$SNP_type == "nonsynonymous",]

model <- aov(rasa ~ daf_bin, data = dat2)
summary(model)
lsd_result <- LSD.test(model, "daf_bin", alpha=0.05,console = TRUE)
print(lsd_result$groups)

model <- aov(rasa ~ daf_bin*plddt_class, data = dat2)
summary(model)
lsd_result <- LSD.test(model, "daf_bin", alpha=0.05,console = TRUE)
print(lsd_result$groups)

model <- lm(rasa ~ daf_bin * plddt_class, contrasts=list(daf_bin=contr.sum, plddt_class =contr.sum), data = dat2)
Anova(model, type=3)


dat %>%
  # Add a new column called 'bin': cut the initial 'carat' in bins
  mutate( bin=cut(daf,c(seq(from = 0, to = 1, by = 0.025)), include.lowest = FALSE) ) %>%  # plot
  ggplot(aes(x=bin, y=rasa, fill=SNP_type) ) +
  geom_boxplot(alpha = 0.35,size=0.1)+
  #facet_grid(RSA_class ~ plddt_class, scales = "free_y") +
  scale_fill_manual(values=c("#66A61E", "#7570B3"))+
  xlab("Derived allele frequency") + ylab("RSA")+theme(axis.line = element_blank(), panel.background = element_blank(),
                                                       panel.border = element_rect(fill=NA,color="black", size=0.1, linetype="solid"), 
                                                       axis.text.x = element_text(angle=270, hjust=0, vjust=1, colour = "black"))+
  annotate("text", x = length(levels(cut(dat$daf, seq(from = 0, to = 1, by = 0.025)))), y = 0, label = "p-value < 2.2e-16", hjust = 1, vjust = 0)

ggsave(filename = "daffreq_0025.pdf",width = 15,height = 5,units = "in",dpi = 1000)


## Different DAF class.
dat$daf_bin=cut(dat$daf,c(0, 0.01, 0.025, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.975,0.99,1))

# Fitting linear model to compare the difference of DAF
m0 <- lm(rasa ~ SNP_type,data = dat)
m1 <- lm(rasa ~ SNP_type + daf_bin,data = dat)
m2 <- lm(rasa ~ SNP_type * daf_bin,data = dat)
anova(m0, m1) # significance of DAF bins
anova(m1, m2) # significance of DAF bins depending on NS

for (class in unique(dat$plddt_class)) {
  
  subset_df <- dat[dat$plddt_class == class & dat$SNP_type =="nonsynonymous" , ]
  model <- aov(rasa ~ daf_bin, data = subset_df)
  print(summary(model))
  lsd_result <- LSD.test(model, "daf_bin", console = TRUE)
}

dat2 <- dat[dat$SNP_type == "nonsynonymous",]
#& dat$RSA_class =="exposed"
model <- aov(rasa ~ daf_bin, data = dat2)
summary(model)
lsd_result <- LSD.test(model, "daf_bin", alpha=0.05,console = TRUE)
print(lsd_result$groups)

model <- aov(rasa ~ daf_bin*plddt_class, data = dat2)
summary(model)
lsd_result <- LSD.test(model, "daf_bin", alpha=0.05,console = TRUE)
print(lsd_result$groups)



dat %>%
  # Add a new column called 'bin': cut the initial 'carat' in bins
  mutate( bin= cut(
    daf,
    breaks = c(0, 0.01, 0.025, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.975,0.99,1)) ) %>%  # plot
  ggplot(aes(x=bin, y=rasa, fill=SNP_type) ) +
  geom_boxplot(alpha = 0.35,size=0.1)+
  facet_grid(.~plddt_class)+
  scale_fill_manual(values=c("#66A61E", "#7570B3"))+
  xlab("Derived allele frequency") + ylab("RSA")+theme(axis.line = element_blank(), panel.background = element_blank(),
                                                       panel.border = element_rect(fill=NA,color="black", size=0.1, linetype="solid"), 
                                                       axis.text.x = element_text(angle=270, hjust=0, vjust=1, colour = "black"))
ggsave(filename = "daffreq_class.pdf",width = 15,height = 5,units = "in",dpi = 1000)


