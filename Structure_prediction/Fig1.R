library(ggplot2)
library(scales)
library(dplyr)
library(viridis)
library(RColorBrewer)

#fig1A
setwd('C:\\Users\\a\\Desktop\\Protein_QG\\Fig\\Fig1')

sub <- read.csv("E:\\3dres\\fig1\\count_aa_indelandsub\\count_diff_NAM.txt",sep = '\t',header = F)
sub$V2 = factor(sub$V2, levels=c("SNP>1&InDel>0","SNP>1&InDel=0","SNP=1&InDel>0","SNP=0&InDel>0","SNP=1&InDel=0","SNP=0&InDel=0"))
ggplot(sub, aes(x=V1, y=V3, fill=V2, pattern = V2)) + geom_bar(stat="identity")+  
  ylab("Proportion of \n total protein sequences") + 
  xlab("NAM founder lines") +
  labs(title = "",fill = "", color = "") +
  scale_fill_manual(values=c("#D6E1F5","#C8C2E4","#9BBBE1","#BDA0D9","#FF9A9A","#EAEAEA"))+
  scale_color_manual(values=c("#D6E1F5","#C8C2E4","#9BBBE1","#BDA0D9","#FF9A9A","#EAEAEA"))+
  theme(plot.title = element_text(size = 16),
        axis.line = element_line(color ="black"), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linetype = "solid"),
        axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5,colour = "black"))
ggsave(filename = "fig1A.pdf",width = 8,height = 6,units = "in",dpi = 1000)


#1B pie plot
df <- data.frame("category" = c('VeryHigh','Confident','Low','VeryLow'),
                 "amount" = c(9724592,5138394,3217285,6319992))
df$category=factor(df$category,levels=c('VeryHigh','Confident','Low','VeryLow'))
ggplot(df, aes(x="", y=amount, fill=category)) +
  geom_col() +
  coord_polar(theta = "y") +
  theme(axis.text = element_blank(),
        axis.title=element_blank(),
        axis.line = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border =element_blank() )+
  labs(fill = "Scores", color = "") +ggtitle("B73")+
  scale_fill_manual(values=c("#41B6C4","#7FCDBB","#C7E9B4","#EDF8B1"),labels = c("VeryHigh" = "plDDT∈[90,100]", "Confident" = "plDDT∈[70,90)", "Low" = "plDDT∈[50,70)", "VeryLow" = "plDDT∈[0,50)"))+
  geom_text(aes(label = percent(amount / sum(df$amount))), position = position_stack(vjust = 0.5),size =5)
ggsave("pie.b73.aa.pdf", width = 4, height = 4)

df <- data.frame("category" = c('VeryHigh','Confident','Low','VeryLow'),
                 "amount" = c(127472530,61140799,36647854,67556174))
df$category=factor(df$category,levels=c('VeryHigh','Confident','Low','VeryLow'))
ggplot(df, aes(x="", y=amount, fill=category)) +
  geom_col() +
  coord_polar(theta = "y") +
  theme(axis.text = element_blank(),
        axis.title=element_blank(),
        axis.line = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border =element_blank() )+
  labs(fill = "Scores", color = "") +ggtitle("NAM")+
  scale_fill_manual(values=c("#41B6C4","#7FCDBB","#C7E9B4","#EDF8B1"),labels = c("VeryHigh" = "plDDT∈[90,100]", "Confident" = "plDDT∈[70,90)", "Low" = "plDDT∈[50,70)", "VeryLow" = "plDDT∈[0,50)"))+
  geom_text(aes(label = percent(amount / sum(df$amount))), position = position_stack(vjust = 0.5),size =4.5)
ggsave("pie.nam.aa.pdf", width = 4, height = 4)


df <- data.frame("category" = c('VeryHigh','Confident','Low','VeryLow','Unresloved'),
                 "amount" = c(24398,18561,16339,8964,4277))
df$category=factor(df$category,levels=c('VeryHigh','Confident','Low','VeryLow','Unresloved'))
ggplot(df, aes(x="", y=amount, fill=category)) +
  geom_col() +
  coord_polar(theta = "y") +
  theme(axis.text = element_blank(),
        axis.title=element_blank(),
        axis.line = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border =element_blank() )+
  labs(fill = "Scores", color = "") +ggtitle("B73")+
  scale_fill_manual(values=c("#41B6C4","#7FCDBB","#C7E9B4","#EDF8B1","#FFFFD9"),labels = c("VeryHigh" = "plDDT∈[90,100]", "Confident" = "plDDT∈[70,90)", "Low" = "plDDT∈[50,70)", "VeryLow" = "plDDT∈[0,50)"))+
  geom_text(aes(label = paste0(round(amount / sum(amount) * 100, 1), "%")), position = position_stack(vjust = 0.5),size =3)
ggsave("pie.b73.protein.pdf", width = 4, height = 4)


df <- data.frame("category" = c('VeryHigh','Confident','Low','VeryLow','Unresloved'),
                 "amount" = c(322737,212242,184549,76121,293593))
df$category=factor(df$category,levels=c('VeryHigh','Confident','Low','VeryLow','Unresloved'))
ggplot(df, aes(x="", y=amount, fill=category)) +
  geom_col() +
  coord_polar(theta = "y") +
  #scale_fill_brewer()+
  theme(axis.text = element_blank(),
        axis.title=element_blank(),
        axis.line = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border =element_blank() )+
  labs(fill = "Scores", color = "") +ggtitle("NAM")+
  scale_fill_manual(values=c("#41B6C4","#7FCDBB","#C7E9B4","#EDF8B1","#FFFFD9"),labels = c("VeryHigh" = "plDDT∈[90,100]", "Confident" = "plDDT∈[70,90)", "Low" = "plDDT∈[50,70)", "VeryLow" = "plDDT∈[0,50)"))+
  geom_text(aes(label = paste0(round(amount / sum(amount) * 100, 1), "%")), position = position_stack(vjust = 0.55),size =3)
ggsave("pie.nam.protein.pdf", width = 4, height = 4)

#B73 
data4 <- read.table("all.residues.plddt.txt",sep = '\t',header = F)
data4$V5 <- ifelse(data4$V4 >= 90, "VeryHigh",
                   ifelse(data4$V4 >= 70, "Confident",
                          ifelse(data4$V4 >= 50, "Low", "VeryLow")))

data4$V5 = factor(data4$V5, levels=c("VeryHigh","Confident","Low","VeryLow"))
data5<- data4[which(data4$V1 == "B73"),]
max(data5$V4)
min(data5$V4)
#[0,20) is ignored
ggplot(data=data5, aes(x=V4)) +
  geom_histogram(aes(fill = V5),binwidth = 2,boundary = 0,color ="white",closed = "left")+
  scale_fill_manual(values=c("#41B6C4","#7FCDBB","#C7E9B4","#EDF8B1"))+
  ylab("Residues counts") + 
  xlab("plDDT") +
  labs(fill = "") +
  scale_x_continuous(breaks = c(20,50,70,90,100),limits = c(20, 100))+
  scale_y_continuous(labels = label_scientific())+
  theme(axis.line = element_line(color ="black"), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15)  ) 
ggsave(filename = "b73_residues_plddt.pdf",width = 8,height = 5,units = "in",dpi = 1000)

#nam
max(data4$V4)
min(data4$V4)
ggplot(data=data4, aes(x=V4)) +
  geom_histogram(aes(fill = V5),binwidth = 2,boundary = 0,color ="white",closed = "left")+
  scale_fill_manual(values=c("#41B6C4","#7FCDBB","#C7E9B4","#EDF8B1"))+
  ylab("Residues counts") + 
  xlab("plDDT") +
  labs(fill = "") +
  scale_x_continuous(breaks = c(20,50,70,90,100),limits = c(20, 100))+
  scale_y_continuous(labels = label_scientific())+
  theme(axis.line = element_line(color ="black"), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave(filename = "nam_residues_plddt.pdf",width = 8,height = 5,units = "in",dpi = 1000)

#B73
dat = read.table("B73.plddtandrsa.txt")
dat$V5 = factor(dat$V5, levels=c("VeryHigh","Confident","Low","VeryLow"))
ggplot(dat, aes(V4, fill = V5,color=V5)) +
  geom_density(aes(y=after_stat(count)),position = "stack",size= 0.1) +
  ylab("Residues counts") + 
  xlab("RSA") +
  xlim(0,1)+
  labs(fill = "", color = "") +
  scale_fill_manual(values=c("#41B6C4","#7FCDBB","#C7E9B4","#EDF8B1"))+
  scale_color_manual(values=c("#41B6C4","#7FCDBB","#C7E9B4","#EDF8B1"))+
  scale_y_continuous(labels = label_scientific())+
  theme(axis.line = element_line(color ="black"), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "white", linetype = "solid"))
ggsave("B73_RSA.pdf", plot = last_plot(), width = 6, height = 3)

dat = read.table("all.plddtandrsa.txt")
dat$V5 = factor(dat$V5, levels=c("VeryHigh","Confident","Low","VeryLow"))
ggplot(dat, aes(V4, fill = V5,color=V5)) +
  geom_density(aes(y=after_stat(count)),position = "stack",size= 0.1) +
  ylab("Residues counts") + 
  xlab("RSA") +
  xlim(0,1)+
  labs(fill = "", color = "") +
  scale_fill_manual(values=c("#41B6C4","#7FCDBB","#C7E9B4","#EDF8B1"))+
  scale_color_manual(values=c("#41B6C4","#7FCDBB","#C7E9B4","#EDF8B1"))+
  theme(axis.line = element_line(color ="black"), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "white", linetype = "solid"))
ggsave("nam_RSA.pdf", plot = last_plot(), width = 6, height = 3)

#1E heatmap 
#b73
data <- read.csv("all.protein.plddt.id.txt",sep = '\t',header = F)
data$V5 <- sapply(strsplit(data$V3, "_"), function(x) x[1])
data$V5<-as.numeric(data$V5)
data1<- data[which(data$V1 == "plddt/B73"),]
breaks_x <- seq(0, 1200, by = 100)
breaks_y <- seq(0, 100, by = 10)
counts <- table(cut(data1$V5, breaks = breaks_x), cut(data1$V4, breaks = breaks_y))
counts_df <- as.data.frame(as.table(counts))
names(counts_df) <- c("x_bin", "y_bin", "count")

ggplot(counts_df, aes(x = x_bin, y = y_bin, fill = count)) +
  geom_tile() +
  #scale_colour_brewer(palette = "Set1")
  scale_fill_gradientn(colours =c(brewer.pal(7,"Reds")))+
  labs(x = "Protein length", y = "Median plDDT score", fill = "Count",title="B73") +
  scale_y_discrete(breaks = seq(0,100,by = 10))+
  scale_x_discrete(expand = c(0, 0)) +   
  scale_y_discrete(expand = c(0, 0)) +  
  theme(plot.title = element_text(size=15,hjust=0.5),axis.line = element_blank(), panel.grid = element_blank(),panel.background = element_blank(),panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"), axis.text.x = element_text( size=5,colour = "black"), axis.text.y = element_text( size=5,colour = "black"))
ggsave(filename = "protein_plddtb73.pdf",width = 6,height = 3,units = "in",dpi = 1000)

# nam
data <- read.csv("all.protein.plddt.id.txt",sep = '\t',header = F)
data$V5 <- sapply(strsplit(data$V3, "_"), function(x) x[1])
data$V5<-as.numeric(data$V5)
#data1<- data[which(data$V1 == "plddt/B73"),]
breaks_x <- seq(0, 1200, by = 100)
breaks_y <- seq(0, 100, by = 10)
counts <- table(cut(data$V5, breaks = breaks_x), cut(data$V4, breaks = breaks_y))
counts_df <- as.data.frame(as.table(counts))
names(counts_df) <- c("x_bin", "y_bin", "count")

ggplot(counts_df, aes(x = x_bin, y = y_bin, fill = count)) +
  geom_tile() +
  scale_fill_gradientn(colours =c(brewer.pal(7,"Reds")))+
  labs(x = "Protein length", y = "Median plDDT score", fill = "Count",title="NAM") +
  scale_y_discrete(breaks = seq(0,100,by = 10))+
  scale_x_discrete(expand = c(0, 0)) +   
  scale_y_discrete(expand = c(0, 0)) +  
  theme(plot.title = element_text(size=15,hjust=0.5),axis.line = element_blank(), panel.grid = element_blank(),panel.background = element_blank(),panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"), axis.text.x = element_text( size=5,colour = "black"), axis.text.y = element_text( size=5,colour = "black"))
ggsave(filename = "protein_plddtnam.pdf",width = 6,height = 3,units = "in",dpi = 1000)

