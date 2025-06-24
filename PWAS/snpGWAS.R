library(tidyverse)
library(data.table)
library(reshape2)
library(foreach)
library(doParallel)
library(rtracklayer)


setwd("/home/song/wangs/snpgwas")

#read
gff_file <- "/home/song/wangs/pwas_pwp/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz"

#pwas <-read.table("./output/GE_GWAS2.assoc.txt",header = T)

cds <- readGFF(gff_file,
               filter=list(type=c("CDS"))) %>% 
  as.data.table()
cds$seqid <- gsub("^chr", "", cds$seqid)
cds$seqid <- as.integer(cds$seqid)    
cds$gene <- unlist(cds$Parent)
cds$gene <- sub("_.*", "", cds$gene)

gff <- readGFF(gff_file,
               filter=list(type=c("mRNA"))) %>% 
  as.data.frame()

gff$gene <- unlist(gff$Parent)

gene <- readGFF(gff_file,
               filter=list(type=c("gene"))) %>% 
  as.data.frame()

gene$gene <- gene$ID


pwas_files <- list.files("./output", pattern = "\\.assoc.txt$", full.names = TRUE)



#run

setkey(cds, seqid, start, end)


pwas_cds <- foreach(pwas_file = pwas_files,.combine = bind_rows, .packages = c("data.table", "dplyr")) %dopar% {
  
  pwas <- fread(pwas_file)  
  snps_dt <- pwas %>%
    mutate(end = ps) %>%
    as.data.table()
  
  setkey(snps_dt, chr, ps, end)
  
  fo <- foverlaps(snps_dt, cds, type = "within")
  
  top_snp_per_gene <- fo %>%
    filter(!is.na(gene)) %>%
    group_by(gene) %>%
    summarise(min_p = min(p_wald, na.rm = TRUE)) %>%
    ungroup()
  merge_data <- merge(top_snp_per_gene,gene,by="gene") 
  
  merge_data$file <- basename(pwas_file)
  merge_data
}


write.table(pwas_cds, file = "all_results.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# read above results.
data <-read.table("/home/song/wangs/snpgwas/all_results.tsv", header = T)


threshold <- (0.05/1421004)

data$file <- gsub("GE_GWAS2.assoc.txt", "DTS_BLUP.Hung_2012_1_nam", data$file)
data$file <- gsub("GE_GWAS3.assoc.txt", "DTA_BLUP.Hung_2012_1_nam", data$file)
data$file <- gsub("GE_GWAS4.assoc.txt", "asi_BLUP.Hung_2012_1_nam", data$file)
data$file <- gsub("GE_GWAS5.assoc.txt", "plantheight_BLUP.Hung_2012_1_nam ", data$file)
data$file <- gsub("GE_GWAS6.assoc.txt", "leaflength_BLUP.Hung_2012_1_nam", data$file)
data$file <- gsub("GE_GWAS7.assoc.txt", "leafwidth_BLUP.Hung_2012_1_nam", data$file)
data$file <- gsub("GE_GWAS8.assoc.txt", "tassprimbranchno_BLUP.Hung_2012_1_nam", data$file)
data$file <- gsub("GE_GWAS9.assoc.txt", "tasslength_BLUP.Hung_2012_1_nam", data$file)
data$file <- gsub("GE_GWAS10.assoc.txt", "upperleafangle_BLUP.Hung_2012_1_nam", data$file)
data$file <- gsub("GE_GWAS11.assoc.txt", "cobdiam_BLUP.Hung_2012_1_nam", data$file)
data$file <- gsub("GE_GWAS12.assoc.txt", "coblength_BLUP.Hung_2012_1_nam", data$file)
data$file <- gsub("GE_GWAS13.assoc.txt", "earrowno_BLUP.Hung_2012_1_nam", data$file)
data$file <- gsub("GE_GWAS14.assoc.txt", "kernelnoperrow_BLUP.Hung_2012_1_nam", data$file)
data$file <- gsub("GE_GWAS15.assoc.txt", "earmass_BLUP.Hung_2012_1_nam", data$file)
data$file <- gsub("GE_GWAS16.assoc.txt", "totalkernelweight_BLUP.Hung_2012_1_nam ", data$file)
data$file <- gsub("GE_GWAS17.assoc.txt", "weight20kernels_BLUP.Hung_2012_1_nam", data$file)
data$file <- gsub("GE_GWAS18.assoc.txt", "totalkernelno_BLUP.Hung_2012_1_nam", data$file)
data$file <- gsub("GE_GWAS19.assoc.txt", "nodenumberbelowear_BLUP.Hung_2012_1_nam", data$file)
data$file <- gsub("GE_GWAS20.assoc.txt", "nodenumberaboveear_BLUP.Hung_2012_1_nam", data$file)
data$file <- gsub("GE_GWAS21.assoc.txt", "numbraceroots_BLUP.Hung_2012_1_nam", data$file)
data$file <- gsub("GE_GWAS22.assoc.txt", "chlorophylla_BLUP.Wallace_2014_nam", data$file)
data$file <- gsub("GE_GWAS23.assoc.txt", "chlorophyllb_BLUP.Wallace_2014_nam", data$file)
data$file <- gsub("GE_GWAS24.assoc.txt", "malate_BLUP.Wallace_2014_nam", data$file)
data$file <- gsub("GE_GWAS25.assoc.txt", "fumarate_BLUP.Wallace_2014_nam", data$file)
data$file <- gsub("GE_GWAS26.assoc.txt", "glutamate_BLUP.Wallace_2014_nam", data$file)
data$file <- gsub("GE_GWAS27.assoc.txt", "aminoacids_BLUP.Wallace_2014_nam", data$file)
data$file <- gsub("GE_GWAS28.assoc.txt", "protein_BLUP.Wallace_2014_nam", data$file)
data$file <- gsub("GE_GWAS29.assoc.txt", "nitrate_BLUP.Wallace_2014_nam", data$file)
data$file <- gsub("GE_GWAS30.assoc.txt", "starch_BLUP.Wallace_2014_nam", data$file)
data$file <- gsub("GE_GWAS31.assoc.txt", "sucrose_BLUP.Wallace_2014_nam", data$file)
data$file <- gsub("GE_GWAS32.assoc.txt", "glucose_BLUP.Wallace_2014_nam", data$file)
data$file <- gsub("GE_GWAS33.assoc.txt", "fructose_BLUP.Wallace_2014_nam", data$file)
data$file <- sub("_B.*", "", data$file)
colnames(data)[14] <- "trait"

data_significance <- data[!is.na(data$min_p) & data$min_p < threshold, ]
data_significance$file <- sub("_B.*", "", data_significance$file)
data_significance2 <- data_significance[,c("file","gene","min_p")]
write.table(data_significance2, file = "data_significance.csv",sep = ',',
             quote = FALSE, row.names = FALSE)

PWAS_file_signficant <- "/home/song/wangs/pwas_pwp/PWAS/QC_hits_max_B73"
pwas<- read.table(PWAS_file_signficant, header = T)
pwas$B73 <- sub("_.*", "", pwas$B73)
pwas$gene <- pwas$B73

pwas_gene <- merge(hits_max,gene_info,by="pangene")

significan_combin_pwas <- merge(pwas,data,by=c("gene","trait"),all.x=TRUE)
significan_combin_pwas <- significan_combin_pwas[,c("trait","gene","Structure","min_p")]
significan_combin_gwas <- merge(data_significance,pwas,by=c("gene","trait"),all.x=TRUE)
significan_combin_gwas <- significan_combin_gwas[,c("trait","gene","Structure","min_p")]

df_merged <- rbind(significan_combin_gwas, significan_combin_pwas) #253
df_unique <- unique(df_merged)

df_unique <- df_unique[complete.cases(df_unique[, 3:4]), ]

df_unique[, 4] <- -log10(df_unique[, 4])  #count:92


p <- ggplot(df_unique,
            aes(x=Structure,
                y=min_p)) +
  geom_point(color="black",alpha=0.6) +
  #geom_abline(intercept=0, slope=1) +
  labs(x=bquote("Structure-based tests :"~-log[10]~"(p)"), 
       y=bquote("SNP-based tests:"~-log[10]~"(p)")) +
  geom_hline(yintercept=-log10(threshold), linetype="dashed") +
  geom_vline(xintercept=-log10(BC_threshold), linetype="dashed") +
  theme_bw(base_size=26)

ggsave(filename = "snpComparison_with_structure-based_tests.pdf",width = 8, height = 7)
