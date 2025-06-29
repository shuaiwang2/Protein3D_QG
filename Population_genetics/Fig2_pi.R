01.split.vcf.sh
cat samples.txt | parallel -j 50 '
    bcftools view -s {} /work/songbaoxing/wangs/alphafold_QG/282/vcf/hmp321_agpv5_chr.biallelic.SNPs_only.vcf.gz -Oz -o ./vcf/{.}.vcf.gz; echo " {} end"
'

#02.vcftofasta_cds.sh
cat samples.txt | /data/wangs/software/bcftools/bin/parallel -j 10 '
python /home/wangs/my_data/software/POPULATION-DIVERSITY/vcf2fa.py ../ref/Zm-B73-REFERENCE-NAM-5.0.fa  ./vcf/{.}.vcf.gz >  ./fasta/{.}.vcf.fa
#~/miniconda3/bin/samtools faidx ./fasta/{.}.vcf.fa
/data/wangs/software/gffread/bin/gffread -x ./cds/{.}.cds -g ./fasta/{.}.vcf.fa ../ref/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3
'

#03.msa.sh
python /home/wangs/my_data/software/POPULATION-DIVERSITY/make_protein_files.py ID.cds.list

```{R}
ls 06concatenatedPi/bin* | /data/wangs/software/bcftools/bin/parallel -j 50 "/usr/bin/time ${bpppopstats} input.sequence.file={} pop.stats=PiN_PiS alphabet='Codon(letter=DNA)' > {}.outfile"
```


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

model <- lm(weighted_pi ~ plddt+plddt:rsa_num, data = df_summary)
summary(model)


ggplot(df_summary,aes(rsa_num,weighted_pi))+
  geom_point(size = 0.8,)+
  geom_line(col="black")+
  #geom_smooth(method = "lm", se = FALSE)  +
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
ggsave("pi.pdf",width=6,height=2.5,units="in",dpi=1000)

Residuals:
      Min        1Q    Median        3Q       Max 
-0.051198 -0.016458 -0.000518  0.015982  0.090405 
 
Coefficients:
                      Estimate Std. Error t value Pr(>|t|)    
(Intercept)            0.50674    0.01745  29.032  < 2e-16 ***
plddt[50,70)          -0.23765    0.02468  -9.627 5.69e-11 ***
plddt[70,90)          -0.37844    0.02468 -15.331 2.69e-16 ***
plddt[90,100]         -0.42401    0.02468 -17.177  < 2e-16 ***
plddt[0,50):rsa_num   -0.04923    0.03027  -1.626    0.114    
plddt[50,70):rsa_num   0.17697    0.03027   5.846 1.69e-06 ***
plddt[70,90):rsa_num   0.17789    0.03027   5.877 1.55e-06 ***
plddt[90,100]:rsa_num  0.22714    0.03027   7.504 1.52e-08 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 
Residual standard error: 0.02749 on 32 degrees of freedom
Multiple R-squared:  0.963,    Adjusted R-squared:  0.955 
F-statistic: 119.1 on 7 and 32 DF,  p-value: < 2.2e-16