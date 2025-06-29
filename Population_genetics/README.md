
# Population genetics using Goodman association panel 


## VCF using B73 V4 as reference and Lift over to V5.
```
#282 download https://datacommons.cyverse.org/browse/iplant/home/shared/panzea/hapmap3/hmp321/unimputed/282_libs_2015

wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/chain_files/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz
gunzip Zm-B73-REFERENCE-NAM-5.0.fa.gz
sed -i 's/>chr/>/g' Zm-B73-REFERENCE-NAM-5.0.fa 

vcftools --gzvcf hmp321_agpv4_chr1.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | vcftools --vcf - --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c > hmp321_agpv4_chr1.biallelic.SNPs_only.vcf.gz &
vcftools --gzvcf hmp321_agpv4_chr2.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | vcftools --vcf - --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c > hmp321_agpv4_chr2.biallelic.SNPs_only.vcf.gz &
vcftools --gzvcf hmp321_agpv4_chr3.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | vcftools --vcf - --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c > hmp321_agpv4_chr3.biallelic.SNPs_only.vcf.gz &
vcftools --gzvcf hmp321_agpv4_chr4.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | vcftools --vcf - --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c > hmp321_agpv4_chr4.biallelic.SNPs_only.vcf.gz &
vcftools --gzvcf hmp321_agpv4_chr5.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | vcftools --vcf - --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c > hmp321_agpv4_chr5.biallelic.SNPs_only.vcf.gz &
vcftools --gzvcf hmp321_agpv4_chr6.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | vcftools --vcf - --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c > hmp321_agpv4_chr6.biallelic.SNPs_only.vcf.gz &
vcftools --gzvcf hmp321_agpv4_chr7.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | vcftools --vcf - --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c > hmp321_agpv4_chr7.biallelic.SNPs_only.vcf.gz &
vcftools --gzvcf hmp321_agpv4_chr8.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | vcftools --vcf - --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c > hmp321_agpv4_chr8.biallelic.SNPs_only.vcf.gz &
vcftools --gzvcf hmp321_agpv4_chr9.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | vcftools --vcf - --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c > hmp321_agpv4_chr9.biallelic.SNPs_only.vcf.gz &
vcftools --gzvcf hmp321_agpv4_chr10.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | vcftools --vcf - --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c > hmp321_agpv4_chr10.biallelic.SNPs_only.vcf.gz &

CrossMap vcf B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain hmp321_agpv4_chr1.biallelic.SNPs_only.vcf.gz Zm-B73-REFERENCE-NAM-5.0.fa hmp321_agpv5_chr1.biallelic.SNPs_only.vcf &
CrossMap vcf B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain hmp321_agpv4_chr2.biallelic.SNPs_only.vcf.gz Zm-B73-REFERENCE-NAM-5.0.fa hmp321_agpv5_chr2.biallelic.SNPs_only.vcf &
CrossMap vcf B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain hmp321_agpv4_chr3.biallelic.SNPs_only.vcf.gz Zm-B73-REFERENCE-NAM-5.0.fa hmp321_agpv5_chr3.biallelic.SNPs_only.vcf &
CrossMap vcf B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain hmp321_agpv4_chr4.biallelic.SNPs_only.vcf.gz Zm-B73-REFERENCE-NAM-5.0.fa hmp321_agpv5_chr4.biallelic.SNPs_only.vcf &
CrossMap vcf B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain hmp321_agpv4_chr5.biallelic.SNPs_only.vcf.gz Zm-B73-REFERENCE-NAM-5.0.fa hmp321_agpv5_chr5.biallelic.SNPs_only.vcf &
CrossMap vcf B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain hmp321_agpv4_chr6.biallelic.SNPs_only.vcf.gz Zm-B73-REFERENCE-NAM-5.0.fa hmp321_agpv5_chr6.biallelic.SNPs_only.vcf &
CrossMap vcf B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain hmp321_agpv4_chr7.biallelic.SNPs_only.vcf.gz Zm-B73-REFERENCE-NAM-5.0.fa hmp321_agpv5_chr7.biallelic.SNPs_only.vcf &
CrossMap vcf B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain hmp321_agpv4_chr8.biallelic.SNPs_only.vcf.gz Zm-B73-REFERENCE-NAM-5.0.fa hmp321_agpv5_chr8.biallelic.SNPs_only.vcf &
CrossMap vcf B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain hmp321_agpv4_chr9.biallelic.SNPs_only.vcf.gz Zm-B73-REFERENCE-NAM-5.0.fa hmp321_agpv5_chr9.biallelic.SNPs_only.vcf &
CrossMap vcf B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain hmp321_agpv4_chr10.biallelic.SNPs_only.vcf.gz Zm-B73-REFERENCE-NAM-5.0.fa hmp321_agpv5_chr10.biallelic.SNPs_only.vcf &

vcftools --vcf hmp321_agpv5_chr1.biallelic.SNPs_only.vcf --freq --out hmp321_agpv5_chr1.biallelic.SNPs_only.af &
vcftools --vcf hmp321_agpv5_chr2.biallelic.SNPs_only.vcf --freq --out hmp321_agpv5_chr2.biallelic.SNPs_only.af &
vcftools --vcf hmp321_agpv5_chr3.biallelic.SNPs_only.vcf --freq --out hmp321_agpv5_chr3.biallelic.SNPs_only.af &
vcftools --vcf hmp321_agpv5_chr4.biallelic.SNPs_only.vcf --freq --out hmp321_agpv5_chr4.biallelic.SNPs_only.af &
vcftools --vcf hmp321_agpv5_chr5.biallelic.SNPs_only.vcf --freq --out hmp321_agpv5_chr5.biallelic.SNPs_only.af &
vcftools --vcf hmp321_agpv5_chr6.biallelic.SNPs_only.vcf --freq --out hmp321_agpv5_chr6.biallelic.SNPs_only.af &
vcftools --vcf hmp321_agpv5_chr7.biallelic.SNPs_only.vcf --freq --out hmp321_agpv5_chr7.biallelic.SNPs_only.af &
vcftools --vcf hmp321_agpv5_chr8.biallelic.SNPs_only.vcf --freq --out hmp321_agpv5_chr8.biallelic.SNPs_only.af &
vcftools --vcf hmp321_agpv5_chr9.biallelic.SNPs_only.vcf --freq --out hmp321_agpv5_chr9.biallelic.SNPs_only.af &
vcftools --vcf hmp321_agpv5_chr10.biallelic.SNPs_only.vcf --freq --out hmp321_agpv5_chr10.biallelic.SNPs_only.af &

python3 /media/shuaiwang/3d/foldingproteins/ProteinSequenceDiversityAnnotation/annotateVariants2.py /media/shuaiwang/3d/282/hmp321_agpv5_chr1.biallelic.SNPs_only.af.frq > chr1.frq & 
python3 /media/shuaiwang/3d/foldingproteins/ProteinSequenceDiversityAnnotation/annotateVariants2.py /media/shuaiwang/3d/282/hmp321_agpv5_chr2.biallelic.SNPs_only.af.frq > chr2.frq & 
python3 /media/shuaiwang/3d/foldingproteins/ProteinSequenceDiversityAnnotation/annotateVariants2.py /media/shuaiwang/3d/282/hmp321_agpv5_chr3.biallelic.SNPs_only.af.frq > chr3.frq & 
python3 /media/shuaiwang/3d/foldingproteins/ProteinSequenceDiversityAnnotation/annotateVariants2.py /media/shuaiwang/3d/282/hmp321_agpv5_chr4.biallelic.SNPs_only.af.frq > chr4.frq & 
python3 /media/shuaiwang/3d/foldingproteins/ProteinSequenceDiversityAnnotation/annotateVariants2.py /media/shuaiwang/3d/282/hmp321_agpv5_chr5.biallelic.SNPs_only.af.frq > chr5.frq & 
python3 /media/shuaiwang/3d/foldingproteins/ProteinSequenceDiversityAnnotation/annotateVariants2.py /media/shuaiwang/3d/282/hmp321_agpv5_chr6.biallelic.SNPs_only.af.frq > chr6.frq & 
python3 /media/shuaiwang/3d/foldingproteins/ProteinSequenceDiversityAnnotation/annotateVariants2.py /media/shuaiwang/3d/282/hmp321_agpv5_chr7.biallelic.SNPs_only.af.frq > chr7.frq & 
python3 /media/shuaiwang/3d/foldingproteins/ProteinSequenceDiversityAnnotation/annotateVariants2.py /media/shuaiwang/3d/282/hmp321_agpv5_chr8.biallelic.SNPs_only.af.frq > chr8.frq & 
python3 /media/shuaiwang/3d/foldingproteins/ProteinSequenceDiversityAnnotation/annotateVariants2.py /media/shuaiwang/3d/282/hmp321_agpv5_chr9.biallelic.SNPs_only.af.frq > chr9.frq & 
python3 /media/shuaiwang/3d/foldingproteins/ProteinSequenceDiversityAnnotation/annotateVariants2.py /media/shuaiwang/3d/282/hmp321_agpv5_chr10.biallelic.SNPs_only.af.frq > chr10.frq &

cat chr*.frq >all.freq 
python3 OverVariantAnnotationsWithPdbsDSSPUsingUpdatedvcf2.py > DSSP.overlap_282.txt    
```
## Download GERP score 
*out_bed.bed* 
 https://datadryad.org/stash/dataset/doi:10.5061/dryad.70t85k2     # download Zea_mays.allCh.rates.txt
```
#V4 liftover to V5
Zea_mays.allCh.rates.bed to  Zea_mays.allCh.rates.txt # space to Tab
awk '{print $1,$2,$2+200,$4}' Zea_mays.allCh.rates.txt > Zea_mays.allCh.rates.bed
# $2 is mistake coordinate to bed format
CrossMap bed /media/shuaiwang/3d/282/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain /media/shuaiwang/3d/outgroup/doi_10_5061_dryad_70t85k2__v20181212/GERP/GERP/Zea_mays.allCh.rates.bed out_bed.bed
```
## derived allele frequency 
*TIL11.gvcf TIL18.gvcf luxurians.gvcf diploperennis.gvcf*
```
ref=Zm-B73-REFERENCE-NAM-5.0.fa
gff=Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3
anchorwave gff2seq -r ${ref} -i ${gff} -o cds.fa
minimap2 -x splice -t 100 -k 12 -a -p 0.4 -N 20 Zv-TIL11-REFERENCE-PanAnd-1.0.fa cds.fa > TIL11.sam
minimap2 -x splice -t 100 -k 12 -a -p 0.4 -N 20 ${ref} cds.fa > ref.sam
anchorwave proali  -i ${gff} -as cds.fa -r ${ref} -a TIL11.sam -ar ref.sam -s Zv-TIL11-REFERENCE-PanAnd-1.0.fa -n TIL11.anchors -o TIL11.maf -f TIL11.f.maf -t 10 -R 1 -Q 1
perl /home/wangs/my_data/soft/tasseladmin-tassel-5-standalone-846381e171c8/run_pipeline.pl -Xmx250G -debug -MAFToGVCFPlugin -referenceFasta Zm-B73-REFERENCE-NAM-5.0.fa -mafFile TIL11.maf -sampleName TIL11 -gvcfOutput TIL11.gvcf -fillGaps false > TIL11.txt
```
## piN/piS
Fig2_pi.R

## Figure Fig. 2 Fig2.R
