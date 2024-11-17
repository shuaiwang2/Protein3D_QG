
# Population genetics using 282 association panel 

## Download, B73 V4 as reference.
```
#282 download https://datacommons.cyverse.org/browse/iplant/home/shared/panzea/hapmap3/hmp321/unimputed/282_libs_2015
```
## Lift over to V5
```
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
python3 /media/shuaiwang/3d/foldingproteins/ProteinSequenceDiversityAnnotation/OverVariantAnnotationsWithPdbsDSSPUsingUpdatedvcf2.py > DSSP.overlap_282.txt    

cp hmp321_agpv5_chr1.biallelic.SNPs_only.vcf hmp321_agpv5.biallelic.SNPs_only.vcf
sed -i s/chr//g hmp321_agpv5.biallelic.SNPs_only.vcf
cat hmp321_agpv5_chr2.biallelic.SNPs_only.vcf | grep -v "^#" >> hmp321_agpv5.biallelic.SNPs_only.vcf
cat hmp321_agpv5_chr3.biallelic.SNPs_only.vcf | grep -v "^#" >> hmp321_agpv5.biallelic.SNPs_only.vcf
cat hmp321_agpv5_chr4.biallelic.SNPs_only.vcf | grep -v "^#" >> hmp321_agpv5.biallelic.SNPs_only.vcf
cat hmp321_agpv5_chr5.biallelic.SNPs_only.vcf | grep -v "^#" >> hmp321_agpv5.biallelic.SNPs_only.vcf
cat hmp321_agpv5_chr6.biallelic.SNPs_only.vcf | grep -v "^#" >> hmp321_agpv5.biallelic.SNPs_only.vcf
cat hmp321_agpv5_chr7.biallelic.SNPs_only.vcf | grep -v "^#" >> hmp321_agpv5.biallelic.SNPs_only.vcf
cat hmp321_agpv5_chr8.biallelic.SNPs_only.vcf | grep -v "^#" >> hmp321_agpv5.biallelic.SNPs_only.vcf
cat hmp321_agpv5_chr9.biallelic.SNPs_only.vcf | grep -v "^#" >> hmp321_agpv5.biallelic.SNPs_only.vcf
cat hmp321_agpv5_chr10.biallelic.SNPs_only.vcf | grep -v "^#" >> hmp321_agpv5.biallelic.SNPs_only.vcf
bcftools sort  hmp321_agpv5.biallelic.SNPs_only.vcf -o hmp321_agpv5.biallelic.SNPs_only.sort.vcf

```
## annot 
```
python3 OverVariantAnnotationsWithPdbsDSSPUsingUpdatedvcf2.py > DSSP.overlap_282.txt                              
```
