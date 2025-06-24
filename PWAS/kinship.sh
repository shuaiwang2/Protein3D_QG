Here are the commands we used to construct the kinship matrix:
$ plink2 --vcf NAM_v5_PHGv1_all_chroms_vcftools_wanted.vcf --max-alleles 2 --make-bed --out NAM_v5_PHGv1_all_chroms_vcftools_wanted
 #assign family ID begin
$ cat NAM_v5_PHGv1_all_chroms_vcftools_wanted.fam | awk '{print($2" "$2" "$3" "$4" "$5" "$6)}' > temp.fam
$ mv temp.fam NAM_v5_PHGv1_all_chroms_vcftools_wanted.fam
 #assign family ID end

 #assign SNP ID begin
$ cat NAM_v5_PHGv1_all_chroms_vcftools_wanted.bim | awk '{print($1"\t"$1"_"$4"\t"$3"\t"$4"\t"$5"\t"$6)}'  > temp.header.bim
$ mv temp.header.bim NAM_v5_PHGv1_all_chroms_vcftools_wanted.bim
 #assign SNP ID end
$ plink --bfile NAM_v5_PHGv1_all_chroms_vcftools_wanted --recode12 --output-missing-genotype 0 --transpose --out snp  #reformat for emmax-kin
$ emmax-kin -v -h -d 10 snp  # generate kinship matrix
We counted the number of SNPs and individuals  for constructing the kinship matrix, using the commands:
$ wc -l snp.tfam
 4984 snp.tfam
$ wc -l NAM_v5_PHGv1_all_chroms_vcftools_wanted.bim
42329376 NAM_v5_PHGv1_all_chroms_vcftools_wanted.bim
$ wc -l snp.map
42,329,376 snp.map

The file NAM_v5_PHGv1_all_chroms_vcftools_wanted.vcf is available at: https://doi.org/10.6084/m9.figshare.28349087 as mentioned in the Data access section.
Moreover, to make it reproducible, we also give the commands we used on GitHub (https://github.com/shuaiwang2/Protein3D_QG/tree/main/PWAS).
