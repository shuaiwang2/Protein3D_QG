## Pre-processing  
Projection: "mergedPanGeneTables" is projection result (pangene * individual) using PHG, the procedure will acqurie protein name of offspring according to NAM founder lines.   
GRM: GRM (genetic relationship matrix) is produced by PlINK, ```plink --bfile NAM_v5_PHGv1_all_chroms_vcftools_wanted --recode12 --output-missing-genotype 0 --transpose --out snp```  


Multiple sequence alignment or multiple structure alignment and similarity matrix for each pan-gene (python)  
   MAFFT ```mafft --maxiterate 10000 --localpair ./100/sequences > ./100/sequences.mafft ```  
   Muscle ```muscle -in ./100/sequences > ./100/sequences.muscle ```  
   T-coffee ```t_coffee  ./100/sequences -output fasta_aln -thread 10 -outfile ./100/sequences.t_coffee.fast_aln ```  
   US-align ```USalign -dir ./100/ ./100/list -suffix .pdb -mm 4 -mol prot -o ./100/sugary1 > ./100/alignment ```  

Similarity matrix for multiple sequence alignment is computed by BLOSUM62 score
Similarity matrix for multiple structure alignment is computed by 
 
     
