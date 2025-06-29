# Structure_prediction

Fig1.R: all code used  in plotting Fig. 1.
count_aa_indelandsub.py: count pairwisely difference of MSA result (MAFFT). 

1. The predicted protein structures were performed using AlphaFold v2.0.0 with parameter preset as reduced_dbs as described in the AlphaFold database publication on the Frontera supercomputer.  
2. Sequences with residue codes ‘*’ (premature stop codon) or ‘X’ (unknown residue) are also excluded.
3. We optimized the AlphaFold2 pipeline to increase its computational efficiency on this large-scale cluster. This optimization does not generate any different folding.
 Deatiled on Methods
