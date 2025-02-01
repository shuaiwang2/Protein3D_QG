# Protein3D_QG

## 1. Protein structure prediction for 26 diverse maize inbred lines

We predicted protein 3D structure for 26 maize NAM founder lines using AlphaFold v2 with "reduced_dbs" parameter.

RSA was computed using DSSP.


## 2. Population genetics in 282 population

We lift over 282 assocation panel and annotated it using structure feature(pLDDT, RSA) and other population parameters (GERP, DAF).

## 3. PWAS (proteome-wide association) in NAM population

We conducted PWAS using sequence and structure similarity matirx aligned by MAFFT(blosum 62), Muscle(blosum 62), T-Coffee (blosum 62), US-align(TM score) as input to PWAS.


## 4. PWP (proteome-wide prediction) in NAM population 

We conducted PWP usingsequence, haplotype and structure similarity matirx aligned by MAFFT(blosum 62), Muscle(blosum 62), T-Coffee(blosum 62), US-align (TM score).

Haplotype similarity matirx is a categorical variable with 26 levels for all 26 fonuder lines.
