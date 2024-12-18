# PWAS: association testing
# Objective : Get significance metric for PWAS
# Created by: Guillaume Ramstein (ramstein@qgg.au.dk)
# Created on: 12/12/2023

#--------------------------------------------------------
# Script parameters
#--------------------------------------------------------
# Working directory

setwd("/home/song/wangs/pwas_pwp")

# Number of permutations
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 1) {
  k_list <- as.integer(args[1])
  cat("Permutation:", k_list, "\n")
} else {
  k_list <- as.integer(args[1]):as.integer(args[2])
  cat("Permutations:", k_list[1], "to", k_list[length(k_list)], "\n")
}

include_family_structure <- TRUE

# Input
GRM_file <- "snp.hBN.kinf"  #E:\NAMpopulation_alphafold\NAM_natural_variation
ID_file <- "snp.tfam"

proteome_dir <- "/home/song/wangs/pwas_pwp/seq_structure/"
proteome_file <- "mergedPanGeneTables"
N <- 26

PPC_file <- "QC/PPC_list.rds"

PRM_file <- "PWP/PRM.rds"
sim_metric <- "haplotype"

pheno_file <-"all_NAM_phenos.csv"
traits <- c("DTS_BLUP;Hung_2012_1_nam",
            "DTA_BLUP;Hung_2012_1_nam",
            "asi_BLUP;Hung_2012_1_nam",
            "plantheight_BLUP;Hung_2012_1_nam",
            "leaflength_BLUP;Hung_2012_1_nam",
            "leafwidth_BLUP;Hung_2012_1_nam",
            "tassprimbranchno_BLUP;Hung_2012_1_nam",
            "tasslength_BLUP;Hung_2012_1_nam",
            "upperleafangle_BLUP;Hung_2012_1_nam",
            "cobdiam_BLUP;Hung_2012_1_nam",
            "coblength_BLUP;Hung_2012_1_nam",
            "earrowno_BLUP;Hung_2012_1_nam",
            "kernelnoperrow_BLUP;Hung_2012_1_nam",
            "earmass_BLUP;Hung_2012_1_nam",
            "totalkernelweight_BLUP;Hung_2012_1_nam",
            "weight20kernels_BLUP;Hung_2012_1_nam",
            "totalkernelno_BLUP;Hung_2012_1_nam",
            "nodenumberbelowear_BLUP;Hung_2012_1_nam",
            "nodenumberaboveear_BLUP;Hung_2012_1_nam",
            "numbraceroots_BLUP;Hung_2012_1_nam",
            "chlorophylla_BLUP;Wallace_2014_nam",
            "chlorophyllb_BLUP;Wallace_2014_nam",
            "malate_BLUP;Wallace_2014_nam",
            "fumarate_BLUP;Wallace_2014_nam",
            "glutamate_BLUP;Wallace_2014_nam",
            "aminoacids_BLUP;Wallace_2014_nam",
            "protein_BLUP;Wallace_2014_nam",
            "nitrate_BLUP;Wallace_2014_nam",
            "starch_BLUP;Wallace_2014_nam",
            "sucrose_BLUP;Wallace_2014_nam",
            "glucose_BLUP;Wallace_2014_nam",
            "fructose_BLUP;Wallace_2014_nam"
)

# traits <- c(
#   "DTS_BLUP;Hung_2012_1_nam",
#   "plantheight_BLUP;Hung_2012_1_nam",
#   "weight20kernels_BLUP;Hung_2012_1_nam",
#   "totalkernelno_BLUP;Hung_2012_1_nam"
# )

# Output
output_dir <- "PWAS/HQK/"
dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)

#--------------------------------------------------------
# Functions
#--------------------------------------------------------
library(tidyverse)
library(foreach)
library(data.table)
library(RcppEigen)
library(Matrix)

library(qgg)

permute <- function(mat, k) {
  
  if (k > 0) {
    
    set.seed(k)
    
    rownames(mat) <- sample(rownames(mat))
    
  }
  
  return(mat)
  
}

to_PD <- function(x) {
  
  as.matrix(nearPD(x)$mat)
  
}

import_sim <- function(fn) {
  
  sim <- fread(fn) %>%
    column_to_rownames(var="V1") %>%
    as.matrix()
  
  colnames(sim) <- rownames(sim)
  
  return(sim)
  
}

match_to_sim <- function(id_vec, sim, sep=";") {
  
  sim_entries <- rownames(sim)
  
  id_list <- lapply(strsplit(id_vec, sep), function(id) {
    
    id[id %in% sim_entries]
    
  })
  
  id_DF <- do.call(rbind, id_list) %>%
    as.data.frame()
  
  out <- foreach(id=id_DF, .combine=`+`) %do% {
    
    sim[match(id, rownames(sim)), match(id, colnames(sim))]
    
  } / ncol(id_DF)
  
  rownames(out) <- colnames(out) <- names(id_vec)
  
  return(out)
  
}

match_to_matrix <- function(id_vec, mat, sep=";") {
  
  entries <- rownames(mat)
  
  id_list <- lapply(strsplit(id_vec, sep), function(id) {
    
    id[id %in% entries]
    
  })
  
  id_list[lengths(id_list) == 0] <- NA
  
  id_DF <- do.call(rbind, id_list) %>%
    as.data.frame()
  
  out <- foreach(id=id_DF, .combine=`+`) %do% {
    
    mat[match(id, entries), , drop=FALSE]
    
  } / ncol(id_DF)
  
  rownames(out) <- names(id_vec)
  
  return(out)
  
}

fastWald.test <- function(X, y, L) {
  
  # OLS
  fit <- fastLmPure(X, y)
  
  # Adjusting for rank deficiency
  idx <- which(is.finite(fit$coefficients))
  L <- L[, idx, drop=FALSE]
  L <- L[rowSums(L != 0) > 0, , drop=FALSE]
  
  X <- X[, idx, drop=FALSE]
  b <- fit$coefficients[idx]
  
  # Variance of coefficient estimates
  s2 <- fit$s ^ 2
  V <- try(chol2inv(chol(crossprod(X)))*s2)
  
  # Wald test (H0: Lb = 0)
  if (any(class(V) == "try-error")) {
    
    p <- NA
    
  } else {
    
    Lb <- L %*% b
    
    X2 <- t(Lb) %*% chol2inv(chol(L %*% V %*% t(L))) %*% Lb
    
    p <- pchisq(c(X2), df=nrow(L), lower.tail=FALSE)
    
  }
  
  return(p)
  
}

mean_impute <- function(x, f) {
  
  stopifnot(is.matrix(x))
  
  f <- as.character(f)
  
  for (j in 1:ncol(x)) {
    
    means <- tapply(x[, j], f, mean, na.rm=TRUE)
    
    idx <- which(is.na(x[, j]))
    
    x[idx, j] <- means[f[idx]]
    
  }
  
  return(x)
  
}

#--------------------------------------------------------
# Data
#--------------------------------------------------------
# Genomic relationships
GRM <- fread(GRM_file) %>%
  as.matrix() %>%
  to_PD()

rownames(GRM) <- fread(ID_file) %>%
  select(V1) %>%
  unlist() %>%
  unname()

colnames(GRM) <- rownames(GRM)

# Proteome
pangenes <- list.dirs(proteome_dir, full.names=FALSE, recursive=FALSE)

proteome <- fread(proteome_file) %>%
  filter(pan_gene_id %in% paste0("pan_gene_", .GlobalEnv$pangenes)) %>%
  column_to_rownames(var="pan_gene_id") %>%
  as.matrix()

proteome <- gsub("_T", "_P", proteome, perl=TRUE)

call_rate <- 1 - rowMeans(is.na(proteome))

proteome <- proteome[call_rate == 1, ]

rownames(proteome) <- sub("pan_gene_", "", rownames(proteome))
pangenes <- rownames(proteome)

# Proteome-wide relationships
PRM_list <- readRDS(PRM_file)
PRM <- PRM_list[[sim_metric]]

# Phenotypes
pheno <- fread(pheno_file) %>%
  filter(! duplicated(Geno_Code)) %>%
  column_to_rownames(var="Geno_Code") %>%
  select(all_of(traits)) %>%
  as.matrix()

# Subsetting by overlapping genotype IDs
geno_ID <- Reduce(intersect, list(rownames(GRM),
                                  rownames(pheno),
                                  colnames(proteome)))

GRM <- GRM[geno_ID, geno_ID]
GRM <- GRM/mean(diag(GRM))

PRM <- PRM[geno_ID, geno_ID]
PRM <- PRM/mean(diag(PRM))

proteome <- proteome[, geno_ID]
pheno <- pheno[geno_ID, ]

NAM_family <- substr(geno_ID, 1, 4)

gc()

# Protein PCs
PPC_list <- readRDS(PPC_file)
PPC_list <- PPC_list[names(PPC_list) != sim_metric]

#--------------------------------------------------------
# Proteome-wide association study
#--------------------------------------------------------
PWAS <- data.frame()

for (trait in traits) {
  
  cat("\n========================\n", trait, "\n========================\n")
  t0 <- Sys.time()
  
  y <- pheno[, trait]
  
  obs <- is.finite(y)
  N <- sum(obs)
  
  y.obs <- y[obs]
  
  GRM.obs <- GRM[obs, obs]
  PRM.obs <- PRM[obs, obs]
  
  Z.obs <- proteome[, obs]
  
  NAM_family.obs <- NAM_family[obs]
  
  if (include_family_structure) {
    X.obs <- model.matrix(~ NAM_family.obs)
  } else {
    X.obs <- matrix(1, nrow=N, ncol=1)
  }
  
  # Null model
  m0 <- try(greml(y=y.obs, X=X.obs, GRM=list(G=GRM.obs, P=PRM.obs)))
  
  if (class(m0) == "try-error") {
    
    print("Error fitting null model: skipping trait")
    
  } else {
    
    # Null variance matrix
    V <- GRM.obs*m0$theta["G"] + PRM.obs*m0$theta["P"] + diag(N)*m0$theta["E"]
    EVD <- eigen(V)
    
    # Whitening transformation
    V.inv.sqrt <- EVD$vectors %*% diag(1/sqrt(EVD$values)) %*% t(EVD$vectors)
    
    y.adj <- V.inv.sqrt %*% y.obs
    
    # Protein associations
    t0 <- Sys.time()
    out <- foreach(sim_metric=names(PPC_list), .combine=rbind) %do% {
      
      print(sim_metric)
      
      PPC <- PPC_list[[sim_metric]]
      PPC <- PPC[lengths(PPC) > 0]
      
      foreach(pangene=intersect(pangenes, names(PPC)), .combine=rbind) %:%
        foreach(k=k_list, .combine=rbind) %do% {
          
          # Adjusted predictors
          PPC_pangene <- permute(PPC[[pangene]], k)
          
          PPC.obs <- match_to_matrix(Z.obs[pangene, ], PPC_pangene)
          
          n_missing <- apply(PPC.obs, 1, function(x) any(is.na(x))) %>%
            sum()
          
          if (n_missing > 0) PPC.obs <- mean_impute(PPC.obs, NAM_family.obs)
          
          X.adj <- V.inv.sqrt %*% cbind(X.obs, PPC.obs)
          L <- diag(ncol(X.adj))[-(1:ncol(X.obs)), , drop=FALSE]
          
          # Wald test
          pval <- fastWald.test(X.adj, y.adj, L)
          
          data.frame(
            trait=trait,
            permutation=k,
            n=N,
            n_missing=n_missing,
            sim_metric=sim_metric,
            pangene=pangene,
            pval=pval
          )
        
      }
      
    }
    t1 <- Sys.time()
    print(t1 - t0)
    
    # Output
    PWAS <- rbind(PWAS, out)
    
    suffix <- ""
    if (length(k_list) == 1 & k_list[1] != 0) suffix <- paste0("_perm=", k_list[1])
    if (length(k_list) > 1) suffix <- paste0("_perm=", k_list[1], "-", k_list[length(k_list)])
    
    fwrite(out, paste0(output_dir, trait, suffix, ".csv"))
    
    saveRDS(PWAS, paste0(output_dir,"PWAS", suffix, ".rds"))
    
  }
  
}
