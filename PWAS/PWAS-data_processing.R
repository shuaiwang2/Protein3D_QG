# PWAS: data processing
# Objective : Make input for PWAS
# Created by: Guillaume Ramstein (ramstein@qgg.au.dk)
# Created on: 12/12/2023

#--------------------------------------------------------
# Script parameters
#--------------------------------------------------------
overwrite <- TRUE

# Working directory
wd <- "/home/song/wangs/pwas_pwp"
dir.create(wd, showWarnings = FALSE)
setwd(wd)

# Input
GRM_file <- "snp.hBN.kinf"
ID_file <- "snp.tfam"

proteome_dir <- "/home/song/wangs/pwas_pwp/seq_structure/"

#proteome_file <- paste0(proteome_dir, "mergedPanGeneTables")
proteome_file <- "mergedPanGeneTables"

#pheno_file <- paste0(proteome_dir, "all_NAM_phenos.csv")
pheno_file <-"all_NAM_phenos.csv"
sim_metrics <- c(
                 "usalignMatrix_TMScore",
                 "mafftSequenceSimilarityMatrix",
                 "tcoffeeSequenceSimilarityMatrix",
                 "muscleSequenceSimilarityMatrix",
                 "haplotype"
                ) # Need to add ESM score similarity matrices

trace_threshold <- 1/25
lambda_threshold <- 0.01

# Output
QC_file <- "QC/QC.rds"
PRM_file <- "QC/PRM_list.rds"
PPC_file <- "QC/PPC_list.rds"

kinship_file <- "PWP/PRM.rds"


#--------------------------------------------------------
# Functions
#--------------------------------------------------------
library(tidyverse)
library(foreach)
library(data.table)
library(Matrix)

library(qgg)

to_PD <- function(x) {
  
  as.matrix(nearPD(x)$mat)
  
}

import_sim <- function(fn) {
  
  sim <- fread(fn)
  
  if (nrow(sim) > 0) {
    
    sim <- column_to_rownames(sim, var="V1") %>%
      as.matrix()
    
    colnames(sim) <- rownames(sim)
    
  } else {
    
    sim <- NULL
    
  }
  
  return(sim)
  
}

# cluster <- function(x, k_seq=1:10) {
#   
#   BIC <- sapply(k_seq, function(k) {
#     
#     fit <- kmeans(x, k, iter.max=100, nstart=10)
#     
#     n <- length(fit$cluster)
#     m <- ncol(fit$centers)
#     D <- fit$tot.withinss
#     
#     D + 2*m*k
#     
#   })
#   
#   kmeans(x, k_seq[which.min(BIC)], iter.max=100, nstart=10)
#   
# }

match_to_matrix <- function(id_vec, sim, sep=";") {
  
  sim_entries <- rownames(sim)
  
  id_list <- lapply(strsplit(id_vec, sep), function(id) {
    
    id[id %in% sim_entries]
    
  })
  
  stopifnot("All protein IDs must be present in similarity matrix"=all(lengths(id_list) > 0))
  
  id_DF <- suppressWarnings(do.call(rbind, id_list)) %>%
    as.data.frame()
  
  out <- foreach(id=id_DF, .combine=`+`) %do% {
    
    sim[match(id, rownames(sim)), match(id, colnames(sim))]
    
  } / ncol(id_DF)
  
  rownames(out) <- colnames(out) <- names(id_vec)
  
  return(out)
  
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

# Phenotypes
pheno <- fread(pheno_file) %>%
  filter(! duplicated(Geno_Code)) %>%
  column_to_rownames(var="Geno_Code") %>%
  as.matrix()

# Subsetting by overlapping genotype IDs
geno_ID <- Reduce(intersect, list(rownames(GRM),
                                  rownames(pheno),
                                  colnames(proteome)))

proteome <- proteome[, geno_ID]

gc()

#--------------------------------------------------------
# Analysis
#--------------------------------------------------------
# Filtering pangenes
if (overwrite || !(file.exists(QC_file))) {
  
  QC <- foreach(pangene=pangenes, .combine=rbind) %do% {
    
    pangene_dir <- paste0(proteome_dir, pangene, "/")
    
    x <- proteome[paste0("pan_gene_", pangene), ] %>%
      strsplit(";") %>%
      unlist() %>%
      unique()
    
    foreach(sim_metric=setdiff(sim_metrics, "haplotype"), .combine=rbind) %do% {
      
      fn <- paste0(pangene_dir, sim_metric)
      
      if (file.exists(fn)) {
        
        sim <- import_sim(fn)
        
        if (is.null(sim)) {
          
          out <- data.frame(pangene=pangene,
                            metric=sim_metric,
                            trace=NA,
                            lambda_max=NA,
                            lambda_min=NA,
                            n=NA,
                            match_rate=NA,
                            call_rate=NA
          )
          
        } else {
          
          N <- nrow(sim)
          H <- diag(N) - matrix(1/N, nrow=N, ncol=N)
          S <- H %*% sim %*% H
          
          lambda <- eigen((S + t(S))/2)$values
          
          out <- data.frame(pangene=pangene,
                            metric=sim_metric,
                            trace=sum(lambda),
                            lambda_max=max(lambda),
                            lambda_min=min(lambda),
                            n=N,
                            match_rate=mean(sub("[.]pdb", "", rownames(sim)) %in% x),
                            call_rate=mean(x %in% sub("[.]pdb", "", rownames(sim)))
          )
          
        }
        
      } else {
        
        out <- data.frame()
        
      }
      
      return(out)
      
    }
    
  }
  
  saveRDS(QC, QC_file)
  
} else {
  
  QC <- readRDS(QC_file)
  
}

# Selected pangenes
pangene_list <- foreach(sim_metric=setdiff(sim_metrics, "haplotype")) %do% {
  
  filter(QC, metric == sim_metric &
           match_rate == 1 & 
           n == 26) %>%
    select(pangene) %>%
    unique() %>%
    unlist()
  
}

selected_pangenes <- Reduce(intersect, pangene_list) %>%
  intersect(unique(QC$pangene[QC$trace > trace_threshold]))

empty_list <- vector("list", length(selected_pangenes)) %>%
  setNames(selected_pangenes)

# List of input matrices at selected pangenes
N <- 26
H <- diag(N) - matrix(1/N, nrow=N, ncol=N)

PRM_list <- list()
PPC_list <- list()

for (sim_metric in sim_metrics) {
  
  PRM_list[[sim_metric]] <- empty_list
  PPC_list[[sim_metric]] <- empty_list
  
  for (pangene in selected_pangenes) {
    
    if (sim_metric == "haplotype") {
      
      fn <- paste0(proteome_dir, pangene, "/usalignMatrix_TMScore")
      
      protein_ID <- sub("[.]pdb", "", rownames(import_sim(fn)))
      
      sim <- diag(length(protein_ID))
      
    } else {
      
      fn <- paste0(proteome_dir, pangene, "/", sim_metric)
      
      sim <- import_sim(fn)
      
      protein_ID <- sub("[.]pdb", "", rownames(sim))
      
    }
    
    # PRM
    S <- H %*% sim %*% H
    S <- (S + t(S))/2
    rownames(S) <- colnames(S) <- protein_ID
    
    PRM_list[[sim_metric]][[pangene]] <- S
    
    # PPC
    EVD <- eigen(S)
    
    idx <- which(EVD$values > lambda_threshold*trace_threshold)
    
    if (length(idx) > 0) {
      
      PC <- EVD$vectors[, idx, drop=FALSE] %*% diag(sqrt(EVD$values[idx]), nrow=length(idx))
      rownames(PC) <- protein_ID
      colnames(PC) <- paste0("PC", idx)
      
      PPC_list[[sim_metric]][[pangene]] <- PC
      
    }
    
  }
  
}

saveRDS(PRM_list, PRM_file)
saveRDS(PPC_list, PPC_file)

PPC_list <- list()
# Proteome-wide kinship matrices
kinship_list <- foreach(sim_metric=sim_metrics) %do% {
  
  print(sim_metric)

  sim <- PRM_list[[sim_metric]]

  #sim <- sim_list[[sim_metric]]
  
  out <- foreach(pangene=pangenes, .combine=`+`) %do% {
    
    ID <- proteome[paste0("pan_gene_", pangene), ]
    
    mat <- sim[[pangene]]
    
    PRM_pangene <- try(match_to_matrix(ID, mat), silent=TRUE)
    
    if (class(PRM_pangene)[1] != "try-error") {
      
      return(PRM_pangene)
      
    } else {
      
      return(0)
      
    }
    
  }
  
  return(to_PD(out))
  
}

names(kinship_list) <- sim_metrics

saveRDS(kinship_list, kinship_file)

gc()
