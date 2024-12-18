# GS validation
# Objective : Validate PWAS models by GS models
# Created by: Guillaume Ramstein (ramstein@qgg.au.dk)
# Created on: 9/27/2022

#--------------------------------------------------------
# Script parameters
#--------------------------------------------------------
overwrite <- FALSE

# Working directory
wd <- "/home/song/wangs/pwas_pwp"
#dir.create(wd, showWarnings = FALSE)
setwd(wd)

# Input
GRM_file <- "snp.hBN.kinf"
ID_file <- "snp.tfam"

#proteome_dir <- "proteome/"
#proteome_file <- paste0(proteome_dir, "mergedPanGeneTables")
proteome_file <- "mergedPanGeneTables"

pheno_file <- "all_NAM_phenos.csv"

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

sim_metrics <- c(
  "usalignMatrix_TMScore",
  "mafftSequenceSimilarityMatrix",
  "tcoffeeSequenceSimilarityMatrix",
  "muscleSequenceSimilarityMatrix",
  "haplotype"
)

metric_list <- list("baseline"=c("GRM"),
                    "sequencemafft"=c("GRM", "mafftSequenceSimilarityMatrix"),
                    "sequencetcoffee"=c("GRM", "tcoffeeSequenceSimilarityMatrix"),
                    "sequencemuscle"=c("GRM", "muscleSequenceSimilarityMatrix"),
                    "haplotype"=c("GRM", "haplotype"),
                    "structure"=c("GRM", "usalignMatrix_TMScore"),
                    "sequencemafft_structure"=c("GRM", "mafftSequenceSimilarityMatrix", "usalignMatrix_TMScore"),
                    "sequencetcoffee_structure"=c("GRM", "tcoffeeSequenceSimilarityMatrix", "usalignMatrix_TMScore"),
                    "sequencemuscle_structure"=c("GRM", "muscleSequenceSimilarityMatrix", "usalignMatrix_TMScore"),
                    "haplotype_structure"=c("GRM", "haplotype", "usalignMatrix_TMScore")
)

n_cores <- 12

sim_file <- paste0("QC/PRM_list.rds")

# Output
output_dir <- "PWP/"
dir.create(output_dir, showWarnings=FALSE)

PRM_file <- paste0(output_dir, "PRM.rds")
PWP_file <- paste0(output_dir, "PWP_analysis.csv")
LOFO_file <- paste0(output_dir, "LOFO_analysis.csv")


PRM_file <-"PWP/PRM.rds"
#--------------------------------------------------------
# Functions
#--------------------------------------------------------
library(tidyverse)
library(foreach)
library(data.table)
library(Matrix)

library(MM4LMM)
library(qgg)

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

# Proteome data
sim_list <- readRDS(sim_file)

pangenes <- Reduce(intersect, lapply(sim_list[sim_metrics], names))

proteome <- fread(proteome_file) %>%
  filter(pan_gene_id %in% paste0("pan_gene_", .GlobalEnv$pangenes)) %>%
  column_to_rownames(var="pan_gene_id") %>%
  as.matrix()

proteome <- gsub("_T", "_P", proteome, perl=TRUE)

call_rate <- 1 - rowMeans(is.na(proteome))

# Proteome-wide relationships
if (file.exists(PRM_file) && (! overwrite)) {
  
  PRM_list <- readRDS(PRM_file)
  
} else {
  
  PRM_list <- foreach(sim_metric=sim_metrics) %do% {
    
    print(sim_metric)
    
    sim <- sim_list[[sim_metric]]
    
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
  
  names(PRM_list) <- sim_metrics
  
  saveRDS(PRM_list, PRM_file)
  
  gc()
  
}

# Phenotypes
pheno <- fread(pheno_file) %>%
  filter(! duplicated(Geno_Code)) %>%
  column_to_rownames(var="Geno_Code") %>%
  select(all_of(traits)) %>%
  as.matrix()

write.csv(pheno,"pheno_32.csv",quote = F)
# Subsetting by overlapping genotype IDs
geno_ID <- Reduce(intersect, list(rownames(GRM),
                                  rownames(pheno),
                                  colnames(proteome)))

proteome <- proteome[, geno_ID]
pheno <- pheno[geno_ID, ]

NAM_families <- substr(geno_ID, 1, 4)

gc()

# Formatting relationship matrices
for (j in 1:length(PRM_list)) {
  
  PRM <- PRM_list[[j]][geno_ID, geno_ID]
  k <- mean(diag(PRM))
  
  PRM_list[[j]] <- PRM/k
  
}

GRM <- GRM[geno_ID, geno_ID]
GRM <- GRM/mean(diag(GRM))

RM_list <- append(list("GRM"=GRM), PRM_list)

#--------------------------------------------------------
# Variance partition
#--------------------------------------------------------
cat("\n************************\nVariance partition\n************************\n")

# Output data frame
pwp <- data.frame(
  trait=character(0),
  input=character(0),
  llik=numeric(0)
)

for (VC_name in c(unique(unlist(metric_list)), "E")) {
  
  pwp[, VC_name] <- numeric(0)
  
} 

# Analysis
for (trait in traits) {
  
  cat("\n========================\n", trait, "\n========================\n")
  t0 <- Sys.time()
  
  y <- pheno[, trait]
  
  obs <- is.finite(y)
  y.obs <- y[obs]
  X.obs <- matrix(1, nrow=sum(obs), ncol=1)
  
  RM_list.obs <- lapply(RM_list, `[`, obs, obs, drop=FALSE)
  
  # Genomic prediction
  out <- foreach(metric_name=names(metric_list), .combine=bind_rows) %do% {
    
    print(metric_name)
    
    # LMM
    m1 <- try(
      MMEst(Y=y.obs,
            X=X.obs,
            VarList=append(RM_list.obs[metric_list[[metric_name]]],
                           list("E"=diag(sum(obs)))),
            CritVar=1e-6,
            CritLogLik=1e-6,
            NbCores=n_cores)
    )
    
    if (class(m1) == "try-error") {
      
      var_info <- NULL
      
      ll1 <- NA
      
    } else {
      
      var_info <- data.frame(t(m1$X1$Sigma2))
      
      ll1 <- m1$X1$`LogLik (Reml)`
      
    }
    
    cbind(trait=trait,
          input=metric_name,
          var_info,
          llik=ll1)
    
  }
  
  t1 <- Sys.time()
  
  print(t1 - t0)
  
  # Output
  rownames(out) <- NULL
  
  pwp <- bind_rows(pwp, out)
  
  fwrite(pwp, PWP_file)
  
}

#--------------------------------------------------------
# Leave-one-family-out validation
#--------------------------------------------------------
cat("\n************************\nLOFO validation\n************************\n")

lofo <- data.frame()

for (trait in traits) {
  
  cat("\n========================\n", trait, "\n========================\n")
  t0 <- Sys.time()
  
  y <- pheno[, trait]
  
  obs <- is.finite(y)
  y.obs <- y[obs]
  X.obs <- matrix(1, nrow=sum(obs), ncol=1)
  
  RM_list.obs <- lapply(RM_list, `[`, obs, obs, drop=FALSE)
  
  # LOFO validation
  VS <- split(1:sum(obs), factor(NAM_families[obs]))
  
  out <- foreach(metric_name=names(metric_list), .combine=rbind) %do% {
    
    print(metric_name)
    
    # LMM
    m1 <- try(greml(y=y.obs, 
                    X=X.obs,
                    GRM=RM_list.obs[metric_list[[metric_name]]],
                    validate=VS,
                    ncores=n_cores))
    
    # Results
    if (class(m1) == "try-error") {
      
      val_info <- data.frame(NAM_family=NA,
                             Corr=NA,
                             R2=NA,
                             intercept=NA,
                             slope=NA,
                             MSE=NA)
      
    } else {
      
      val_info <- cbind(NAM_family=rownames(m1$accuracy),
                        foreach(DF=m1$validation, .combine=rbind) %do% {
                          
                          yobs <- DF[, "yobs"]
                          ypred <- DF[, "ypred"]
                          
                          fit <- lm(yobs ~ ypred)
                          
                          r <- cor(yobs, ypred, use="complete")
                          
                          data.frame(
                            Corr=r,
                            R2=r^2,
                            intercept=coef(fit)[1],
                            slope=coef(fit)[2],
                            MSE=mean((yobs-ypred)^2, na.rm=TRUE)
                          )
                          
                        })
      
    }
    
    cbind(trait=trait,
          input=metric_name,
          val_info)
    
  }
  
  t1 <- Sys.time()
  
  print(t1 - t0)
  
  # Output
  rownames(out) <- NULL
  
  lofo <- bind_rows(lofo, out)
  
  fwrite(lofo, LOFO_file)
  
}
