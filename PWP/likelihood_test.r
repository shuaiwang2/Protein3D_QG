
results <- data.frame(lambda = numeric(nrow(llik)), p_value = numeric(nrow(llik)))
rownames(results) <- rownames(llik)

for (i in 1:nrow(llik)) {

  logL0_val <- llik$structure_esm[i]
  logL1_val <- llik$structure[i]
  
  lambda <- -2 * (logL0_val - logL1_val)
  p_value <- pchisq(lambda, df = 1, lower.tail = FALSE)
  

  results$lambda[i] <- lambda
  results$p_value[i] <- p_value
}

results$p_value <0.01


results <- data.frame(lambda = numeric(nrow(llik)), p_value = numeric(nrow(llik)))
rownames(results) <- rownames(llik)

for (i in 1:nrow(llik)) {
  
  logL0_val <- llik$sequencemafft[i]
  logL1_val <- llik$sequencemafft_structure[i]
  
  lambda <- -2 * (logL0_val - logL1_val)
  p_value <- pchisq(lambda, df = 1, lower.tail = FALSE)
  
  
  results$lambda[i] <- lambda
  results$p_value[i] <- p_value
}

results$p_value <0.05 #yes


for (i in 1:nrow(llik)) {
  
  logL0_val <- llik$sequencemuscle[i]
  logL1_val <- llik$sequencemuscle_structure[i]
  
  lambda <- -2 * (logL0_val - logL1_val)
  p_value <- pchisq(lambda, df = 1, lower.tail = FALSE)
  
  
  results$lambda[i] <- lambda
  results$p_value[i] <- p_value
}

results$p_value <0.05 #yes

for (i in 1:nrow(llik)) {
  
  logL0_val <- llik$sequencetcoffee[i]
  logL1_val <- llik$sequencetcoffee_structure[i]
  
  lambda <- -2 * (logL0_val - logL1_val)
  p_value <- pchisq(lambda, df = 1, lower.tail = FALSE)
  
  
  results$lambda[i] <- lambda
  results$p_value[i] <- p_value
}

results$p_value <0.05 #yes


for (i in 1:nrow(llik)) {
  
  logL0_val <- llik$haplotype[i]
  logL1_val <- llik$haplotype_structure[i]
  
  lambda <- -2 * (logL0_val - logL1_val)
  p_value <- pchisq(lambda, df = 1, lower.tail = FALSE)
  
  
  results$lambda[i] <- lambda
  results$p_value[i] <- p_value
}

results[results$p_value >0.05,] #yes
