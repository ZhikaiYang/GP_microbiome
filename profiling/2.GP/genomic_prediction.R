setwd('/common/jyanglab/zhikaiyang/projects/GP_microbiome')
library(data.table)
phenos = fread("data/geno_trait.txt",data.table = F)
genos = fread("largedata/hmp321_282_agpv4_maf005_miss02_pruned_matrix.txt",data.table = F)

itraits = 40

gp = function(itrait = itraits, geno = genos, pheno = phenos){
  
  na.index <-  which(is.na(pheno[,itrait]))
  # length(na.index)
  pheno_i <- pheno[-na.index,c(1,itrait)]
  
  length(which(geno$genotype %in% pheno_i$genotype))
  length(which(pheno_i$genotype %in% geno$genotype))
  
  # merge genotype and phenotype
  pheno_geno = merge(pheno_i, geno, by.x = "genotype", by.y ="genotype")
  dim(pheno_geno)   
  
  # phenotypes 
  y <- pheno_geno[,2]
  y <- matrix(y, ncol=1)
  
  # markers 
  geno2 <- pheno_geno[,-c(1,2)] # 
  
  set.seed(1234)
  id = sort(sample(1:ncol(geno2),50000))
  geno2 = geno2[,id]
  
  # Missing marker imputation
  Z <- matrix(0, ncol=ncol(geno2), nrow=nrow(geno2))
  for (j in 1:ncol(geno2)){
    #cat("j = ", j, '\n')
    Z[,j] <- ifelse(is.na(geno2[,j]), mean(geno2[,j], na.rm=TRUE), geno2[,j])
  }
  
  # SNP Matrix standardization
  Zs <- scale(Z, center = TRUE, scale = TRUE)
  # dimensions 
  n <- nrow(Zs)
  m <- ncol(Zs)
  
  # Calcualte genomic relationship
  G <- tcrossprod(Zs) / ncol(Zs)  # G <- Zs %*% t(Zs) / ncol(Zs)
  G <- G + diag(n)*0.001
  
  #Fit GBLUP by using the mixed.solve function in the [rrBLUP]
  library(rrBLUP)
  fit <- mixed.solve(y = y, K=G)
  # additive genetic variance
  fit$Vu
  # residual variance
  fit$Ve
  # intercept 
  fit$beta
  # additive genetic values
  head(fit$u)
  tail(fit$u)
  # genomic h2
  fit$Vu / (fit$Vu + fit$Ve)
  # ratio of variance components 
  fit$Ve / fit$Vu
  #accuracy
  cor.test(y,fit$u)
  
  #Fit rrBLUP by using the mixed.solve function in the [rrBLUP]
  fit2 <- mixed.solve(y = y, Z=Zs)
  # marker additive genetic variance
  fit2$Vu
  # residual variance
  fit2$Ve
  # intercept 
  fit2$beta
  # marker additive genetic effects
  head(fit2$u)
  tail(fit2$u)
  # ratio of variance components 
  fit2$Ve / fit2$Vu
  
  # accuracy
  y.hat2 <- Zs %*% matrix(fit2$u)
  cor.test(y, y.hat2)
  
  
  
  # K-fold validation
  K = 10
  d = floor(n/K)
  set.seed(1234)
  i.mix = sample(1:n)
  folds = vector(mode="list",length=K)
  
  for (k in 1:K) {
    folds[[k]] = i.mix[((k-1)*d+1):(k*d)]
  }
  
  p_K10_rr <- rep(0.0, 10)
  
  for (k in 1:K) {
    cat("Fold",k,"\n")
    
    i.trn = unlist(folds[-k])
    i.tst = folds[[k]]
    
    
    y.trn = y[i.trn]    # training responses
    
    y.tst = y[i.tst]    # testing responses
    
    Zs.trn <- Zs[i.trn,]
    Zs.tst <- Zs[i.tst,]
    
    fit2.trn <- mixed.solve(y = y.trn, Z=Zs.trn)
    
    
    # prediction
    y.hat2 <- Zs.tst %*% matrix(fit2.trn$u)
    p_K10_rr[k] <- cor(y.hat2, y.tst)
    
  }
  
  
  p_K10 <- data.frame( trait = c(rep(colnames(pheno)[itrait],10)) , method = c(rep("rrBLUP", 10)),accuracy= c(p_K10_rr), fold = 1:10, N = nrow(y) )
  return(p_K10)
  
}

p_k10_normal = gp(itrait = itraits, geno = genos, pheno = phenos)

fwrite(p_k10_normal, "largedata/prediction_accuracy_normal.txt", sep = "\t", quote = F, row.names = F, col.names = T)


################################################################################################################


microbiomes = fread("data/blup_N_300_tax_groups.txt", data.table = F)
#principal components
pcs = fread("data/hmp321_282_agpv4_maf005_miss03_pruned.eigenvec",data.table = F)
pcs = pcs[,2:5]
colnames(pcs) = c("genotype", paste0(rep("pc",3),1:3))

gp_microbiome = function(itrait = itraits, N = 1, geno = genos, pheno = phenos, microbiome = microbiomes, pc = pcs){
  
  na.index <-  which(is.na(pheno[,itrait]))
  # length(na.index)
  pheno_i <- pheno[-na.index,c(1,itrait)]
  
  length(which(geno$genotype %in% pheno_i$genotype))
  length(which(pheno_i$genotype %in% geno$genotype))
  
  # merge genotype and phenotype
  pheno_geno = merge(pheno_i, geno, by.x = "genotype", by.y ="genotype")
  
  dim(pheno_geno)   
  
  if (N == 0) {
    microbiome = microbiome[,1:151]
  }else{
    microbiome = microbiome[,c(1,152:301)]
  }
  
  na.index <-  which(is.na(microbiome[,2]))
  # length(na.index)
  microbiome <- microbiome[-na.index,]
  
  
  pheno_geno_microbiome = merge(pheno_geno, microbiome, by.x = "genotype", by.y ="genotype")
  # phenotypes 
  y <- pheno_geno_microbiome[,2]
  y <- matrix(y, ncol=1)
  
  # markers 
  geno2 <- pheno_geno_microbiome[,3:ncol(pheno_geno)] # 
  
  set.seed(1234)
  id = sort(sample(1:ncol(geno2),50000))
  geno2 = geno2[,id]
  
  # Missing marker imputation
  Z <- matrix(0, ncol=ncol(geno2), nrow=nrow(geno2))
  for (j in 1:ncol(geno2)){
    #cat("j = ", j, '\n')
    Z[,j] <- ifelse(is.na(geno2[,j]), mean(geno2[,j], na.rm=TRUE), geno2[,j])
  }
  
  # SNP Matrix standardization
  Zs <- scale(Z, center = TRUE, scale = TRUE)
  # dimensions 
  n <- nrow(Zs)
  m <- ncol(Zs)
  
  
  #Fit rrBLUP by using the mixed.solve function in the [rrBLUP]
  fit2 <- mixed.solve(y = y, Z=Zs)
  # marker additive genetic variance
  fit2$Vu
  # residual variance
  fit2$Ve
  # intercept 
  fit2$beta
  # marker additive genetic effects
  head(fit2$u)
  tail(fit2$u)
  # ratio of variance components 
  fit2$Ve / fit2$Vu
  
  # accuracy
  y.hat2 <- Zs %*% matrix(fit2$u)
  cor.test(y, y.hat2)
  
  
  
  # K-fold validation
  K = 10
  d = floor(n/K)
  set.seed(1234)
  i.mix = sample(1:n)
  folds = vector(mode="list",length=K)
  
  for (k in 1:K) {
    folds[[k]] = i.mix[((k-1)*d+1):(k*d)]
  }
  
  p_K10_rr <- rep(0.0, 10)
  
  for (k in 1:K) {
    cat("Fold",k,"\n")
    
    i.trn = unlist(folds[-k])
    i.tst = folds[[k]]
    
    
    y.trn = y[i.trn]    # training responses
    
    y.tst = y[i.tst]    # testing responses
    
    Zs.trn <- Zs[i.trn,]
    Zs.tst <- Zs[i.tst,]
    
    fit2.trn <- mixed.solve(y = y.trn, Z=Zs.trn)
    
    
    # prediction
    y.hat2 <- Zs.tst %*% matrix(fit2.trn$u)
    p_K10_rr[k] <- cor(y.hat2, y.tst)
    
  }
  
  
  p_K10 <- data.frame( trait = c(rep(colnames(pheno)[itrait],10)) , method = c(rep("rrBLUP_less", 10)),accuracy= c(p_K10_rr), fold = 1:10 , N = nrow(y))
  
  #GP with microbiome and PCs
  
  microbiome2 = pheno_geno_microbiome[,(ncol(pheno_geno)+1):ncol(pheno_geno_microbiome)]
  
  pheno_pc = merge(pheno_geno_microbiome[,1:3], pc, by.x="genotype", by.y = "genotype")
  
  pc2 = pheno_pc[,4:ncol(pheno_pc)]
  
  pc2 = as.matrix(pc2)
  
  # Missing microbiome imputation
  Z_microbiome <- matrix(0, ncol=ncol(microbiome2), nrow=nrow(microbiome2))
  for (j in 1:ncol(microbiome2)){
    #cat("j = ", j, '\n')
    Z_microbiome[,j] <- ifelse(is.na(microbiome2[,j]), mean(microbiome2[,j], na.rm=TRUE), microbiome2[,j])
  }
  
  
  
  
  #Fit rrBLUP with PCs and microbiome
  fit2_microbiome <- mixed.solve(y = y, Z=Z_microbiome, X=pc2)
  # marker additive genetic variance
  fit2_microbiome$Vu
  # residual variance
  fit2_microbiome$Ve
  # intercept 
  fit2_microbiome$beta
  # marker additive genetic effects
  head(fit2_microbiome$u)
  tail(fit2_microbiome$u)
  # ratio of variance components 
  fit2_microbiome$Ve / fit2_microbiome$Vu
  
  # accuracy
  y.hat2 <- pc2 %*% matrix(fit2_microbiome$beta) + Z_microbiome %*% matrix(fit2_microbiome$u)
  cor.test(y.hat2, y)
  
  
  # K-fold validation
  K = 10
  d = floor(n/K)
  set.seed(1234)
  i.mix = sample(1:n)
  folds = vector(mode="list",length=K)
  
  for (k in 1:K) {
    folds[[k]] = i.mix[((k-1)*d+1):(k*d)]
  }
  
  p_K10_rr_microbiome <- rep(0.0, 10)
  
  for (k in 1:K) {
    cat("Fold",k,"\n")
    
    i.trn = unlist(folds[-k])
    i.tst = folds[[k]]
    
    
    y.trn = y[i.trn]    # training responses
    
    y.tst = y[i.tst]    # testing responses
    
    Z_microbiome.trn <- Z_microbiome[i.trn,]
    Z_microbiome.tst <- Z_microbiome[i.tst,]
    
    pc2.trn = pc2[i.trn,]
    pc2.tst = pc2[i.tst,]
    
    fit2_microbiome.trn <- mixed.solve(y = y.trn, Z=Z_microbiome.trn, X=pc2.trn)
    
    
    # prediction
    y.hat2 <- pc2.tst %*% matrix(fit2_microbiome.trn$beta) + Z_microbiome.tst %*% matrix(fit2_microbiome.trn$u)
    p_K10_rr_microbiome[k] <- cor(y.hat2, y.tst)
    
  }
  
  p_K10_microbiome <- data.frame( trait = c(rep(colnames(pheno)[itrait],10)) , method = c(rep("rrBLUP_less_microbiome", 10)),accuracy= c(p_K10_rr_microbiome), fold = 1:10, N = nrow(y) )
  
  
  
  #Fit rrBLUP with Geno and Microbiome
  Z_geno_microbiome = cbind(Zs,Z_microbiome)
  K_geno_microbiome = diag(c(rep(fit2$Vu,ncol(Zs)),rep(fit2_microbiome$Vu,ncol(Z_microbiome))))
  fit2_geno_microbiome <- mixed.solve(y = y, Z=Z_geno_microbiome, K=K_geno_microbiome)
  
  # accuracy
  y.hat2 <- Z_geno_microbiome %*% matrix(fit2_geno_microbiome$u)
  cor.test(y.hat2, y)
  
  # K-fold validation
  K = 10
  d = floor(n/K)
  set.seed(1234)
  i.mix = sample(1:n)
  folds = vector(mode="list",length=K)
  
  for (k in 1:K) {
    folds[[k]] = i.mix[((k-1)*d+1):(k*d)]
  }
  
  p_K10_rr_geno_microbiome <- rep(0.0, 10)
  
  for (k in 1:K) {
    cat("Fold",k,"\n")
    
    i.trn = unlist(folds[-k])
    i.tst = folds[[k]]
    
    
    y.trn = y[i.trn]    # training responses
    
    y.tst = y[i.tst]    # testing responses
    
    Z_geno_microbiome.trn <- Z_geno_microbiome[i.trn,]
    Z_geno_microbiome.tst <- Z_geno_microbiome[i.tst,]
    
    fit2_geno_microbiome.trn <- mixed.solve(y = y.trn, Z=Z_geno_microbiome.trn, K=K_geno_microbiome)
    
    
    # prediction
    y.hat2 <- Z_geno_microbiome.tst %*% matrix(fit2_geno_microbiome.trn$u)
    p_K10_rr_geno_microbiome[k] <- cor(y.hat2, y.tst)
    
  }
  
  
  p_K10_geno_microbiome <- data.frame( trait = c(rep(colnames(pheno)[itrait],10)) , method = c(rep("rrBLUP_less_geno_microbiome", 10)),accuracy= c(p_K10_rr_geno_microbiome), fold = 1:10, N = nrow(y) )
  
  p_K10_all = rbind(p_K10, p_K10_microbiome, p_K10_geno_microbiome)
  
  return(p_K10_all)
  
}

p_k10_all = gp_microbiome(itrait = itraits, N = 1, geno = genos, pheno = phenos, microbiome = microbiomes, pc = pcs)

fwrite(p_k10_all, "largedata/prediction_accuracy_all_again.txt", sep = "\t", quote = F, row.names = F, col.names = T)



