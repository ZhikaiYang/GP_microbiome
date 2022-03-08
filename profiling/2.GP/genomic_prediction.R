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
  
  
  p_K10 <- data.frame( trait = c(rep(colnames(pheno)[itrait],10)) , method = c(rep("rrBLUP", 10)),accuracy= c(p_K10_rr) )
  return(p_K10)
  
}

p_k10 = gp(itrait = itraits, geno = genos, pheno = phenos)

fwrite(p_K10, "largedata/prediction_accuracy.txt", sep = "\t", quote = F, row.names = F, col.names = T)



microbiome = 



