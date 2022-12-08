setwd('/common/jyanglab/zhikaiyang/projects/GP_microbiome')

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

id_seed <- as.numeric(as.character(args[1]))
id_date <- as.numeric(as.character(args[2]))
id_trait <- as.numeric(as.character(args[3]))

print(id_seed)
print(id_date)
print(id_trait)

library(data.table)
library(dplyr)
library(tidyr)
library(rrBLUP)

sinfo_microbiome_geno = fread("largedata/Zhikai/sinfo_micobiome_geno.txt",header = T, data.table = F)


# information from column 1 to 6, and column 3633

# ceil_id : used for merging with corresponding phenotype
# row : the row number of the observation for one genotype
# nitroge : level of nitrogen fertilizer treatment
# block : level of block factor
# sp and spb :  split plot and split plot block, which are used to account for location variation, I didn't used the information in analysis
# genotype (column 3633) : genotype name, used for merging with pcs



# microbiome from column 7 to 3632

# asv_###### : microbe id and its log(relative abundance (from 0  to 1) + 0.001 (a small pseudocount, so that log transformation works))



# genome from 3634 to 53633

# 1_25630 to 10_149467739 : snps id and their vaules (0-2)



#phenotype
#vegetation index
vis = fread("data/ppj220030-sup-0003-tables2.csv",header = T, data.table = F)
# information from column 1 to 6

# Row : the row number of the observation for one genotype
# Row.Numbers, Genotype, Pedigree, Treatment : won't be used 
# date: the date of data collected

# vegetation indexes from column 7 to 14



# spread the combined data for trait "canopy coverage" with measurements from date 6-July to 5-Sep
vis_canopy = vis[,c(1:6,id_trait)]
table(vis_canopy$date)
vis_canopy = spread(vis_canopy, key = "date", value = colnames(vis)[id_trait]) 

# ceil_id : used for merging with corresponding phenotype
vis_canopy$ceil_id = paste("c",ceiling(vis_canopy$Row/2),sep = "_")
vis_canopy = vis_canopy[,c(18, 14, 6:13, 15:17)]
vis_canopy = vis_canopy[,c(1,3,7,13)]
sinfo_pheno_microbiome_geno = merge(vis_canopy[,c(1,id_date)], sinfo_microbiome_geno, by.x = "ceil_id", by.y = "ceil_id")


#principal components
pcs = fread("data/hmp321_282_agpv4_maf005_miss03_pruned.eigenvec")
pcs = pcs[,2:5]
colnames(pcs) = c("genotype", paste0(rep("pc",3),1:3))

sinfo_pheno_microbiome_geno_pc = merge(sinfo_pheno_microbiome_geno, pcs, by.x = "genotype", by.y = "genotype")

id_na = which(is.na(sinfo_pheno_microbiome_geno_pc[, which(colnames(sinfo_pheno_microbiome_geno_pc) == colnames(vis_canopy)[id_date])]))
if (length(id_na) != 0) {
  sinfo_pheno_microbiome_geno_pc = sinfo_pheno_microbiome_geno_pc[-id_na, ]
}


#####################################################################

  
  y = matrix(sinfo_pheno_microbiome_geno_pc[,which(colnames(sinfo_pheno_microbiome_geno_pc) == colnames(vis_canopy)[id_date])],ncol = 1)
  
  
  #rrBLUP
  
  X = cbind(rep(1,nrow(y)), 
            ifelse(sinfo_pheno_microbiome_geno_pc$nitrogen == "HN", 1, 0), 
            ifelse(sinfo_pheno_microbiome_geno_pc$block == "N", 1, 0))
  
  id_1st_snp = which(colnames(sinfo_pheno_microbiome_geno_pc) == "1-25630")  #  3635th column is the 1st SNP
  
  Z = as.matrix(sinfo_pheno_microbiome_geno_pc[,id_1st_snp:(id_1st_snp+49999)])
  
  # SNP Matrix standardization
  Zs <- scale(Z, center = TRUE, scale = TRUE)
  # dimensions 
  n <- nrow(Zs)
  m <- ncol(Zs)
  
  
  #Fit rrBLUP by using the mixed.solve function in the [rrBLUP]
  fit2 <- mixed.solve(y = y, Z=Zs, X=X)
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
  y.hat2 <- X %*% matrix(fit2$beta)+ Zs %*% matrix(fit2$u)
  cor.test(y, y.hat2)
  
  
  # K-fold validation
  K = 5
  d = floor(n/K)
  set.seed(id_seed)
  i.mix = sample(1:n)
  folds = vector(mode="list",length=K)
  
  for (k in 1:K) {
    folds[[k]] = i.mix[((k-1)*d+1):(k*d)]
  }
  
  p_K10_rr <- rep(0.0, 5)
  
  for (k in 1:K) {
    cat("Fold",k,"\n")
    
    i.trn = unlist(folds[-k])
    i.tst = folds[[k]]
    
    
    y.trn = y[i.trn]    # training responses
    
    y.tst = y[i.tst]    # testing responses
    
    Zs.trn <- Zs[i.trn,]
    Zs.tst <- Zs[i.tst,]
    
    X.trn <- X[i.trn,]
    X.tst <- X[i.tst,]
    
    fit2.trn <- mixed.solve(y = y.trn, Z=Zs.trn, X=X.trn)
    
    
    # prediction
    y.hat2 <- X.tst %*% matrix(fit2.trn$beta) + Zs.tst %*% matrix(fit2.trn$u)
    p_K10_rr[k] <- cor(y.hat2, y.tst)
    
  }
  
  
  p_K10 <- data.frame(seed = c(rep(id_seed, 5)), date = c(rep(colnames(sinfo_pheno_microbiome_geno_pc)[3],5)) , method = c(rep("rrBLUP", 5)),accuracy= c(p_K10_rr), fold = 1:5 , N = nrow(y))
  
  
  #rrBLUP with microbiome and PCs
  X_microbiome = cbind(rep(1,nrow(y)), 
                       ifelse(sinfo_pheno_microbiome_geno_pc$nitrogen == "HN", 1, 0), 
                       ifelse(sinfo_pheno_microbiome_geno_pc$block == "N", 1, 0),
                       sinfo_pheno_microbiome_geno_pc$pc1,
                       sinfo_pheno_microbiome_geno_pc$pc2,
                       sinfo_pheno_microbiome_geno_pc$pc3)
  
  id_1st_asv = which(colnames(sinfo_pheno_microbiome_geno_pc) == "asv_000048")
  
  Z_microbiome = as.matrix(sinfo_pheno_microbiome_geno_pc[,id_1st_asv:(id_1st_asv+3625)])
  
  #Fit rrBLUP with PCs and microbiome
  fit2_microbiome <- mixed.solve(y = y, Z=Z_microbiome, X=X_microbiome)
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
  y.hat2 <- X_microbiome %*% matrix(fit2_microbiome$beta) + Z_microbiome %*% matrix(fit2_microbiome$u)
  cor.test(y.hat2, y)
  
  
  # K-fold validation
  K = 5
  d = floor(n/K)
  set.seed(id_seed)
  i.mix = sample(1:n)
  folds = vector(mode="list",length=K)
  
  for (k in 1:K) {
    folds[[k]] = i.mix[((k-1)*d+1):(k*d)]
  }
  
  p_K10_rr_microbiome <- rep(0.0, 5)
  
  for (k in 1:K) {
    cat("Fold",k,"\n")
    
    i.trn = unlist(folds[-k])
    i.tst = folds[[k]]
    
    
    y.trn = y[i.trn]    # training responses
    
    y.tst = y[i.tst]    # testing responses
    
    Z_microbiome.trn <- Z_microbiome[i.trn,]
    Z_microbiome.tst <- Z_microbiome[i.tst,]
    
    X_microbiome.trn = X_microbiome[i.trn,]
    X_microbiome.tst = X_microbiome[i.tst,]
    
    fit2_microbiome.trn <- mixed.solve(y = y.trn, Z=Z_microbiome.trn, X=X_microbiome.trn)
    
    
    # prediction
    y.hat2 <- X_microbiome.tst %*% matrix(fit2_microbiome.trn$beta) + Z_microbiome.tst %*% matrix(fit2_microbiome.trn$u)
    p_K10_rr_microbiome[k] <- cor(y.hat2, y.tst)
    
  }
  
  p_K10_microbiome <- data.frame( seed = c(rep(id_seed, 5)), date = c(rep(colnames(sinfo_pheno_microbiome_geno_pc)[3],5)) , method = c(rep("rrBLUP_PC_microbiome", 5)),accuracy= c(p_K10_rr_microbiome), fold = 1:5, N = nrow(y) )
  
  
  
  
  #Fit rrBLUP with Geno and Microbiome
  X_geno_microbiome = cbind(rep(1,nrow(y)), 
                            ifelse(sinfo_pheno_microbiome_geno_pc$nitrogen == "HN", 1, 0), 
                            ifelse(sinfo_pheno_microbiome_geno_pc$block == "N", 1, 0))
  
  Z_geno_microbiome = cbind(Zs,Z_microbiome)
  
  
  K_geno_microbiome = diag(c(rep(fit2$Vu,ncol(Zs)),rep(fit2_microbiome$Vu,ncol(Z_microbiome))))
  
  rm(geno_names)
  rm(geno)
  rm(sinfo_microbiome)
  rm(sinfo_microbiome_avg)
  rm(sinfo_pheno_microbiome_geno)
  rm(X)
  rm(X_microbiome)
  rm(Z_microbiome)
  rm(Z)
  
  fit2_geno_microbiome <- mixed.solve(y = y, Z=Z_geno_microbiome, K=K_geno_microbiome, X=X_geno_microbiome)
  
  # accuracy
  y.hat2 <- X_geno_microbiome %*% matrix(fit2_geno_microbiome$beta) + Z_geno_microbiome %*% matrix(fit2_geno_microbiome$u)
  cor.test(y.hat2, y)
  
  # K-fold validation
  K = 5
  d = floor(n/K)
  set.seed(id_seed)
  i.mix = sample(1:n)
  folds = vector(mode="list",length=K)
  
  for (k in 1:K) {
    folds[[k]] = i.mix[((k-1)*d+1):(k*d)]
  }
  
  p_K10_rr_geno_microbiome <- rep(0.0, 5)
  
  for (k in 1:K) {
    cat("Fold",k,"\n")
    
    i.trn = unlist(folds[-k])
    i.tst = folds[[k]]
    
    
    y.trn = y[i.trn]    # training responses
    
    y.tst = y[i.tst]    # testing responses
    
    Z_geno_microbiome.trn <- Z_geno_microbiome[i.trn,]
    Z_geno_microbiome.tst <- Z_geno_microbiome[i.tst,]
    
    X_geno_microbiome.trn <- X_geno_microbiome[i.trn,]
    X_geno_microbiome.tst <- X_geno_microbiome[i.tst,]
    
    
    fit2_geno_microbiome.trn <- mixed.solve(y = y.trn, Z=Z_geno_microbiome.trn, K=K_geno_microbiome, X = X_geno_microbiome.trn)
    
    
    # prediction
    y.hat2 <- X_geno_microbiome.tst %*% matrix(fit2_geno_microbiome.trn$beta) + Z_geno_microbiome.tst %*% matrix(fit2_geno_microbiome.trn$u)
    p_K10_rr_geno_microbiome[k] <- cor(y.hat2, y.tst)
    
  }
  
  
  p_K10_geno_microbiome <- data.frame( seed = c(rep(id_seed, 5)), date = c(rep(colnames(sinfo_pheno_microbiome_geno_pc)[3],5)) , method = c(rep("rrBLUP_geno_microbiome", 5)),accuracy= c(p_K10_rr_geno_microbiome), fold = 1:5, N = nrow(y) )
  
  p_K10_all = rbind(p_K10, p_K10_microbiome, p_K10_geno_microbiome)
  
  outname = paste0("largedata/outputs/prediction_accuracy_seed_", id_seed, "_trait_",colnames(vis)[id_trait], "_date_",colnames(vis_canopy)[id_date],"_jyang.txt")
  fwrite(p_K10_all, outname, sep = "\t", quote = F, row.names = F, col.names = T)
  rm(sinfo_pheno_microbiome_geno_pc)






