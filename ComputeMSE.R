# input path 
cor.name <- 'pearson'
input.path <- paste0('/CCAS/home/mzhu32/10X_seqFISH/', cor.name)

# parameters
n.fold <- 100
n.pred <- 100
n.gene <- 347
n.pred.gene <- 35

# compute mse undetected and undetected
compute.mse.pred2pred <- function(original.cor, prediction.cor){
  error <- (original.cor - prediction.cor)^2
  sse <- sum(error[1:n.pred.gene, 1:n.pred.gene])/2
  mse <- sse/((n.pred.gene^2 - n.pred.gene)/2)
  return(mse)
}

compute.mse.pred2origin <- function(original.cor, prediction.cor){
  error <- (original.cor - prediction.cor)^2
  sse <- sum(error[1:n.pred.gene, (n.pred.gene+1):n.gene])
  mse <- sse/(n.pred.gene*(n.gene-n.pred.gene))
  return(mse)
}


# compute mse detected and undetected

# create mse matrix
mse.pred2pred <- matrix(data = NA, nrow = n.fold, ncol = 4)
colnames(mse.pred2pred) <- c('SP', 'MP_Mean', 'MP_Median', 'MP_Fisher')
mse.pred2origin <- matrix(data = NA, nrow = n.fold, ncol = 4)
colnames(mse.pred2origin) <- c('SP', 'MP_Mean', 'MP_Median', 'MP_Fisher')

for (k in 1:n.fold){
  # load original correlation  
  file <- paste0(input.path, '/Fold_', k, '_Original_', cor.name, '.txt')
  original.cor <- read.table(file, sep = '\t')

  # load single prediction correlation
  file <- paste0(input.path, '/Fold_', k, '_SI_', cor.name, '.txt')
  si.cor <- read.table(file, sep = '\t')
  
  # load multiple prediction correlations
  mi.cor.list <- vector('list', n.pred)
  for (i in 1:n.pred){
    file <- paste0(input.path, '/Fold_', k, '_MI_', i,'_', cor.name, '.txt')
    mi.cor <- read.table(file, sep = '\t')
    mi.cor.list[[i]] <- as.matrix(mi.cor)
  }
  
  
  
  si.mse.pred2pred <- compute.mse.pred2pred(original.cor, si.cor)
  mi.mean.mse.pred2pred <- compute.mse.pred2pred(original.cor, mi.cor.mean)
  mi.median.mse.pred2pred <- compute.mse.pred2pred(original.cor, mi.cor.median)
  mi.fisher.mse.pred2pred <- compute.mse.pred2pred(original.cor, mi.cor.fisher)
  mse.pred2pred[k,] <- c(si.mse.pred2pred, mi.mean.mse.pred2pred, 
               mi.median.mse.pred2pred, mi.fisher.mse.pred2pred)
  
  si.mse.pred2origin <- compute.mse.pred2origin(original.cor, si.cor)
  mi.mean.mse.pred2origin <- compute.mse.pred2origin(original.cor, mi.cor.mean)
  mi.median.mse.pred2origin <- compute.mse.pred2origin(original.cor, mi.cor.median)
  mi.fisher.mse.pred2origin <- compute.mse.pred2origin(original.cor, mi.cor.fisher)
  mse.pred2origin[k,] <- c(si.mse.pred2origin, mi.mean.mse.pred2origin, 
                         mi.median.mse.pred2origin, mi.fisher.mse.pred2origin)
}

# output path
output.path <- '/CCAS/home/mzhu32/10X_seqFISH/'
file <- paste0(output.path, 'MSE_pred2pred_', cor.name, '.txt')
write.table(mse.pred2pred, file, quote = F, sep = '\t', row.names = F, col.names = T)
file <- paste0(output.path, 'MSE_pred2origin_', cor.name, '.txt')
write.table(mse.pred2origin, file, quote = F, sep = '\t', row.names = F, col.names = T)
