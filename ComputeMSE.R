# input path 
cor.name <- 'pearson'
input.path <- paste0('Your working directory', cor.name)

# parameters
n.fold <- 100
n.pred <- 100
n.gene <- 347
n.pred.gene <- 35

# compute mse
compute.mse.pred2pred <- function(original.cor, prediction.cor){
  error <- (original.cor - prediction.cor)^2
  sse <- sum(error[1:n.pred.gene, 1:n.pred.gene])/2
  mse <- sse/((n.pred.gene^2 - n.pred.gene)/2)
  return(mse)
}

# create mse matrix
mse.pred2pred <- matrix(data = NA, nrow = n.fold, ncol = 4)
colnames(mse.pred2pred) <- c('SP', 'MP_Mean', 'MP_Median', 'MP_Fisher')

for (k in 1:n.fold){
  # load original correlation  
  file <- paste0(input.path, '/Fold_', k, '_Original_', cor.name, '.txt')
  original.cor <- read.table(file, sep = '\t')

  # load single imputation correlation
  file <- paste0(input.path, '/Fold_', k, '_SI_', cor.name, '.txt')
  si.cor <- read.table(file, sep = '\t')
  
  # load multiple imputation correlations
  file <- paste0(input.path, '/Fold_', k, '_MI_', cor.name, '_mean.txt')
  mi.cor.mean <- read.table(file, sep = '\t')

  file <- paste0(input.path, '/Fold_', k, '_MI_', cor.name, '_median.txt')
  mi.cor.median <- read.table(file, sep = '\t')
  
  file <- paste0(input.path, '/Fold_', k, '_MI_', cor.name, '_fisher.txt')
  mi.cor.fisher <- read.table(file, sep = '\t')
  
  si.mse.pred2pred <- compute.mse.pred2pred(original.cor, si.cor)
  mi.mean.mse.pred2pred <- compute.mse.pred2pred(original.cor, mi.cor.mean)
  mi.median.mse.pred2pred <- compute.mse.pred2pred(original.cor, mi.cor.median)
  mi.fisher.mse.pred2pred <- compute.mse.pred2pred(original.cor, mi.cor.fisher)
  mse.pred2pred[k,] <- c(si.mse.pred2pred, mi.mean.mse.pred2pred, 
               mi.median.mse.pred2pred, mi.fisher.mse.pred2pred)
}

# output path
output.path <- 'Your working directory'
file <- paste0(output.path, 'MSE_pred2pred_', cor.name, '.txt')
write.table(mse.pred2pred, file, quote = F, sep = '\t', row.names = F, col.names = T)
