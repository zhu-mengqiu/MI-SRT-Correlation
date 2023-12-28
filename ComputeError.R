# input path 
cor.name <- 'pearson'
input.path <- paste0('Your working directory', cor.name)

# parameters
n.fold <- 100
n.pred <- 100
n.gene <- 347
n.pred.gene <- 35

# compute classification error
compute.error.rate <- function(original.net, prediction.net){
  error <- original.net != prediction.net 
  
  pred2pred <- sum(error[1:n.pred.gene, 1:n.pred.gene])/2
  pred2pred.num <- (n.pred.gene^2 - n.pred.gene)/2
  
  error.rate.pred2pred <- pred2pred / pred2pred.num
  
  return(error.rate.pred2pred)
}

error.pred2pred <- matrix(data = NA, nrow = n.fold, ncol = 2)
colnames(error.pred2pred) <- c('SI', 'MI_Median')

for (k in 1:n.fold){
  # load original network
  file <- paste0(input.path, '/Fold_', k, '_Original_', cor.name, '_n_b', '.txt')
  original.b <- read.table(file, sep = '\t')
  
  # load single imputation network
  file <- paste0(input.path, '/Fold_', k, '_SI_', cor.name, '_n_b', '.txt')
  si.b <- read.table(file, sep = '\t')
  
  # load multiple imputation network
  file <- paste0(input.path, '/Fold_', k, '_MI_', cor.name, '_n_b_mode','.txt')
  mi.b <- read.table(file, sep = '\t')
  
  error.rate.si <- compute.error.rate(original.b, si.b)
  error.pred2pred[k, 1] <- error.rate.si
  
  error.rate.mi <- compute.error.rate(original.b, mi.b)
  error.pred2pred[k, 2] <- error.rate.mi
}

# output path
output.path <- paste0("E:/seqFISH_10X/")

file <- paste0(output.path, 'Error_pred2pred_b_', cor.name, '.txt')
write.table(error.pred2pred, file, quote = F, sep = '\t', row.names = F, col.names = T)