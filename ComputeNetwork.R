# input path
cor.name <- 'pearson'
input.path <- paste0('Your working directory', cor.name)

# parameters
n.fold <- 100
n.pred <- 100
n.gene <- 347
n.pred.gene <- 35

# function to compute connection
compute.connection <- function(cor, cutoff){
  connection <- 0
  
  if (cor >= cutoff){
    connection <- 1
  } elif (cor <= -cutoff){
    connection <- -1
  } else {
    connection <- 0
  }
  
  return (connection)
}

compute.network <- function(matrix.cor, cutoff){
  matrix.network <- apply(matrix.cor, c(1,2), compute.connection(x, cutoff))
  return (matrix.network)
}

compute.mode <- function(x){
  ux <- c(-1,0,1)
  ux[which.max(tabulate(match(x, ux)))]
}

cutoff <- 0.6

# compute network for original, SI, MI, and MI_combined
for (k in 1:n.fold){
  # load original correlation
  file <- paste0(input.path, '/Fold_', k, '_Original_', cor.name, '.txt')
  
  original.cor <- read.table(file, sep = '\t')
  original.network <- compute.network(original.cor, cutoff)
  
  file <- paste0(input.path, '/Fold_', k, '_Original_', cor.name, '_n_b', '.txt')
  write.table(original.network, file, sep = '\t', quote = F)
  
  # load single prediction correlation
  file <- paste0(input.path, '/Fold_', k, '_SI_', cor.name, '.txt')
  
  si.cor <- read.table(file, sep = '\t')
  si.network <- compute.network(si.cor, cutoff)
  
  file <- paste0(input.path, '/Fold_', k, '_SI_', cor.name, '_n_b', '.txt')
  write.table(si.network, file, sep = '\t', quote = F)
  
  # load multiple prediction correlations
  mi.network.list <- vector('list', n.pred)
  for (i in 1:n.pred){
    file <- paste0(input.path, '/Fold_', k, '_MI_', i, '_', cor.name, '.txt')
    
    mi.cor <- read.table(file, sep = '\t')
    mi.network <- compute.network(mi.cor, cutoff)
    mi.network.list[[i]] <- as.matrix(mi.network)
    
    file <- paste0(input.path, '/Fold_', k, '_MI_', i, '_', cor.name, '_n_b', '.txt')
    write.table(mi.network, file, sep = '\t', quote = F)
  }
  
  # combine network by mode
  mi.network.arr <- array(unlist(mi.network.list), c(n.gene, n.gene, n.pred))
  mi.network.mode <- apply(mi.network.arr, 1:2, compute.mode)
  
  file <- paste0(input.path, '/Fold_', k, '_MI_', cor.name, '_n_b_mode', '.txt')
  write.table(mi.network.mode, file, sep = '\t', quote = F)
}
  