# input path 
cor.method <- 'pearson'
input.path <- paste0('Your working directory', cor.name)

# parameters
n.fold <- 100
n.pred <- 100
n.gene <- 347
n.pred.gene <- 35

for (k in 1:n.fold){
  # load single prediction correlation
  file <- paste0(input.path, '/Fold_', k, '_SI_', cor.method, '.txt')
  si.cor <- read.table(file, sep = '\t')
  
  # load multiple prediction correlations
  mi.cor.list <- vector('list', n.pred)
  for (i in 1:n.pred){
    file <- paste0(input.path, '/Fold_', k, '_MI_', i,'_', cor.method, '.txt')
    mi.cor <- read.table(file, sep = '\t')
    mi.cor.list[[i]] <- as.matrix(mi.cor)
  }
  
  mi.cor.arr <- array(unlist(mi.cor.list), c(n.gene, n.gene, n.pred))
  
  # combine correlations by mean  
  mi.cor.mean <- apply(mi.cor.arr, 1:2, mean)
  rownames(mi.cor.mean) <- rownames(si.cor)
  colnames(mi.cor.mean) <- colnames(si.cor)
  
  output.path <- input.path
  file <- paste0(input.path, '/Fold_', k, '_MI_', cor.name, '_mean.txt')
  write.table(mi.cor.mean, file, quote=F, sep = '\t', row.names = T, col.names = T)  
  
  # combine correlations by median
  mi.cor.median <- apply(mi.cor.arr, 1:2, median)
  rownames(mi.cor.median) <- rownames(si.cor)
  colnames(mi.cor.median) <- colnames(si.cor)
  
  output.path <- input.path
  file <- paste0(input.path, '/Fold_', k, '_MI_', cor.name, '_median.txt')
  write.table(mi.cor.median, file, quote=F, sep = '\t', row.names = T, col.names = T)
  
  # combine correlations by Fisher transformation
  mi.fisher.list <- vector('list', n.pred)
  for (i in 1:n.pred){
    mi.cor <- mi.cor.list[[i]]
    mi.cor.fisher <- 0.5*log((1+mi.cor)/(1-mi.cor))
    mi.fisher.list[[i]] <- mi.cor.fisher
  }
  mi.fisher.arr <- array(unlist(mi.fisher.list), c(n.gene, n.gene, n.pred))
  mi.fisher.mean <- apply(mi.fisher.arr, 1:2, mean)
  mi.cor.fisher <- (exp(2*mi.fisher.mean) - 1) / (exp(2*mi.fisher.mean) + 1)
  diag(mi.cor.fisher) <- 1
  
  rownames(mi.cor.fisher) <- rownames(si.cor)
  colnames(mi.cor.fisher) <- colnames(si.cor)
  
  output.path <- input.path
  file <- paste0(input.path, '/Fold_', k, '_MI_', cor.name, '_fisher.txt')
  write.table(mi.cor.fiser, file, quote=F, sep = '\t', row.names = T, col.names = T)  
}
