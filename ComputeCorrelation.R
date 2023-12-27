# define input and output path
input.path <- 'Your working directory'
cor.method <- 'spearman'
output.path <- paste0('Your working directory', cor.method)

# define log-normalize function
log.normalize <- function(data, factor = 1e4){
  data <- sweep(data, 1, rowSums(data), '/') * factor
  data <- log(data + 1)
  return (data)
}

n.fold <- 100
n.pred <- 100

for (k in 1:n.fold){
  # load original data
  file <- paste0(input.path, '/Fold_', k, '_Original.txt')
  original <- read.table(file, sep = '\t', header = TRUE)

  # load SI data
  file <- paste0(input.path, '/Fold_', k, '_SI.txt')
  si <- read.table(file, sep = '\t', header = TRUE)
  
  # find indices of predicted genes
  original.cols <- colnames(original)
  predict.cols <- colnames(si)
  predict.idx <- sapply(predict.cols, function(x) which(original.cols == x)[1])
  
  # correlation of original data
  original.new <- cbind(original[,predict.idx], original[,-predict.idx])
  original.cor <- cor(log.normalize(original.new), method = cor.method)
  file <- paste0(output.path, '/Fold_', k, '_Original_', cor.method, '.txt')
  write.table(original.cor, file, sep = '\t', quote = F)
  
  si.new <- cbind(si, original[,-predict.idx])
  si.cor <- cor(log.normalize(si.new), method = cor.method)
  file <- paste0(output.path, '/Fold_', k, '_SI_', cor.method, '.txt')
  write.table(si.cor, file, sep = '\t', quote = F)
  
  # MI data
  for (i in 1:n.pred){
    file <- paste0(input.path, '/Fold_', k, '_MI_', i,'.txt')
    mi <- read.table(file, sep = '\t', header = TRUE)
    mi.new <- cbind(mi, original[,-predict.idx])
    mi.cor <- cor(log.normalize(mi.new), method = cor.method)
    file <- paste0(output.path, '/Fold_', k, '_MI_', i, '_', cor.method, '.txt')
    write.table(mi.cor, file, sep = '\t', quote = F)
  }
}
