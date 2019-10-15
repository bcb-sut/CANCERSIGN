
motif_num_to_motif_str <- function(num,k_mer)
{
  num <- num - 1
  nuc <- c('A','C','G','T')
  mut <- c('(C > A)',
           '(C > G)',
           '(C > T)',
           '(T > A)',
           '(T > C)',
           '(T > G)')
  if(k_mer == 5)
  {
    mutation <- mut[((num %/% 16) %/% 16)+1]
    F2 <- nuc[(((num %/% 16) %% 16) %/% 4) +1]
    F3 <- nuc[(((num %/% 16) %% 16) %% 4) +1]
    F1 <- nuc[((num %% 16) %/% 4) +1]
    F4 <- nuc[((num %% 16) %% 4) +1]
    return(paste(F1,F2,mutation,F3,F4))
  }
  if(k_mer == 3)
  {
    mutation <- mut[(num %/% 16)+1]
    F1 <- nuc[((num %% 16) %/% 4) +1]
    F2 <- nuc[((num %% 16) %% 4) +1]
    return(paste(F1,mutation,F2))
  }
}


Samples_IDs <- as.character(unlist(read.table('result/Samples_IDs.txt')))


if(motif_or_signature=='motif') {
  features_indx <- sort(as.numeric(unlist(read.table(paste0('result/selected_',as.character(k_mer),'mer_motifs_for_clst.txt')))))
  features_str <- motif_num_to_motif_str(features_indx,k_mer)
  M3mer <- read.table(file=paste0('result/M',as.character(k_mer),'mer.csv'),sep=",")
  features_matrix <- M3mer[features_indx,]
  row.names(features_matrix) <- features_str
  colnames(features_matrix) <- Samples_IDs
  
  if(count_or_proportion=='proportion') {
    sums <- colSums(features_matrix)
    zero_indxs <- which(sums==0)       # Throw out the samples which have zero values for all features
    if(length(zero_indxs)!=0) {
      features_matrix <- features_matrix[,-zero_indxs]
      sums <- sums[-zero_indxs]
    }
    features_matrix <- (t(t(features_matrix) / sums))*100
    row.names(features_matrix) <- features_str
  }
  
} else if(motif_or_signature=='signature') {
  N_opt <- read.table(paste0('output/signatures/',as.character(k_mer),'_mer/N_opt.txt'),sep = ',')[1,1]
  features_matrix <- read.table(paste0('output/signatures/',as.character(k_mer),'_mer/E-n-',as.character(N_opt),'.txt'),sep = ' ')
  row.names(features_matrix) <- paste('Signature',c(1:dim(features_matrix)[1]))
  colnames(features_matrix) <- Samples_IDs
  
  sums <- colSums(features_matrix)
  zero_indxs <- which(sums==0)       # Throw out the samples which have zero values for all features
  if(length(zero_indxs)!=0) {
    features_matrix <- features_matrix[,-zero_indxs]
    sums <- sums[-zero_indxs]
  }
  features_matrix <- (t(t(features_matrix) / sums))*100
  
}

features_matrix <- cbind(data.frame(features = row.names(features_matrix)),features_matrix)

write.table(features_matrix,file="result/features_matrix_for_clst.csv",sep=",", col.names=TRUE, row.names=FALSE)


setProgress(1, detail = paste0("Preparation done!"))













