
#setProgress(0.5, detail = paste0("Extracting the selected counts of 5-mer motifs from 5-mer catalog"))

library(reshape2)

cat("Transform selected 3-mer motifs into 5-mer motifs...\n")

selected_motifs <- sort(as.numeric(unlist(read.table('result/selected_3_mer_motifs.txt'))))
selected_motifs <- 16*(selected_motifs-1)
motif_indxs <- selected_motifs
for(i in 1:15)
{
  motif_indxs <- rbind(motif_indxs,selected_motifs+i)
}
motif_indxs <- motif_indxs + 1


num2str <- function(num)
{
  nuc <- c('A','C','G','T')
  mut <- c('(C > A)','(C > G)','(C > T)','(T > A)','(T > C)','(T > G)')

  mutation <- mut[((num %/% 16) %/% 16)+1]

  F2 <- nuc[(((num %/% 16) %% 16) %/% 4) +1]
  F3 <- nuc[(((num %/% 16) %% 16) %% 4) +1]

  F1 <- nuc[((num %% 16) %/% 4) +1]
  F4 <- nuc[((num %% 16) %% 4) +1]

  return(paste(F1,F2,mutation,F3,F4))
}

motifs <- t(apply(motif_indxs-1,c(0,1),num2str))

M_5_mer <- read.table(file=paste0('result/M5mer.csv'),sep=",",header = F)

data5 <- as.data.frame(M_5_mer[melt(motif_indxs)$value,1])
if(dim(M_5_mer)[2] > 1)
{
  for(i in 2:dim(M_5_mer)[2])
  {
    data5 <- cbind(data5,as.data.frame(M_5_mer[melt(motif_indxs)$value,i]))
  }
}
data5 <- data.frame(data5)
rownames(data5) <- melt(motifs)$value
colnames(data5) <- paste0('g',1:dim(M_5_mer)[2])

write.table(data.frame(data5),
            file="result/selectedM5mer.csv",
            sep=",", col.names=FALSE, row.names=TRUE)
    


#setProgress(1, detail ="Finished")

