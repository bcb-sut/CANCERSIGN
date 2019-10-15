
library("BSgenome.Hsapiens.UCSC.hg19")
Standard_columns <- c('sample_id','chromosome','position','reference','mutated_to')
setProgress(0.0001, detail = paste0("Opening data table"))

catalog <- read.table('result/input_table.csv',sep = ',',header = T)
if('MT' %in% levels(catalog$chromosome)) {catalog$chromosome[which(catalog$chromosome=='MT')] <- 'M'}

if(method == 1)
{
  setProgress(0.5, detail = paste0("Randomizing..."))
  
  num <- c(0,1,2,3)                 
  names(num) <- c('A','C','G','T')              # make a dictionary
  
  catalog$num_reference <- as.numeric(num[catalog$reference])             # convert reference column to numeric values
  catalog$rands <- floor(runif(dim(catalog)[1])*4)                        # generate random numbers between 0 to 3
  
  #catalog$num_alt2 <- (catalog$num_reference + catalog$rands) %% 4        # generate random numeric values for ALT column 
  catalog$num_alt2 <- (catalog$rands) %% 4        # generate random numeric values for ALT column 
  
  catalog$alt2 <- c('A','C','G','T')[catalog$num_alt2+1]             # reconstruct character values from numeric values
  catalog$mutated_to <- catalog$alt2                                 # replace old ALT column with new random ALT column
  catalog <- catalog[,-c(6:9)]                               # delete temporarily added columns to catalog matrix
  
  write.table(catalog,file=paste0('result/input_table_random_1.csv'),sep=",",col.names=TRUE,row.names=FALSE)
  
  setProgress(1, detail = paste0("Finished!"))
  
} else if (method == 2) {
  
  ids <- unique(catalog$sample_id)

  n_samples <- length(ids)
  
  increment_progress_by <- 0.9999 / n_samples
  
  iii <- 0
  for(id in ids)
  { iii <- iii + 1
    cat(paste('sample',as.character(iii),'started...\n'))
    pat_catalog <- catalog[which(catalog$sample_id == id),]
    chr <- unique(pat_catalog$chromosome)
    for(ch in chr)
    {
      chr_catalog <- pat_catalog[which(pat_catalog$chromosome==ch),]
      N_A <- length(which(chr_catalog$reference=='A'))
      N_C <- length(which(chr_catalog$reference=='C'))
      N_G <- length(which(chr_catalog$reference=='G'))
      N_T <- length(which(chr_catalog$reference=='T'))
      
      L <- as.numeric(seqlengths(Hsapiens)[paste0('chr',ch)])
      
      rand_chr_cat <- data.frame(position=integer(),reference=character())
      
      while(T)
      {
        num_muts <- dim(chr_catalog)[1] - dim(rand_chr_cat)[1]
        if(num_muts == 0){break}
        
        num_A <- N_A - length(which(rand_chr_cat$reference=='A'))
        num_C <- N_C - length(which(rand_chr_cat$reference=='C'))
        num_G <- N_G - length(which(rand_chr_cat$reference=='G'))
        num_T <- N_T - length(which(rand_chr_cat$reference=='T'))
        num_N <- 0
        
        position <- sample(L,3*num_muts) # generate random positionitions...
        rand_cat <- data.frame(position=position,reference=as.character(BSgenome::getSeq(Hsapiens,paste0('chr',ch),position,position)))
        
        A_rows <- which(rand_cat$reference=='A')
        C_rows <- which(rand_cat$reference=='C')
        G_rows <- which(rand_cat$reference=='G')
        T_rows <- which(rand_cat$reference=='T')
        N_rows <- which(rand_cat$reference=='N')
        
        A_diff <- num_A - length(A_rows)
        C_diff <- num_C - length(C_rows)
        G_diff <- num_G - length(G_rows)
        T_diff <- num_T - length(T_rows)
        N_diff <- num_N - length(N_rows)
        
        extra_rows <- c()
        if(A_diff < 0){
          extra_rows <- c(extra_rows,A_rows[1:-A_diff])
          A_diff <- 0}
        if(C_diff < 0){
          extra_rows <- c(extra_rows,C_rows[1:-C_diff])
          C_diff <- 0}
        if(G_diff < 0){
          extra_rows <- c(extra_rows,G_rows[1:-G_diff])
          G_diff <- 0}
        if(T_diff < 0){
          extra_rows <- c(extra_rows,T_rows[1:-T_diff])
          T_diff <- 0}
        if(N_diff < 0){
          extra_rows <- c(extra_rows,N_rows[1:-N_diff])
          N_diff <- 0} 
        rand_cat <- rand_cat[-extra_rows,]
        
        rand_chr_cat <- rbind(rand_chr_cat,rand_cat)
        
        rand_chr_cat <- rand_chr_cat[!duplicated(rand_chr_cat$position),]
      } 
      
      chr_catalog$position[which(chr_catalog$reference=='A')] <- rand_chr_cat$position[which(rand_chr_cat$reference=='A')]
      chr_catalog$position[which(chr_catalog$reference=='C')] <- rand_chr_cat$position[which(rand_chr_cat$reference=='C')]
      chr_catalog$position[which(chr_catalog$reference=='G')] <- rand_chr_cat$position[which(rand_chr_cat$reference=='G')]
      chr_catalog$position[which(chr_catalog$reference=='T')] <- rand_chr_cat$position[which(rand_chr_cat$reference=='T')]
      
      pat_catalog[which(pat_catalog$chromosome==ch),] <- chr_catalog
    }
    catalog[which(catalog$sample_id == id),] <- pat_catalog
    cat(paste('sample',as.character(iii),'finished!\n'))
    
    incProgress(increment_progress_by, detail = paste0('Randomizing... Sample no. ',iii,'/',n_samples,' processed.'))
    
  }
  write.table(catalog,file=paste0('result/input_table_random_2.csv'),sep=",",col.names=TRUE,row.names=FALSE)
  
  setProgress(1, detail = paste0("Finished!"))
}
