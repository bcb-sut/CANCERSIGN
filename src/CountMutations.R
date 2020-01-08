

# args = commandArgs(trailingOnly=TRUE)
# 
# input_file_path <- args[1]
# output_dir <- args[2]

CountMutations <- function(input_file_path, output_dir){
    
    start_time <- Sys.time()
    cat('loading libraries')
    s <- Sys.time()
    suppressMessages({
        library(BSgenome.Hsapiens.UCSC.hg19)
        library(data.table)
    })
    f <- Sys.time()
    cat(paste0(' (',format(f-s),')\n'))        
    #------------------------------------------------------------------------------------
    cat('preprocessing')
    s <- Sys.time()
    catalog <- fread(input_file_path)
    catalog$chromosome <- paste0('chr', toupper(gsub('chr', '', as.character(catalog$chromosome))))
    catalog$chromosome <- gsub('chrMT', 'chrM', catalog$chromosome)
    catalog <- catalog[which(catalog$chromosome %in% paste0('chr', c(1:22,'X','Y','M')))]
    
    catalog$`position+2` <- catalog$position + 2
    catalog$position <- catalog$position - 2
    colnames(catalog) <- c('sample_id','chr','pos-2','ref','alt','pos+2')
    catalog$strand <- '+'
    ref_G_or_A_rows <- which(catalog$ref == 'G' | catalog$ref == 'A')
    catalog$strand[ref_G_or_A_rows] <- '-'
    nuc <- c(1:4)
    names(nuc) <- c('A','C','G','T')
    rev_nuc <- c(1:4)
    names(rev_nuc) <- c('T','G','C','A')
    catalog$alt[ref_G_or_A_rows] <- names(rev_nuc[nuc[catalog$alt[ref_G_or_A_rows]]])
    f <- Sys.time()
    cat(paste0(' (',format(f-s),')\n'))        
    #------------------------------------------------------------------------------------
    cat('making Granges object')
    s <- Sys.time()
    gr <- makeGRangesFromDataFrame(df = catalog,
                                   keep.extra.columns = T,
                                   ignore.strand = F,
                                   seqnames.field = 'chr',
                                   start.field = 'pos-2',
                                   end.field = 'pos+2',
                                   strand.field = 'strand')
    f <- Sys.time()
    cat(paste0(' (',format(f-s),')\n'))        
    #------------------------------------------------------------------------------------
    cat('extracting contexts')
    s <- Sys.time()
    L2_L1_ref_R1_R2 <- as.character(BSgenome::getSeq(Hsapiens,gr))
    f <- Sys.time()
    cat(paste0(' (',format(f-s),')\n'))        
    #------------------------------------------------------------------------------------
    cat('matching motifs')
    s <- Sys.time()
    L2_L1_ref_R1_R2_alt <- paste0(L2_L1_ref_R1_R2,mcols(gr)$alt)
    motifs <- c()            
    for(mut in c('CA','CG','CT','TA','TC','TG')) {
        for(left1 in c('A','C','G','T')) {
            for(right1 in c('A','C','G','T')) {
                for(left2 in c('A','C','G','T')) {
                    for(right2 in c('A','C','G','T')) {
                        motifs <- c(motifs,paste0(left2,left1,substr(mut,1,1),right1,right2,substr(mut,2,2))) #L2_L1_ref_R1_R2_alt
                    }
                }
            }
        }
    }
    motifs_nums <- match(L2_L1_ref_R1_R2_alt,motifs)
    f <- Sys.time()
    cat(paste0(' (',format(f-s),')\n'))        
    #------------------------------------------------------------------------------------
    cat('making motif matrixes')
    s <- Sys.time()
    sample_ids <- mcols(gr)$sample_id
    sample_names <- unique(sample_ids)
    sample_motif_list <- lapply(sample_names,function(s){
        cata <- rep(0,16*16*6)
        tt <- table(motifs_nums[which(sample_ids == s)])
        cata[as.numeric(names(tt))] <- as.numeric(tt)
        return(cata)
    })
    f <- Sys.time()
    cat(paste0(' (',format(f-s),')\n'))        
    #------------------------------------------------------------------------------------
    cat('saving the results')
    s <- Sys.time()
    M5mer <- simplify2array(sample_motif_list)
    M3mer <- t(simplify2array(lapply(c(0:95),function(i) colSums(M5mer[c((i*16+1):(i*16+16)),]))))
    
    
    
    motifs_3mer <- c()
    for(mut in c('C>A','C>G','C>T','T>A','T>C','T>G')) {
        for(left1 in c('A','C','G','T')) {
            for(right1 in c('A','C','G','T')) {
                motifs_3mer <- c(motifs_3mer,paste0(left1,'[',mut,']',right1))
            }
        }
    }
    motifs_5mer <- c()
    for(mut in c('C>A','C>G','C>T','T>A','T>C','T>G')) {
        for(left1 in c('A','C','G','T')) {
            for(right1 in c('A','C','G','T')) {
                for(left2 in c('A','C','G','T')) {
                    for(right2 in c('A','C','G','T')) {
                        motifs_5mer <- c(motifs_5mer,paste0(left2,left1,'[',mut,']',right1,right2))
                    }
                }
            }
        }
    }
    
    if(!dir.exists(output_dir)){
        dir.create(output_dir, recursive = T)
    }
    
    M3mer <- as.data.table(M3mer)
    colnames(M3mer) <- sample_names
    M3mer_complete <- cbind(data.table(motif = motifs_3mer), M3mer)
    M3mer_path <- paste0(output_dir,'/','motif-matrix-3mer.csv')
    fwrite(M3mer_complete, file=M3mer_path,sep=",",col.names=T,row.names=F)
    
    M5mer <- as.data.table(M5mer)
    colnames(M5mer) <- sample_names
    M5mer_complete <- cbind(data.table(motif = motifs_5mer), M5mer)
    M5mer_path <- paste0(output_dir,'/','motif-matrix-5mer.csv')
    fwrite(M5mer_complete, file=M5mer_path,sep=",",col.names=T,row.names=F)
    
    f <- Sys.time()
    cat(paste0(' (',format(f-s),')\n'))        
    #------------------------------------------------------------------------------------
    finish_time <- Sys.time()
    cat(paste0('total elapsed time: (',format(finish_time-start_time),')'))
    
}    
    
