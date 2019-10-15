

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
cat('importing input_table.csv')
s <- Sys.time()
catalog <- fread(paste0('result/',input_table_file_name,'.csv'),sep = ',',header = T)

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

if(save_M5mer == T){fwrite(data.frame(M5mer),file=paste0('result/',M5mer_file_name,'.csv'),sep=",",col.names=FALSE,row.names=FALSE)}
if(save_M3mer == T){fwrite(data.frame(M3mer),file=paste0('result/',M3mer_file_name,'.csv'),sep=",",col.names=FALSE,row.names=FALSE)}

write.table(sample_names,'result/Sample_IDs.txt',row.names = F,col.names = F)

f <- Sys.time()
cat(paste0(' (',format(f-s),')\n'))        
#------------------------------------------------------------------------------------
finish_time <- Sys.time()
cat(paste0('total elapsed time: (',format(finish_time-start_time),')'))


