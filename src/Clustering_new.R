



library(data.table)
library(ggplot2)
library(NbClust)

fm <- fread('result/features_matrix_for_clst.csv')
features <- fm$features # <<<<<
fm <- fm[,-1]
sample_IDs <- colnames(fm) # <<<<<
fm <- data.table::transpose(fm)



nb <- suppressWarnings(NbClust(fm, distance = "euclidean", min.nc = 2, max.nc = 10, method = "kmeans"))


cluster_assignments <- nb$Best.partition
N_opt_for_clst <- max(cluster_assignments)


centroids <- lapply(1:max(cluster_assignments), function(i) colMeans(fm[which(cluster_assignments == i)]))



# ------------------------------------------------------------------------------------------------------------------------
# Plots of clustering results --------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------

# PCA plot ------------------------------------------
# S <- data.frame(SP)
S <- fm
pr_comp <- prcomp(S)

res <-as.data.frame(as.matrix(S) %*% pr_comp[[2]])
if(dim(res)[2] < 2){res$PC2 <- 0}

# if(max(abs(res$PC1)) != 0){res$PC1 <- res$PC1 / max(abs(res$PC1))}
# if(max(abs(res$PC2)) != 0){res$PC2 <- res$PC2 / max(abs(res$PC2))}

# - - - - - - - - - - - - - - - 
clst_asgns <- factor(cluster_assignments)
tt <- table(clst_asgns)
new_levels <- as.list(names(tt))
names(new_levels) <- paste0('Cluster ',names(tt),paste0(' (',as.numeric(tt),' members)'))
levels(clst_asgns) <- new_levels
# - - - - - - - - - - - - - - - 

res$Cluster <- clst_asgns  

pca_plt <- ggplot(res, aes(PC1,PC2),fill = factor(Cluster))+
    geom_point(aes( colour = Cluster),size = 4)+
    theme_bw()+
    labs(x='Principal component 1',y='Principal component 2')+
    ggtitle("PCA of input samples in the space of selected features")+
    theme(axis.title=element_text(size=14,face="bold"),
          legend.title=element_blank(),
          legend.text=element_text(size=14,face="bold"),
          plot.title = element_text(size = 14,face="bold",hjust = 0.5))

pdf('output/clustering/PCA_of_input_samples.pdf',width = 9,height = 7)
plot(pca_plt)
dev.off()




# Proportions plot ----------------------------------
if(count_or_proportion == 'proportion' | motif_or_signature == 'signature')
{
    
    clusters <- c()
    for(i in 1:N_opt_for_clst){clusters <- c(clusters,rep(i,length(features)))}
    cluster_names <- paste('Cluster',clusters)
    dafr <- data.table(cluster=cluster_names, feature=rep(features,N_opt_for_clst), amount=unlist(centroids))
    
    plt1 <- ggplot(dafr, aes(x = factor(cluster_names), y = amount, fill = feature)) +
        geom_bar(position = position_stack(), stat = "identity", width = .7) +
        geom_text(aes(label = ''), position = position_stack(vjust = 0.5), size = 2)+
        scale_fill_discrete(name = paste0(as.character(k_mer),'-mer ',motif_or_signature,'s:'))+
        theme_bw()+
        ggtitle(paste0("Average proportions of mutational ",motif_or_signature,"s in each cluster (%)"))+
        theme(axis.title.x=element_blank(),
              axis.text.x=element_text(size=14,face="bold"),
              axis.title.y=element_blank(),
              axis.text.y=element_text(size = 14),
              legend.title=element_text(size = 14),
              legend.text=element_text(size = 14),
              plot.title = element_text(size = 14,face="bold",hjust = 0.5))+
    pdf('output/clustering/Comparing clusters.pdf',width = 9,height = 7)
    plot(plt1)
    dev.off()
}   




