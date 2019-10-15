library(ggfortify)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)


# ------------------------------------------------------------------------------------------------------------------------
# Loading the features matrix and its row names as features names --------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
features_matrix <- read.table(file='result/features_matrix_for_clst.csv',sep=",",header = T)
features_str <- as.character(features_matrix[,1])

if(dim(features_matrix)[2] - 1 < 2)
{
  shinyjs::alert(paste0('There are too few samples in which the selected motifs are recorded.',
                        ' Almost all samples have no mutation in the form of selected motifs.',
                        '\nNo result has been generated!'))
} else {

features_matrix <- features_matrix[,-1]
row.names(features_matrix) <- features_str



# ------------------------------------------------------------------------------------------------------------------------
# Clustering the samples based on the features ---------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
SP <- t(features_matrix)
max_N <- min(10,dim(unique(SP))[1])   # Maximum number of clusters
eval_results <- data.frame(avg_silh_widths=rep(0,max_N))
clst_results <- data.frame(matrix(0,dim(SP)[1],max_N))
row.names(clst_results) <- row.names(SP)
colnames(clst_results) <- paste0('clustering_N_',1:max_N)
centroids_results <- list()
num_of_repeats_of_clutering_for_each_N <- 20   # To make the clustering results reproducible and stable
increment_progress_by <- 0.999/(max_N * num_of_repeats_of_clutering_for_each_N)

start_time <- Sys.time()
for(N in 1:max_N)
{
  
  best_avg_silh_widths <- -10
  best_clst <- rep(0,dim(SP)[1])
  
  for(num_repeat in 1:num_of_repeats_of_clutering_for_each_N)
  {
        # > > > > > > > > > > Clustering the genomes based on selected mutational motifs
        euc_dis <- function(a,b)
        {
          return(sqrt(sum((a-b)^2)))
        } 
        dist_from_cen <- function(i,j) return(euc_dis(SP[i,],centroids[j,]))       # return distance between object  and centroid j
        min_dist  <- function(i) return(which.min(mapply(dist_from_cen,i,c(1:N)))) # return the nearest centroid to the object i
        
        num_iter <- 0
        while(TRUE)
        {
          centroids <- matrix(0,1,dim(SP)[2])        # intialize centroids using K-means++ algorithm
          cen_indxs <- ceiling(runif(1)*(dim(SP)[1]))    # choose the first centroid randomely
          rem_indxs <- setdiff(c(1:(dim(SP)[1])),cen_indxs)
          centroids[1,] <- SP[cen_indxs[1],]
          min_distance_from_centroids <- function(i) # return min distance of object i from currently selected centroids
          {
            dists <- sapply(c(1:(dim(centroids)[1])),function(j) euc_dis(SP[i,],centroids[j,]))
            return(min(dists))
          }
          if(N > 1)
          {
            for(i in 1:(N-1)) # repeat for N-1 times to complete the set of N initial centroids
            {
              D_probs <- sapply(rem_indxs,min_distance_from_centroids)
              if(all(D_probs == 0)){D_probs <- D_probs + 1}
              if(length(D_probs) == 1){
                new_indx <- rem_indxs
              } else {
                new_indx <- sample(rem_indxs,1,prob=D_probs)
              }
              centroids <- rbind(centroids,SP[new_indx,])
              cen_indxs <- c(cen_indxs,new_indx)
              rem_indxs <- setdiff(rem_indxs,new_indx)
            }
          }
          #plot(rbind(SP,centroids),pch=19,col=c(rep('red',dim(SP)[1]),rep('blue',dim(centroids)[1])))
          
          clst0 <- rep(0,dim(SP)[1])     #SP[i,] <-> clst0[i]  (clusters indices: 1,2,3,...,N)
          clst1 <- rep(0,dim(SP)[1])     #SP[i,] <-> clst1[i]  (clusters indices: 1,2,3,...,N)
          
          while(TRUE)
          {
            clst1 <- sapply(c(1:(dim(SP)[1])),min_dist)
            
            for(i in c(1:N))
            {
              members_indices <- which(clst1 == i)
              if(length(members_indices) != 0)
              {
                centroids[i,] <- colMeans(matrix(SP[members_indices,],length(members_indices))) # calculate centroids
              }
            }
            
          
            # cost function: sum-of-squares
            #cost <- sum(sapply(c(1:(dim(SP)[1])),function(i) euc_dis(SP[i,],centroids[clst1[i],])))
            if(identical(clst1,clst0)){break}
            num_iter <- num_iter + 1
            clst0 <- clst1
          }
          
          if(dim(table(clst1)) == N) {break} # if there is any single centroid, run the clustering procedure again!
        }
        #plot(rbind(SP),pch=19,col=clst1)
        
        #costs <- c(costs,cost)
    
        # > > > > > > > > > > Evaluate the clustering based on average silhouette widths
        avg_dist_from_clst <- function(i,j)   # average distance between object i and all members of cluster j
        {
          members <- which(clst1==j)
          return((sum(sapply(members,function(t) euc_dis(SP[i,],SP[t,]))))/length(members))
        }
        silh_width <- function(i)  # silhouette width of cluster i
        {
          a <- avg_dist_from_clst(i,clst1[i])
          out_avgs <- sapply(setdiff(c(1:N),clst1[i]),function(t) avg_dist_from_clst(i,t))
          if(N == 1)
          {
            b <- 1
          } else {
            b <- min(out_avgs)  
          }
          return((b-a)/(max(a,b)))
        }
        avg_silh_widths <- mean(sapply(c(1:(dim(SP)[1])),silh_width))  # average silhouette width of the whole clustering
        
        if(avg_silh_widths > best_avg_silh_widths)
        {
          best_avg_silh_widths <- avg_silh_widths
          best_clst <- clst1
          best_centroids <- centroids
        }

        incProgress(increment_progress_by, detail = paste0("N = ",as.character(N),
                                                           ", at repeat no.",as.character(num_repeat),'/',
                                                           as.character(num_of_repeats_of_clutering_for_each_N)))
        
        cat(paste0("N = ",as.character(N)," : Clustering done with ",
                   as.character(num_iter)," iterations at repeat no.",as.character(num_repeat),'/',
                   as.character(num_of_repeats_of_clutering_for_each_N),
                   " with avg_silh_widths = ",as.character(avg_silh_widths),"\n"))
  }
  
  # save the best result for this N
  eval_results$avg_silh_widths[N] <- best_avg_silh_widths
  clst_results[,N] <- best_clst
  colnames(best_centroids) <- features_str
  centroids_results[[N]] <- best_centroids
  
  cat(paste0("N = ",as.character(N)," : Evaluation done\n\n"))
    
}

N_opt_for_clst <- order(-eval_results$avg_silh_widths)[1]  # obtain the most efficient N for clustering


setProgress(1, detail = paste0("Clustering done!"))


stop_time <- Sys.time()
print(stop_time - start_time)


# save the clustering results
eval_results <- cbind(data.frame(number_of_clusters=c(1:max_N)),eval_results)
write.table(eval_results,file=paste0('output/clustering/Clustering_Evaluations.csv'),sep=",",col.names=TRUE,row.names=FALSE)
clst_results <- cbind(data.frame(Sample_ID=row.names(clst_results)),clst_results)
write.table(clst_results,file=paste0('output/clustering/Clustering_Results.csv'),sep=",",col.names=TRUE,row.names=FALSE)
clst_results <- clst_results[,-1]



# ------------------------------------------------------------------------------------------------------------------------
# Plots of clustering results --------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
withProgress(message = 'Generating the plots', value = 0.5,{
    # PCA plot ------------------------------------------
    S <- data.frame(SP)
    pr_comp <- prcomp(S)
    if(length(which(pr_comp$sdev != 0)) > 2){
      S$'Cluster' <- factor(clst_results[,N_opt_for_clst])
      pca_plt <- autoplot(pr_comp,data=S,colour='Cluster',frame=F,size=4)+
        theme_bw()+
        labs(x='Principal component 1',y='Principal component 2')+
        ggtitle("PCA of input samples in the space of selected features")+
        
        scale_color_hue(name="Clusters:",
                        labels =  paste0('cluster ',c(1:N_opt_for_clst),' (',
                                         as.character(table(clst_results[,N_opt_for_clst])),
                                         ' member',ifelse(table(clst_results[,N_opt_for_clst]) == 1,'','s'),')'))+
        
        theme(axis.title=element_text(size=14,face="bold"),
              legend.title=element_blank(),
              legend.text=element_text(size=14,face="bold"),
              plot.title = element_text(size = 14,face="bold",hjust = 0.5))
      pdf('output/clustering/PCA_of_input_samples.pdf',width = 9,height = 7)
      plot(pca_plt)
      dev.off()
    } else {
      pr_comp <- prcomp(S)
      res <-as.data.frame(as.matrix(S) %*% pr_comp[[2]])
      if(dim(res)[2] < 2){res$PC2 <- 0}
      
      if(max(abs(res$PC1)) != 0){res$PC1 <- res$PC1 / max(abs(res$PC1))}
      if(max(abs(res$PC2)) != 0){res$PC2 <- res$PC2 / max(abs(res$PC2))}
        
      res$Cluster <- factor(clst_results[,N_opt_for_clst])  
      
      pca_plt <- ggplot(res, aes(PC1,PC2),fill = factor(Cluster))+
        geom_point(aes( colour = Cluster),size = 4)+
        theme_bw()+
        labs(x='Principal component 1',y='Principal component 2')+
        ggtitle("PCA of input samples in the space of selected features")+
        theme(axis.title=element_text(size=14,face="bold"),
              legend.title=element_text(size=14,face="bold"),
              legend.text=element_text(size=14,face="bold"),
              plot.title = element_text(size = 14,face="bold",hjust = 0.5))
      pdf('output/clustering/PCA_of_input_samples.pdf',width = 9,height = 7)
      plot(pca_plt)
      dev.off()
  }
  
  
  # Proportions plot ----------------------------------
  if(count_or_proportion == 'proportion' | motif_or_signature == 'signature')
  {
      clusters <- c()
      for(i in 1:N_opt_for_clst){clusters <- c(clusters,rep(i,length(features_str)))}
      cluster_names <- paste('Cluster',clusters)
      dafr <- data.frame(cluster=cluster_names,feature=rep(features_str,N_opt_for_clst),
                         amount=melt(t(centroids_results[[N_opt_for_clst]]))$value)
      
      plt1 <- ggplot(dafr, aes(x = factor(cluster_names), y = amount, fill = feature)) +
        geom_bar(position = position_stack(), stat = "identity", width = .7) +
        geom_text(aes(label = ''), position = position_stack(vjust = 0.5), size = 2)+
        scale_fill_discrete(name = paste0(as.character(k_mer),'-mer ',motif_or_signature,'s:'))+
        theme_bw()+
        theme(axis.title.x=element_blank(),
              axis.text.x=element_text(size=14,face="bold"),
              axis.title.y=element_text(size = 14,face="bold"),
              axis.text.y=element_text(size = 14),
              legend.title=element_text(size = 14),
              legend.text=element_text(size = 14))+
        labs(y = paste0("Average proportions of mutational ",motif_or_signature,"s in each cluster (%)"))
      pdf('output/clustering/Comparing clusters.pdf',width = 9,height = 7)
      plot(plt1)
      dev.off()
  }   
  
    
  ####################################################################################################################   
  ####################################################################################################################   
  ####################################################################################################################   
  # Multiple plot function
  #
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  #
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }    
  ####################################################################################################################   
  ####################################################################################################################   
  ####################################################################################################################   
    

  # Counts plot -----------------------------------
  if(count_or_proportion == 'count')
  {
    get_plts <- function(same_scale)
    {
        plts <- list()
        for(i in 1:dim(centroids_results[[N_opt_for_clst]])[1])
        {
            me <- melt(centroids_results[[N_opt_for_clst]][i,])
            me$features <- factor(row.names(me))
            max_y <- max(me$value)
            if(same_scale){max_y <- max(centroids_results[[N_opt_for_clst]])}
            plts[[i]] <- ggplot(me,aes(x=features,y=value,fill=features))+geom_col()+
              scale_fill_discrete(name = paste0(as.character(k_mer),'-mer ','motif','s:'))+
              theme_bw()+
              theme(axis.title.x=element_blank(),
                    axis.text.x=element_text(size=14,face="bold",angle = 90),
                    axis.title.y=element_text(size = 14,face="bold"),
                    axis.text.y=element_text(size = 14),
                    legend.title=element_text(size = 14),
                    legend.text=element_text(size = 14),
                    title =element_text(size = 15),
                    legend.position="none")+
              labs(y = paste0("Average counts of mutational ",'motifs'," within this cluster"))+
              ggtitle(paste0('Cluster ',i,
                             ifelse(same_scale,'\n(in the same scale compared to other plots)','')))+
              theme(plot.title = element_text(hjust = 0.5))+
              ylim(0,max_y)
        }
        return(plts)
    }
    
    num_plts <- dim(centroids_results[[N_opt_for_clst]])[1]
    pdf(paste0('output/clustering/Comparing clusters (1).pdf'),width = (10+num_plts*5),height = 8)
    multiplot(plotlist =get_plts(same_scale=F),cols = num_plts)
    dev.off()
    pdf(paste0('output/clustering/Comparing clusters (2).pdf'),width = (10+num_plts*5),height = 8)
    multiplot(plotlist =get_plts(same_scale=T),cols = num_plts)
    dev.off()
  }
  setProgress(1, detail = paste0("Finished!"))
})


}