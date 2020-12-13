

Clustering <- function(input_file,
                       output_dir,
                       basis,
                       selected_motifs = NULL
                       ){
    
    suppressMessages(suppressWarnings({
        library(data.table)
        library(NbClust)
        library(ggplot2)
    }))
    
    if(basis == 'signatures'){
        
        E <- fread(input_file)
        sampleIDs <- E[['Sample_ID']]
        E <- E[ , - 'Sample_ID']
        E <- as.matrix(E)
        rownames(E) <- sampleIDs
        feature_matrix <- E
        
        if(dim(feature_matrix)[2] == 1){
            stop('cannot perform clustering based on less than two signatures',call. = F)
        }
        
        colnames(feature_matrix) <- gsub('Exposure_to_','',colnames(feature_matrix))
        colnames(feature_matrix) <- gsub('_',' ',colnames(feature_matrix))
        colnames(feature_matrix) <- gsub('No','No.',colnames(feature_matrix))
        
        rs <- rowSums(feature_matrix)
        rs[which(rs == 0)] <- 1
        
        feature_matrix <- feature_matrix / rs
        
    } else if(basis == 'motifs'){
        
        M <- fread(input_file)
        features <- M[['motif']]
        M <- M[ , - 'motif']
        M <- t(M)
        colnames(M) <- features
        feature_matrix <- M
        
        if(length(selected_motifs) < 2){
            stop('The number of selected motifs should not be less than two',call. = F)
        }
        
        selected_motifs <- gsub('\\[','\\\\[',gsub('\\]','\\\\]',selected_motifs))
        selected_cols <- which(sapply(colnames(feature_matrix), function(m) any(sapply(selected_motifs, function(sm) grepl(sm, m)))))
        
        if(length(selected_cols) != length(selected_motifs)){
            stop('Could not match some of the selected motifs to the available motifs',call. = F)
        }
        
        feature_matrix <- feature_matrix[ , selected_cols]
        
    } else {
        stop('Invalid basis type!',call. = F)
    }
    
    
    ########################################
    if(dir.exists(output_dir)){
        unlink(output_dir, recursive = T)
    }
    dir.create(output_dir)
    ########################################
    
    
    #---------------------------------------
    pdf(file = NULL)
    NbClust_result <- tryCatch({
        suppressWarnings(
                NbClust(feature_matrix, distance = "euclidean", min.nc = 2, max.nc = 10, method = "kmeans"))
    }, error = function(err) {
        return('error')
    }
    )
    dev.off()
    #---------------------------------------
    
    if(identical(NbClust_result, 'error')){
        stop('Could not find more than 1 cluster in the input data',call. = F)
    }
    
    cluster_assignments <- NbClust_result$Best.partition

    fwrite(data.table(`Sample ID` = names(cluster_assignments), `Cluster Assignment` = as.numeric(cluster_assignments)), 
           paste0(output_dir,'/cluster_assignments.tsv'), sep = '\t', col.names = T)
    
    
    PCA_plot <- function(feature_matrix, do_scale, plot_title, file_name, output_dir){
        
        pr_comp <- prcomp(feature_matrix, scale. = do_scale)
    
        res <-as.data.frame(feature_matrix %*% pr_comp$rotation)
        if(dim(res)[2] < 2){res$PC2 <- 0}
    
        # - - - - - - - - - - - - - - -
        clst_asgns <- factor(cluster_assignments)
        tt <- table(clst_asgns)
        new_levels <- as.list(names(tt))
        names(new_levels) <- paste0('Cluster ',names(tt),paste0(' (',as.numeric(tt),' members)'))
        levels(clst_asgns) <- new_levels
        # - - - - - - - - - - - - - - -
    
        res$Cluster <- clst_asgns
        
        pca_plt <- ggplot()+
            geom_point(data = res, aes(x=PC1, y=PC2, colour = Cluster),size = 3, stroke = 0, shape=16)+
            theme_bw()+
            labs(x='Principal component 1',y='Principal component 2')+
            ggtitle(plot_title)+
            theme(axis.title=element_text(size=14,face="bold"),
                  legend.title=element_blank(),
                  legend.text=element_text(size=14,face="bold"),
                  plot.title = element_text(size = 14,face="bold",hjust = 0.5),
                  axis.text.x = element_text(size = 9),
                  axis.text.y = element_text(size = 9),
                  panel.grid = element_blank())
    
        pdf(paste0(output_dir,'/',file_name),width = 8.5,height = 6)
        plot(pca_plt)
        dev.off()
    }
    
    PCA_plot(feature_matrix = feature_matrix, 
             do_scale = FALSE,
             plot_title = "PCA of samples",
             file_name = "PCA.pdf",
             output_dir = output_dir)
    
    # PCA_plot(feature_matrix = feature_matrix, 
    #          do_scale = TRUE,
    #          plot_title = "PCA of samples (scaled)",
    #          file_name = "PCA_scaled.pdf",
    #          output_dir = output_dir)
    
    
    
    
    
    
    box_plot <- function(feature_matrix, plot_title, file_name, output_dir){
        # library(ggpubr)
        S <- as.data.table(feature_matrix)
        S <- 100*S
        S$Cluster <- paste0('Cluster ',cluster_assignments)
        
        S_melt <- data.table::melt(S, id.vars = 'Cluster', variable.name = 'Feature')
        
        assignment_list <- list()
        for(n in unique(cluster_assignments)){
            assignment_list[[paste0('Cluster ',n)]] <- paste0('Cluster ',n,' (', length(which(cluster_assignments == n)), ' members)')
        }
        
        box_plt <- ggplot(data = S_melt, aes(x = Cluster, y = value, fill = Cluster))+theme_bw()+
            # geom_point(position=position_jitterdodge(),alpha=0.5, pch = 21)+
            geom_boxplot(outlier.size = 2, outlier.stroke = 0)+
            facet_grid(.~Feature)+
            scale_fill_discrete(labels = assignment_list)+
            # stat_compare_means(method = "anova", size = 3, label.x = 1.4, label.y = 78)+
            ylab('Contribution (%)')+
            xlab('Clusters')+
            ggtitle(plot_title)+
            theme(axis.text.x = element_blank(),
                  #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                  axis.ticks.x = element_blank(),
                  legend.title = element_blank(),
                  legend.text=element_text(size=14,face="bold"),
                  plot.title = element_text(size = 14,face="bold",hjust = 0.5),
                  axis.title=element_text(size=14,face="bold"),
                  panel.grid = element_blank(),
                  strip.text.x = element_text(face="bold"))
        
        pdf(paste0(output_dir,'/',file_name),width = 19,height = 6)
        plot(box_plt)
        dev.off()
    }
    
    if(basis == 'signatures'){
        box_plot(feature_matrix = feature_matrix, 
                 plot_title = "Contribution of mutational signatures to the mutational profile of samples in each cluster",
                 file_name = "box_plots.pdf",
                 output_dir = output_dir)
    } else if(basis == 'motifs'){
        box_plot(feature_matrix = feature_matrix, 
                 plot_title = "Contribution of mutational motifs to the mutational profile of samples in each cluster",
                 file_name = "box_plots.pdf",
                 output_dir = output_dir)
    }
    
    

}

