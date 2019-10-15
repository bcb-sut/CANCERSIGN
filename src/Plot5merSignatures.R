library(ggplot2)
library(grid)

#rm(list=ls(all=TRUE))
#graphics.off()

col1 <- "#F6766D"         #"deepskyblue2"
col2 <- "#7CAD03"         #"black"
col3 <- "#01BEC4"         #"red"
col4 <- "#C77CFD"         #"darkgray"
col5 <- "#FFC313"         #"limegreen"


plot_signatures_5mer <- function(N)
{
  
  
  motifs <- matrix(read.table("result/selectedM5mer.csv", header = FALSE, sep = ",")[,1])
  
  
  M <- fread(paste0(destination_folder,'5mer_Signatures_(N=',N,').tsv'))
  M <- 100*M
  
  
  M$motifs5mer <- factor(motifs)
  M$mut_types <- factor(sort(rep(c(1:(length(motifs)/16)),16)))
  M$mut_indxs <- factor(c(1:dim(M)[1]))

  
  
  get_color <- function(n)
  {
    nn <- ceiling(n/5)
    all_colors <- rep(c(col1,col2,col3,col4,col5),nn)
    return(all_colors[1:n])
  }
  
  plot_sig_n_5mer <- function(i)
  {
      #top <- 0.3
      top <- max(M[,i,with=F])*1.1
      
      fill_palette <- get_color(length(motifs)/16)
      label_palette <- c()
      for(j in 1:length(fill_palette))
      {
        label_palette <- c(label_palette,rep(fill_palette[j],16))
      }
      
      
          
      g <- ggplot(M, aes_string(x='mut_indxs', y = paste0('Signature_No_',i), fill = 'mut_types', group = 'mut_types'))+
          geom_col(width = 0.7)+
          theme_bw()+
          theme(legend.position="none")+
          ylab("Contribution (%)")+
          labs(title = paste0('(5-mer) Signature No.',as.character(i)))+
          theme(plot.title = element_text(size=20, face="bold",hjust = 0.5),
          panel.border = element_rect(color = "gray"))+
      
          theme(
                axis.text.y=element_text(size = 18,angle = 90, hjust = 0.5,colour = 'black'),
                axis.title.y=element_text(size = 18, margin = margin(t = 0, r = 20, b = 0, l = 0))
          )+
          
          scale_fill_manual(values=fill_palette)+
      
          theme(panel.grid.major = element_blank())+
      
          scale_y_continuous(expand = c(0,0),limits = c(0,top+0.01))+
      
          scale_x_discrete(labels=M$motifs5mer)+
          theme(
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                         colour = label_palette,
                                         size=9,face="bold"),
              axis.ticks.x=element_blank(),
              axis.title.x=element_blank())+
          theme(plot.margin=unit(c(1.5,1,1.5,1),"cm"))
  
      return(g)
  }

  plots <- list()
  for(k in 1:N){plots[[k]] <- plot_sig_n_5mer(k)}
  
  return(plots)
}


