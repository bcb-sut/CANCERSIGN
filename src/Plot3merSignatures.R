library(ggplot2)
library(grid)
library(ggpubr)
library(gridExtra)

#rm(list=ls(all=TRUE))
#graphics.off()

plot_signatures_3mer <- function(N)
{
  
  M <- fread(paste0(destination_folder,'3mer_Signatures_(N=',N,').tsv'))
  M <- 100*M
  
  A1 <- c()
  for(i in c(0:5)) A1 <- c(A1,c(1:4)+4*0+16*i)  
  C1 <- c()
  for(i in c(0:5)) C1 <- c(C1,c(1:4)+4*1+16*i)  
  G1 <- c()
  for(i in c(0:5)) G1 <- c(G1,c(1:4)+4*2+16*i)  
  T1 <- c()
  for(i in c(0:5)) T1 <- c(T1,c(1:4)+4*3+16*i)  
  
  C21 <- c(1:16)+16*0
  C22 <- c(1:16)+16*1
  C23 <- c(1:16)+16*2
  T21 <- c(1:16)+16*3
  T22 <- c(1:16)+16*4
  T23 <- c(1:16)+16*5
  
  A3 <- c(0:23)*4 + 1
  C3 <- c(0:23)*4 + 2
  G3 <- c(0:23)*4 + 3
  T3 <- c(0:23)*4 + 4
  
  M$mut_types <- factor(c(rep('C > A',16),rep('C > G',16),rep('C > T',16),
                          rep('T > A',16),rep('T > C',16),rep('T > G',16)))
  M$mut_indxs <- factor(c(1:96))
  
  
  col1 <- "#F6766D"         #"deepskyblue2"
  col2 <- "#7CAD03"         #"black"
  col3 <- "#01BEC4"         #"red"
  col4 <- "#C77CFD"         #"darkgray"
  col5 <- "#FFC313"         #"limegreen"
  col6 <- "palevioletred1"  #"palevioletred1"
  
  
  plot_sig_n_3mer <- function(n)
  {
    
    #top <- 0.3
    top <- max(M[,n,with=F])
    
    gg <- ggplot(M, aes_string(x='mut_indxs', y = paste0('Signature_No_',n), fill = 'mut_types', group = 'mut_types'))+
          geom_col(width = 0.7)+
          geom_text(aes(label=mut_id), vjust=-1.5, colour="black",size = 2)+
          theme_bw()+
          theme(legend.position="none")+
          ylab("Contribution (%)")+
          labs(title = paste0('Signature No. ',as.character(n)))+
          theme(plot.title = element_text(size=28, 
                                          #face="bold",
                                          hjust = 0.5,
                                          vjust = 5),
                panel.border = element_rect(color = "gray"))+

          
          scale_fill_manual(values=c(col1,col2,col3,col4,col5,col6))+
          
          
          theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                
                axis.text.y=element_text(size = 25,angle = 90, hjust = 0.5,colour = 'black'),
                axis.title.y=element_text(size = 25, margin = margin(t = 0, r = 20, b = 0, l = 0))
                )+
      
          theme(panel.grid.minor.y = element_blank(),
                panel.grid.major.x = element_blank())+
        
          theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
          
    
          scale_y_continuous(expand = c(0.025,0),limits = c(-3*top/30-top/30,top*1.3))+
    
          annotate(geom="text", x=(A1+0.02), y=-3*top/25, label='A',
                   color="black",angle = 90,size = 3,family='sans')+
          annotate(geom="text", x=(C1+0.02), y=-3*top/25, label='C',
                   color="black",angle = 90,size = 3,family='sans')+
          annotate(geom="text", x=(G1+0.02), y=-3*top/25, label='G',
                   color="black",angle = 90,size = 3,family='sans')+
          annotate(geom="text", x=(T1+0.02), y=-3*top/25, label='T',
                   color="black",angle = 90,size = 3,family='sans')+
    
    
          annotate(geom="text", x=(C21+0.02), y=-3*top/30+top/45, label='C',
                   color=col1,angle = 90,size = 3,fontface = 'bold',family='sans')+
          annotate(geom="text", x=(C22+0.02), y=-3*top/30+top/45, label='C',
                   color=col2,angle = 90,size = 3,fontface = 'bold',family='sans')+
          annotate(geom="text", x=(C23+0.02), y=-3*top/30+top/45, label='C',
                   color=col3,angle = 90,size = 3,fontface = 'bold',family='sans')+
    
          annotate(geom="text", x=(T21+0.02), y=-3*top/30+top/45, label='T',
                   color=col4,angle = 90,size = 3,fontface = 'bold',family='sans')+
          annotate(geom="text", x=(T22+0.02), y=-3*top/30+top/45, label='T',
                   color=col5,angle = 90,size = 3,fontface = 'bold',family='sans')+
          annotate(geom="text", x=(T23+0.02), y=-3*top/30+top/45, label='T',
                   color=col6,angle = 90,size = 3,fontface = 'bold',family='sans')+
    
    
          annotate(geom="text", x=(A3+0.02), y=-3*top/30+2*top/31.5, label='A',
                   color="black",angle = 90,size = 3,family='sans')+
          annotate(geom="text", x=(C3+0.02), y=-3*top/30+2*top/31.5, label='C',
                   color="black",angle = 90,size = 3,family='sans')+
          annotate(geom="text", x=(G3+0.02), y=-3*top/30+2*top/31.5, label='G',
                   color="black",angle = 90,size = 3,family='sans')+
          annotate(geom="text", x=(T3+0.02), y=-3*top/30+2*top/31.5, label='T',
                   color="black",angle = 90,size = 3,family='sans')+
      
          geom_rect(xmin= 1-0.33, xmax=16+0.33, ymin=top*1.1, ymax=top*1.19,fill = col1)+
          geom_rect(xmin=17-0.33, xmax=32+0.33, ymin=top*1.1, ymax=top*1.19,fill = col2)+
          geom_rect(xmin=33-0.33, xmax=48+0.33, ymin=top*1.1, ymax=top*1.19,fill = col3)+
          geom_rect(xmin=49-0.33, xmax=64+0.33, ymin=top*1.1, ymax=top*1.19,fill = col4)+
          geom_rect(xmin=65-0.33, xmax=80+0.33, ymin=top*1.1, ymax=top*1.19,fill = col5)+
          geom_rect(xmin=81-0.33, xmax=96+0.33, ymin=top*1.1, ymax=top*1.19,fill = col6)+
      
      
          annotate(geom="text", x=8, y=top*1.26, label='C > A',
                   color="black",size = 9,
                   #fontface = 'bold',
                   family='sans')+
      
          annotate(geom="text", x=8+16, y=top*1.26, label='C > G',
                   color="black",size = 9,
                   #fontface = 'bold',
                   family='sans')+
      
          annotate(geom="text", x=8+16*2, y=top*1.26, label='C > T',
                   color="black",size = 9,
                   #fontface = 'bold',
                   family='sans')+
      
          annotate(geom="text", x=8+16*3, y=top*1.26, label='T > A',
                   color="black",size = 9,
                   #fontface = 'bold',
                   family='sans')+
      
          annotate(geom="text", x=8+16*4, y=top*1.26, label='T > C',
                   color="black",size = 9,
                   #fontface = 'bold',
                   family='sans')+
      
          annotate(geom="text", x=8+16*5, y=top*1.26, label='T > G',
                   color="black",size = 9,
                   #fontface = 'bold',
                   family='sans')
    
    return(gg)    
  }
  
  M$mut_id <- ''
  
  plots <- list()
  for(i in 1:N){plots[[i]] <- plot_sig_n_3mer(i)}
  
  return(plots)
}


