suppressMessages(library(shiny))
suppressMessages(library(shinyjs))

suppressMessages(library(shinyWidgets))

suppressMessages(library(rhandsontable))
suppressMessages(library(parallel))
suppressMessages(library(data.table))

suppressMessages(library(ggplot2))


rm(list=ls(all=TRUE))

# Set the path for working directory
setwd('..')


cat('\n--------------------------------------------------------------------------------------------\n')
cat('--------------------------------------------------------------------------------------------\n')
cat('---------------------------------------- CANCERSIGN ----------------------------------------\n')
cat('--------------------------------------------------------------------------------------------\n')
cat('--------------------------------------------------------------------------------------------\n')

cat(paste0('\nCurrent working directory: ',getwd(),'\n'))


write.table(Sys.getpid(),file = 'temp/pid.txt',row.names = FALSE,col.names = FALSE,append = TRUE)

############################################################################################################################################
###########################################################                      ###########################################################
###########################################################   Global variables   ###########################################################
###########################################################                      ###########################################################
############################################################################################################################################

# ICGC_required_columns <- c('icgc_sample_id','chromosome','chromosome_start',
#                            'chromosome_end','reference_genome_allele','mutated_to_allele')

Standard_columns <- c('sample_id','chromosome','position','reference','mutated_to')


# a <- read.table('data/Breast_whole.csv',sep = ',',header = F)
# a <- a[,-c(1,5,8)]
# colnames(a) <- c('sample_id','chromosome','position','reference','mutated_to')
# write.table(a,'data/Breast_whole.csv',sep = ',',row.names = F,col.names = T)




num2str <- function(num)
{
  nuc <- c('A','C','G','T')
  mut <- c('C','C','C','T','T','T')
  mutation <- mut[(num %/% 16)+1]
  Fl <- nuc[((num %% 16) %/% 4) +1]
  Fr <- nuc[((num %% 16) %% 4) +1]
  return(paste(Fl,mutation,Fr))
}
c_a <- (16*0+1):(16*1)
names(c_a) <- num2str(c_a-1)
c_g <- (16*1+1):(16*2)
names(c_g) <- num2str(c_g-1)
c_t <- (16*2+1):(16*3)
names(c_t) <- num2str(c_t-1)
t_a <- (16*3+1):(16*4)
names(t_a) <- num2str(t_a-1)
t_c <- (16*4+1):(16*5)
names(t_c) <- num2str(t_c-1)
t_g <- (16*5+1):(16*6)
names(t_g) <- num2str(t_g-1)


initial_table_5mer_motifs_for_clst <- data.frame(Selected = FALSE,
                                                 Flanking1 = factor(rep('',10),levels = c('A','C','G','T')),
                                                 Flanking2 = factor(rep('',10),levels = c('A','C','G','T')),
                                                 Mutation = factor(rep('',10),levels = c('C > A',
                                                                                         'C > G',
                                                                                         'C > T',
                                                                                         'T > A',
                                                                                         'T > C',
                                                                                         'T > G')),
                                                 Flanking3 = factor(rep('',10),levels = c('A','C','G','T')),
                                                 Flanking4 = factor(rep('',10),levels = c('A','C','G','T'))) 


initial_table_3mer_motifs_for_clst <- data.frame(Selected = FALSE,
                                                 Flanking1 = factor(rep('',10),levels = c('A','C','G','T')),
                                                 Mutation = factor(rep('',10),levels = c('C > A',
                                                                                         'C > G',
                                                                                         'C > T',
                                                                                         'T > A',
                                                                                         'T > C',
                                                                                         'T > G')),
                                                 Flanking2 = factor(rep('',10),levels = c('A','C','G','T'))) 







#############################################################################################################################################
###########################################################                       ###########################################################
###########################################################   UI configurations   ###########################################################
###########################################################                       ###########################################################
#############################################################################################################################################

ui <- fluidPage(

  useShinyjs(),
  
  div(id = "preprocessing_ui",
      navbarPage(inverse=TRUE,"CANCERSIGN",
                 tabPanel("Preprocessing input file",

                          wellPanel(id="panelA",
                            tags$div(class="header", checked=NA,tags$h4(
                              strong('Put your data file inside "data" folder. Then type the name of input file (with extension) below:'),
                              style="color:#0060DB")),
                            
                            
                            
                            textInput("input_data_file_name", label=NA,width = '50%'),
                            hr(),
                            
                            tags$div(class="header", checked=NA,tags$h4(
                              strong("How many CPU cores do you want to use?"),style="color:#0060DB")),
                            
                            sliderInput("CPU_cores", "Number of CPU cores for parallelization",
                                         min =1 , max = detectCores(),
                                         width = '50%',
                                         value = ceiling(detectCores()/2), step = 1),
                            
                            hr(),
                            
                            actionButton("see_input_file_name", "Import new data and perform preprocessing",class="btn-success"),
                            
                            textOutput("success_for_input_file_name"),
                            tags$head(tags$style("#success_for_input_file_name{color: #1BDB00;font-size: 100%;}")),
                            
                            hr(),
                            
                            actionButton("skip_preprocessing", "Don't import new data\nand skip preprocessing",class="btn-warning")
                            ),
                         
                          
                          shinyjs::hidden(
                            wellPanel(id="panelB",
                                      tags$div(class="header", checked=NA,
                                               tags$h4(strong(paste0('The tool is going to delete all contents of "output" and "result" folders.', 
                                                                     ' You can make a backup of them before continuing.')),
                                                       style="color:#ff8421")),
                                      hr(),
                                      actionButton("continue_with_preprocessing", "Continue",class="btn-success")
                                      )
                            ),
                          
                          
                          shinyjs::hidden(
                            wellPanel(id="panelC",
                        
                                      tags$div(class="header", checked=NA,
                                               tags$h4(strong('Select one of the options for "output/signatures" folder:'),
                                                       style="color:#0060DB")),
                                      
                                      radioButtons(inputId="keep_or_not_signatures",
                                                   label=NA,
                                                   choices = list("Clear the previously obtained results for signatures " = 1,
                                                                  "Keep the previously obtained results for 3-mer signatures only" = 2,
                                                                  "Keep the previously obtained results for 5-mer signatures only" = 3,
                                                                  "Keep the previously obtained results for both types of signatures" = 4),
                                                   selected = 4),
                                      
                                      hr(),
                        
                                      tags$div(class="header", checked=NA,
                                               tags$h4(strong('Select one of the options for "output/clustering" folder:'),
                                                       style="color:#0060DB")),
                        
                                      radioButtons(inputId="keep_or_not_clustering",
                                                   label=NA,
                                                   choices = list("Clear the previously obtained results for clustering " = 1, 
                                                                  "Keep the previously obtained results for clustering" = 2), 
                                                   selected = 2),
                                      
                                      hr(),hr(),
                                      
                                      actionButton("continue", "Continue",class="btn-success")
                                      
                                      )
                            )
                          ),
                 
                 tags$style(type = 'text/css',
                            '.navbar {background-color: #383838;}',
                            '.navbar-default .navbar-brand {background-color: #383838;color: #FF8300;}')
                 )
      ),
  
  #-----------------------------------------------------------------------------------------------------------------------------------------
  
  hidden(div(id = "processing_ui",
  navbarPage(inverse=TRUE,"CANCERSIGN",
  tabPanel("3-mer signatures",
           
           titlePanel("Extract mutational signatures based on 3-mer motifs"),
           
           
           tabsetPanel(
               tabPanel(h4("New Analysis"), 
                        wellPanel(id="panel_3mer_options",
                                  tags$div(class="header", checked=NA,tags$h4(strong("Select the accuracy level:"),
                                                                              style="color:#0060DB")),
                                  radioButtons(inputId="accuracy_3mer",
                                               label=NA,
                                               choices=c('Moderate','High','Very High')),
                                  
                                  tags$div(class="header", checked=NA,tags$h4(strong("Select N-Max (the tool will search for optimum number of signatures from 1 to N-Max):"),
                                                                              style="color:#0060DB")),
                                  sliderInput("NMF_Max_N_3mer", "",
                                              width = '50%',
                                              min = 2, max = 50,
                                              value = 10, step = 1),
                                  
                                  fluidRow(column(6,actionButton("start_3mer", "Start",class="btn-success", width = '150px')))
                        )
                        ),
               tabPanel(h4("Results"),
                        wellPanel(
                            fluidRow(column(11,
                                     tags$div(class="header", checked=NA,tags$h4(
                                         strong('Get the evaluation diagram for determining the optimum number of clusters:'),
                                         style="color:#0060DB"))
                            ),
                                     column(1,
                                            downloadButton("download_3mer_eval", "",class="btn-info", width = '100%'),
                                            tags$style(type='text/css', "#download_3mer_eval {height: 35px; margin-bottom:12px}"))),
                            hr(),
                            fluidRow(column(9,
                                            tags$div(class="header", checked=NA,tags$h4(
                                                strong('Get the plots of mutational signatures (select the optimum number of signatures):'),
                                                style="color:#0060DB"))
                                            ),
                                     column(2,align="right",
                                            numericInput('num_3mer_sigs_to_download', label=NULL, value=0, min = 2, max = 50, step = 1,width = '50%'),
                                            tags$style(type='text/css', "#num_3mer_sigs_to_download {height: 35px; margin-bottom:2px}")),
                                     column(1,
                                            downloadButton("download_3mer_signatures", "",class="btn-info",width = '100%'),
                                            tags$style(type='text/css', "#download_3mer_signatures {height: 35px; margin-bottom:2px}"))),
                            hr(),
                            fluidRow(column(11,
                                            tags$div(class="header", checked=NA,tags$h4(
                                                strong('Get output files:'),
                                                style="color:#0060DB"))
                                            ),
                                     column(1,
                                            downloadButton("download_files_3mer", "",class="btn-info", width = '100%'),
                                            tags$style(type='text/css', "#download_files_3mer {height: 35px; margin-bottom:12px}")))
                            )
                        )
           )
  ),

  #-----------------------------------------------------------------------------------------------------------------------------------------
  
  tabPanel("5-mer signatures",
           
           titlePanel("Extract mutational signatures based on 5-mer motifs"),
           
           tabsetPanel(
               tabPanel(h4("New Analysis"), 
                        
                   wellPanel(id="panel_5mer_options",
                     
                             tags$div(class="header", checked=NA,tags$h4(strong("Select the 3-mer motifs which you want to expand:"),
                                                                         style="color:#0060DB")),
                             
                             fluidRow(column(2,checkboxGroupInput('C_A_for_sig',' N (C > A) N',choices = c_a)),
                                      column(2,checkboxGroupInput('C_G_for_sig',' N (C > G) N',choices = c_g)),
                                      column(2,checkboxGroupInput('C_T_for_sig',' N (C > T) N',choices = c_t)),
                                      column(2,checkboxGroupInput('T_A_for_sig',' N (T > A) N',choices = t_a)),
                                      column(2,checkboxGroupInput('T_C_for_sig',' N (T > C) N',choices = t_c)),
                                      column(2,checkboxGroupInput('T_G_for_sig',' N (T > G) N',choices = t_g))),
                   
                             textOutput("select_motif_error_5mer"),
                             tags$head(tags$style("#select_motif_error_5mer{color: red;font-size: 100%;}")),
                             
                             hr(),
        
                             
                             tags$div(class="header", checked=NA,tags$h4(strong("Select the accuracy level:"),
                                                                         style="color:#0060DB")),
                             radioButtons(inputId="accuracy_5mer",
                                          label=NA,
                                          choices=c('Moderate','High','Very High')), 
        
                             
                             
                             tags$div(class="header", checked=NA,tags$h4(strong("Select N-Max (the tool will search for optimum number of signatures from 1 to N-Max):"),
                                                                         style="color:#0060DB")),
                             sliderInput("NMF_Max_N_5mer", "",
                                         width = '50%',
                                         min = 2, max = 50,
                                         value = 10, step = 1),
                             
                             
                             fluidRow(column(6,actionButton("start_5mer", "Start",class="btn-success", width = '150px')))
                             
                   )
               ),
               tabPanel(h4("Results"),
                        wellPanel(
                            fluidRow(column(11,
                                            tags$div(class="header", checked=NA,tags$h4(
                                                strong('Get the evaluation diagram for determining the optimum number of clusters:'),
                                                style="color:#0060DB"))
                            ),
                            column(1,
                                   downloadButton("download_5mer_eval", "",class="btn-info", width = '100%'),
                                   tags$style(type='text/css', "#download_5mer_eval {height: 35px; margin-bottom:12px}"))),
                            hr(),
                            fluidRow(column(9,
                                            tags$div(class="header", checked=NA,tags$h4(
                                                strong('Get the plots of mutational signatures (select the optimum number of signatures):'),
                                                style="color:#0060DB"))
                            ),
                            column(2,align="right",
                                   numericInput('num_5mer_sigs_to_download', label=NULL, value=0, min = 2, max = 50, step = 1,width = '50%'),
                                   tags$style(type='text/css', "#num_5mer_sigs_to_download {height: 35px; margin-bottom:2px}")),
                            column(1,
                                   downloadButton("download_5mer_signatures", "",class="btn-info",width = '100%'),
                                   tags$style(type='text/css', "#download_5mer_signatures {height: 35px; margin-bottom:2px}"))),
                            hr(),
                            fluidRow(column(11,
                                            tags$div(class="header", checked=NA,tags$h4(
                                                strong('Get output files:'),
                                                style="color:#0060DB"))
                            ),
                            column(1,
                                   downloadButton("download_files_5mer", "",class="btn-info", width = '100%'),
                                   tags$style(type='text/css', "#download_files_5mer {height: 35px; margin-bottom:12px}")))
                        )
               )
           )
  ),
  
  #-----------------------------------------------------------------------------------------------------------------------------------------
  
  tabPanel("Clustering",
           
           titlePanel("Clustering the samples based on mutational motifs or proportion of signatures"),
           wellPanel(id="panel_clustering_options_1",
                     
                     tags$div(class="header", checked=NA,tags$h4(strong("Based on which method do you want to perform the clustering?"),
                                                                 style="color:#0060DB")),
                     radioButtons(inputId="clustering_options_1",
                                  label=NA,
                                  choices=c('Cluster based on 3-mer mutational motifs',
                                            'Cluster based on 5-mer mutational motifs',
                                            'Cluster based on proportions of 3-mer signatures',
                                            'Cluster based on proportions of 5-mer signatures'),
                                  selected = character(0)),
                     
                     hr(),
                     
                     conditionalPanel("input.clustering_options_1 == 'Cluster based on 3-mer mutational motifs'",
                                      
                                      tags$div(class="header", checked=NA,
                                               tags$h4(strong("Select the 3-mer motifs for which you want to perform clustering:"),
                                                                                    style="color:#0060DB")),
                                      rHandsontableOutput("select_3mer_for_clst"),
                                        
                                      hr(),
                                        
                                      tags$div(class="header", checked=NA,
                                                tags$h4(strong("Cluster based on count values or proportions of mutational motifs?"),
                                                        style="color:#0060DB")),
                                        
                                      radioButtons(inputId="count_or_proportion_3mer",
                                                   label=NA,
                                                   choices=c('Cluster based on "count" values of 3-mer mutational motifs',
                                                            'Cluster based on "proportions" of 3-mer mutational motifs')),
                                      hr()
                                    ),
                     
                     conditionalPanel("input.clustering_options_1 == 'Cluster based on 5-mer mutational motifs'",
                                        
                                      tags$div(class="header", checked=NA,
                                               tags$h4(strong("Select the 5-mer motifs for which you want to perform clustering:"),
                                                       style="color:#0060DB")),
                                        
                                      rHandsontableOutput("select_5mer_for_clst"),
                                        
                                      hr(),
                                        
                                      tags$div(class="header", checked=NA,
                                               tags$h4(strong("Cluster based on count values or proportions of mutational motifs?"),
                                                       style="color:#0060DB")),
                                        
                                      radioButtons(inputId="count_or_proportion_5mer",
                                                   label=NA,
                                                   choices=c('Cluster based on "count" values of 5-mer mutational motifs',
                                                             'Cluster based on "proportions" of 5-mer mutational motifs')),
                                      hr()
                                    ),
                     
                     conditionalPanel('input.clustering_options_1',
                                      actionButton("start_clustering", "Start clustering",class="btn-success", width = '150px')
                                    ),
                     
                     textOutput("select_motif_error_3mer_clst"),
                     tags$head(tags$style("#select_motif_error_3mer_clst{color: red;font-size: 100%;}")),
                     
                     textOutput("select_motif_error_5mer_clst"),
                     tags$head(tags$style("#select_motif_error_5mer_clst{color: red;font-size: 100%;}")),
                     
                     textOutput("error_sig_proportion_3mer_clst"),
                     tags$head(tags$style("#error_sig_proportion_3mer_clst{color: red;font-size: 100%;}")),
                     
                     textOutput("error_sig_proportion_5mer_clst"),
                     tags$head(tags$style("#error_sig_proportion_5mer_clst{color: red;font-size: 100%;}"))
 
          ),
          
          shinyjs::hidden(
            wellPanel(id='final_message_clustering',
                      tags$div(class="header", checked=NA,
                               tags$h4(strong('You can find the results in this directory: output/clustering/'),
                                       style="color:#0060DB"))
            )
          )
           
  )
  )
  )
  )
)






############################################################################################################################################
###########################################################                      ###########################################################
###########################################################   Server functions   ###########################################################
###########################################################                      ###########################################################
############################################################################################################################################

server <- function(input, output, session) {
  
    observeEvent(input$see_input_file_name,{
        if(isolate(input$input_data_file_name) == ''){
            sendSweetAlert(
                session = session,
                title = "Error",
                text = 'Type the input file name.',
                type = "error"
            )
        } else {
            input_file_name <- isolate(input$input_data_file_name)
            if(!file.exists(paste0('data/',input_file_name))){
                sendSweetAlert(
                    session = session,
                    title = "Error",
                    text = paste0('The file "',input_file_name,'" does not exist in the ./data/ directory.'),
                    type = "error"
                )
            } else {
                shinyjs::disable("panelA")
                withProgress(message = 'Importing data', value = 0.5,{
                    input_table <<- fread(paste0('data/',input_file_name),sep = '\t',header = T)  # A global variable
                    input_table <<- data.frame(input_table)
                    if(length(dim(input_table)) >= 2  & dim(input_table)[2] == 1){
                        input_table <<- fread(paste0('data/',input_file_name),sep = ',',header = T)
                        input_table <<- data.frame(input_table)
                    }
                    setProgress(1, detail = paste0("finished"))
                })
                if(length(dim(input_table)) < 2){
                    sendSweetAlert(
                        session = session,
                        title = "Error",
                        text = 'File format is not valid.',
                        type = "error"
                    )
                    shinyjs::enable("panelA")
                } else {
                    required_col_names <- Standard_columns
                    if(length(intersect(colnames(input_table),required_col_names)) != length(required_col_names)){
                        sendSweetAlert(
                            session = session,
                            title = "Error",
                            text = 'Some required columns are missing in the input data.',
                            type = "error"
                        )
                        shinyjs::enable("panelA")
                    } else {
                        
                        shinyjs::show("success_for_input_file_name")
                        output$success_for_input_file_name <- renderText({paste0('The input data imported successfully.')})
                        
                        shinyjs::disable("panelA")
                        shinyjs::show("panelB")
                    }
                }
            }
        } 
    })
  
            
  

    observeEvent(input$continue_with_preprocessing,{
        shinyjs::disable("panelB")
              
        # Clearing the contents of "result" and "output" folders ------------------
        setwd('result')
            unlink(list.files(pattern = "\\.*$"),recursive = TRUE)
        setwd('..')
        setwd('output')
            setwd('signatures')
                setwd('5_mer')
                    unlink(list.files(pattern = "\\.*$"),recursive = TRUE)
                setwd('..')
                setwd('3_mer')
                    unlink(list.files(pattern = "\\.*$"),recursive = TRUE)
                setwd('..')
            setwd('..')
            
            setwd('clustering')
                unlink(list.files(pattern = "\\.*$"),recursive = TRUE)
            setwd('..')
        setwd('..')
        # -------------------------------------------------------------------------
        
    
    
    number_of_cpu_cores <<- isolate(input$CPU_cores) # a global variable
    
    
    
    
    
    withProgress(message = 'Preprocessing input data', value = 0,{
      
        setProgress(0.1, detail = 'Preparing the input table...')
        
        required_col_names <- Standard_columns
        input_table <<- input_table[,required_col_names]
        invalid_rows <- c()
        invalid_rows <- c(invalid_rows,    which(input_table[,'reference'] == '-' | input_table[,'mutated_to'] == '-')    )
        if(length(invalid_rows) !=0 ){input_table <<- input_table[-invalid_rows,]}    # remove
        input_table <<- unique(input_table)
        
        write.table(input_table,'result/input_table.csv',sep = ',',col.names = T,row.names = F)
        
        # Now the input_table is ready for counting mutations:
        input_table_file_name <- 'input_table'
        save_M5mer <- T
        M5mer_file_name <- 'M5mer'
        save_M3mer <- T
        M3mer_file_name <- 'M3mer'
        source('src/CountMutations.R',local = TRUE) # After this step, counts of 3mer and 5mer are placed in "result" folder...
    })
    
    
    shinyjs::hide("preprocessing_ui")    
    shinyjs::show("processing_ui")
            
  }) 
            
            
  observeEvent(input$skip_preprocessing,{
    shinyjs::disable("panelA")
    shinyjs::show("panelC")
  })
  
  
  
  observeEvent(input$continue,{
    shinyjs::disable("panelC")
   
    setwd('output')
    if(input$keep_or_not_signatures != 4){      # This means that at least one of the folders inside output folder will be deleted
      setwd('signatures')
      if(input$keep_or_not_signatures != 3){  # This means that 5_mer folder must be deleted 
        setwd('5_mer')
            unlink(list.files(pattern = "\\.*$"),recursive = TRUE)
        setwd('..')}
      if(input$keep_or_not_signatures != 2){  # This means that 3_mer folder must be deleted 
        setwd('3_mer')
            unlink(list.files(pattern = "\\.*$"),recursive = TRUE)
        setwd('..')}
      setwd('..')}
    
    if(input$keep_or_not_clustering != 2){
      setwd('clustering')
          unlink(list.files(pattern = "\\.*$"),recursive = TRUE)
      setwd('..')}
    setwd('..')
    
    number_of_cpu_cores <<- isolate(input$CPU_cores)
    
    shinyjs::hide("preprocessing_ui")    
    shinyjs::show("processing_ui")
  })

  
  
  
  
  # ------------------------------------------------------------------------------------------------------------------
  # 3-mer tab --------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------------
  
  observeEvent(input$start_3mer,{
      confirmSweetAlert(
          session = session,
          inputId = "confirm_start_3mer",
          type = "warning",
          title = "The previous results for 3-mer signatures will be deleted. Do you want to continue?",
          btn_labels = c("Cancel", "Continue"),
          danger_mode = TRUE
      )
  })
  
  observeEvent(input$confirm_start_3mer, {
      if (isTRUE(input$confirm_start_3mer)) {
          
          shinyjs::disable('panel_3mer_options')

          # first, delete the contents or output/signatures/3_mer
          setwd('output/signatures/3_mer')
          unlink(list.files(pattern = "\\.*$"),recursive = TRUE)
          setwd('../../..')
          
          
          accuracy <- input$accuracy_3mer
          
          if(accuracy == 'Moderate'){
              NMF_iters <- 1e4
              NMF_total_max <- 5e4
              NMF_conv <- 1e-4
              
              Boot_total_max <- 300
              Boot_conv <- 1e-1
          } else if(accuracy == 'High'){
              NMF_iters <- 1e4
              NMF_total_max <- 1e5
              NMF_conv <- 5e-5
              
              Boot_total_max <- 540
              Boot_conv <- 5e-2
          } else if(accuracy == 'Very High'){
              NMF_iters <- 1e4
              NMF_total_max <- 5e5
              NMF_conv <- 1e-5
              
              Boot_total_max <- 780
              Boot_conv <- 1e-2
          }
          
          NMF_max_epoches <- ceiling(NMF_total_max/NMF_iters)
          Boot_iters <- max(20,number_of_cpu_cores)
          Boot_max_epoches <- ceiling(Boot_total_max/Boot_iters)  
          
          start_N <- 1
          Max_N <- isolate(input$NMF_Max_N_3mer)
          
          k_mer <- 3
          file_name <- 'M3mer'
          destination_folder <- paste0("output/signatures/",as.character(k_mer),"_mer/")
          withProgress(message = 'Extracting mutational signatures', value = 0,{source('src/DecipherSignatures.R',local = TRUE)})
          
          sendSweetAlert(
              session = session,
              title = "Finished",
              text = 'See the "Results" tab.',
              type = "success"
          )
          shinyjs::enable('panel_3mer_options')
      }
  }, ignoreNULL = TRUE)
  
  
  plot_eval_diagram_3mer <- function(){
      k_mer <- 3
      destination_folder <- paste0("output/signatures/",as.character(k_mer),"_mer/")
      e <- read.table(paste0(destination_folder,"Evaluation.txt"))
      n <- e[,1]
      repro <- e[,2]
      frobe <- e[,3]
      frobe <- frobe/max(frobe)
      
      par(mar=c(5, 5, 4, 6) + 0.5)
      
      ymin <- 0.1
      
      plot(n, repro, pch=20, axes=FALSE, ylim=c(0,1),xlim = c(n[1],n[length(n)]+0.25), xlab="", ylab="",
           type="p",col="red", main="Evaluation for N",frame.plot = FALSE)
      
      grid(lwd = 2)
      
      par(new=TRUE)
      
      plot(n, frobe, pch=20,  xlab="", ylab="", ylim=c(0,1),xlim = c(n[1],n[length(n)]+0.25),
           axes=FALSE, type="p", col="blue",frame.plot = FALSE)
      
      mtext("Relative Frobenius Reconstruction Error",side=4,col="blue",line=2.5)
      axis(4,lwd = 2, ylim=range(frobe), col="blue",col.axis="blue",las=1)
      
      par(new=TRUE)
      
      plot(n, repro, pch=20, axes=FALSE, ylim=c(0,1),xlim = c(n[1],n[length(n)]+0.25), xlab="", ylab="",
           type="p",col="red", main="Evaluation for N",frame.plot = FALSE)
      lines(n, repro,col='red')
      axis(2, lwd = 2,ylim=range(repro),col="red",col.axis="red",las=1)
      mtext("Signatures Reproducibility",side=2,col ="red",line=3.75)
      
      axis(1,(n[1]-1):(n[length(n)]+1))
      mtext("Number of mutational signatures",side=1,col="black",line=2.5)
  }    
  output$download_3mer_eval <- downloadHandler(
          filename = function() { paste0('Evaluation_diagram(for 3-mer signatures).pdf') },
          content = function(file) {
              withProgress(message = 'Downloading...', value = 0.5,{
                  
                  k_mer <- 3
                  destination_folder <- paste0("output/signatures/",as.character(k_mer),"_mer/")
                  if( ! file.exists(paste0(destination_folder,"Evaluation.txt"))){
                      sendSweetAlert(
                          session = session,
                          title = "Error",
                          text = paste0('Evaluation file does not exist. Make sure that the required file is inside ./output/',
                                        k_mer,'_mer/ directory.'),
                          type = "error"
                      )
                  } else {
                      ggsave(file, plot = plot_eval_diagram_3mer(), device = "pdf")
                      setProgress(1)
                  } 
              })
          }
  )
  
  
  
  
  
  plot_sigs_3mer <- function(selected_N){
      k_mer <- 3
      destination_folder <- paste0("output/signatures/",as.character(k_mer),"_mer/")
      source('src/Plot3merSignatures.R',local = TRUE)
      plt <- plot_signatures_3mer(selected_N)
      for(p in plt){
          plot(p)
      }
  }
  output$download_3mer_signatures <- downloadHandler(
      filename = function() { paste0('3mer_Signatures_(N=',as.numeric(input$num_3mer_sigs_to_download),').pdf') },
      content = function(file) {
          withProgress(message = 'Downloading...', value = 0.5,{
              
              k_mer <- 3
              destination_folder <- paste0("output/signatures/",as.character(k_mer),"_mer/")
              selected_N <- as.numeric(input$num_3mer_sigs_to_download)
              if( ! file.exists(paste0(destination_folder,'3mer_Signatures_(N=',selected_N,').tsv'))){
                  sendSweetAlert(
                      session = session,
                      title = "Error",
                      text = paste0('Signature file for N = ',selected_N,
                                    ' does not exist. Make sure that the required file is inside ./output/',k_mer,'_mer/ directory.'),
                      type = "error"
                  )
              } else {
                  No_of_core_muts <- dim(fread(paste0(destination_folder,'3mer_Signatures_(N=',selected_N,').tsv')))[1]/16
                  ggsave(file, plot = plot_sigs_3mer(selected_N), device = "pdf", width = 3*No_of_core_muts, height = 5)
                  setProgress(1)
              }
          })
      }
  )

  output$download_files_3mer <- downloadHandler(
      filename = function() {
          paste("3mer_signature_analysis_results", "zip", sep=".")
      },
      content = function(fname) {
          k_mer <- 3
          if(length(list.files(paste0('output/signatures/',k_mer,'_mer/'))) == 0){
              sendSweetAlert(
                  session = session,
                  title = "Error",
                  text = paste0('directory ./output/signatures/',k_mer,'_mer/ is empty.'),
                  type = "error"
              )
          } else {
              setwd(paste0('output/signatures/',k_mer,'_mer/'))
              fs <- list.files()
              zip(zipfile=fname, files=fs)
              setwd('../../..')
          }
      },
      contentType = "application/zip"
  )

  
  
  # ------------------------------------------------------------------------------------------------------------------
  # 5-mer tab --------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------------
  observeEvent(input$start_5mer,{
    
    No_of_selected_motifs <- length(input$C_A_for_sig)+
                             length(input$C_G_for_sig)+
                             length(input$C_T_for_sig)+
                             length(input$T_A_for_sig)+
                             length(input$T_C_for_sig)+
                             length(input$T_G_for_sig)
    
    if(No_of_selected_motifs == 0){
      
      shinyjs::show("select_motif_error_5mer")
      output$select_motif_error_5mer <- renderText({paste("Select at least one 3-mer motif.")})
    
      } else if(No_of_selected_motifs > 10) {
        
        shinyjs::show("select_motif_error_5mer")
        output$select_motif_error_5mer <- renderText({paste("Number of selected 3-mer motifs can be at most 10.")})      
      
      } else {
        
          confirmSweetAlert(
              session = session,
              inputId = "confirm_start_5mer",
              type = "warning",
              title = "The previous results for 5-mer signatures will be deleted. Do you want to continue?",
              btn_labels = c("Cancel", "Continue"),
              danger_mode = TRUE
          )
      
      }
  })
  
  observeEvent(input$confirm_start_5mer, {
      if (isTRUE(input$confirm_start_5mer)) {
          
          shinyjs::hide('select_motif_error_5mer')
          shinyjs::disable('panel_5mer_options')
          
          
          # first, delete the contents or output/signatures/3_mer
          setwd('output/signatures/5_mer')
          unlink(list.files(pattern = "\\.*$"),recursive = TRUE)
          setwd('../../..')

          
          selected_3mer_motifs <- c(as.numeric(input$C_A_for_sig), 
                                    as.numeric(input$C_G_for_sig),
                                    as.numeric(input$C_T_for_sig),
                                    as.numeric(input$T_A_for_sig),
                                    as.numeric(input$T_C_for_sig),
                                    as.numeric(input$T_G_for_sig))
          
          write.table(selected_3mer_motifs,file = 'result/selected_3_mer_motifs.txt',row.names = F,col.names = F)
          
          withProgress(message = 'Preprocess data', value = 0,{source('src/MakeSelectedM5mer.R',local = TRUE)})
          
          
          accuracy <- input$accuracy_5mer
          
          if(accuracy == 'Moderate'){
              NMF_iters <- 1e4
              NMF_total_max <- 5e4
              NMF_conv <- 1e-4
              
              Boot_total_max <- 300
              Boot_conv <- 1e-1
          } else if(accuracy == 'High'){
              NMF_iters <- 1e4
              NMF_total_max <- 1e5
              NMF_conv <- 5e-5
              
              Boot_total_max <- 540
              Boot_conv <- 5e-2
          } else if(accuracy == 'Very High'){
              NMF_iters <- 1e4
              NMF_total_max <- 5e5
              NMF_conv <- 1e-5
              
              Boot_total_max <- 780
              Boot_conv <- 1e-2
          }
          
          NMF_max_epoches <- ceiling(NMF_total_max/NMF_iters)
          Boot_iters <- max(20,number_of_cpu_cores)
          Boot_max_epoches <- ceiling(Boot_total_max/Boot_iters)  
          
          start_N <- 1
          Max_N <- isolate(input$NMF_Max_N_5mer)
          
          k_mer <- 5
          file_name <- 'selectedM5mer'
          destination_folder <- paste0("output/signatures/",as.character(k_mer),"_mer/")
          withProgress(message = 'Extracting mutational signatures', value = 0,{source('src/DecipherSignatures.R',local = TRUE)})
          
          sendSweetAlert(
              session = session,
              title = "Finished",
              text = 'See the "Results" tab.',
              type = "success"
          )
          shinyjs::enable('panel_5mer_options')
      }
  }, ignoreNULL = TRUE)
  
  
  plot_eval_diagram_5mer <- function(){
      k_mer <- 5
      destination_folder <- paste0("output/signatures/",as.character(k_mer),"_mer/")
      e <- read.table(paste0(destination_folder,"Evaluation.txt"))
      n <- e[,1]
      repro <- e[,2]
      frobe <- e[,3]
      frobe <- frobe/max(frobe)
      
      par(mar=c(5, 5, 4, 6) + 0.5)
      
      ymin <- 0.1
      
      plot(n, repro, pch=20, axes=FALSE, ylim=c(0,1),xlim = c(n[1],n[length(n)]+0.25), xlab="", ylab="",
           type="p",col="red", main="Evaluation for N",frame.plot = FALSE)
      
      grid(lwd = 2)
      
      par(new=TRUE)
      
      plot(n, frobe, pch=20,  xlab="", ylab="", ylim=c(0,1),xlim = c(n[1],n[length(n)]+0.25),
           axes=FALSE, type="p", col="blue",frame.plot = FALSE)
      
      mtext("Relative Frobenius Reconstruction Error",side=4,col="blue",line=2.5)
      axis(4,lwd = 2, ylim=range(frobe), col="blue",col.axis="blue",las=1)
      
      par(new=TRUE)
      
      plot(n, repro, pch=20, axes=FALSE, ylim=c(0,1),xlim = c(n[1],n[length(n)]+0.25), xlab="", ylab="",
           type="p",col="red", main="Evaluation for N",frame.plot = FALSE)
      lines(n, repro,col='red')
      axis(2, lwd = 2,ylim=range(repro),col="red",col.axis="red",las=1)
      mtext("Signatures Reproducibility",side=2,col ="red",line=3.75)
      
      axis(1,(n[1]-1):(n[length(n)]+1))
      mtext("Number of mutational signatures",side=1,col="black",line=2.5)
  }
  output$download_5mer_eval <- downloadHandler(
      filename = function() { paste0('Evaluation_diagram(for 5-mer signatures).pdf') },
      content = function(file) {
          withProgress(message = 'Downloading...', value = 0.5,{
              
              k_mer <- 5
              destination_folder <- paste0("output/signatures/",as.character(k_mer),"_mer/")
              if( ! file.exists(paste0(destination_folder,"Evaluation.txt"))){
                  sendSweetAlert(
                      session = session,
                      title = "Error",
                      text = paste0('Evaluation file does not exist. Make sure that the required file is inside ./output/',
                                    k_mer,'_mer/ directory.'),
                      type = "error"
                  )
              } else {
                  ggsave(file, plot = plot_eval_diagram_5mer(), device = "pdf")
                  setProgress(1)
              } 
          })
      }
  )
  
  
  
  
  plot_sigs_5mer <- function(selected_N){
      k_mer <- 5
      destination_folder <- paste0("output/signatures/",as.character(k_mer),"_mer/")
      source('src/Plot5merSignatures.R',local = TRUE)
      plt <- plot_signatures_5mer(selected_N)
      for(p in plt){
          plot(p)
      }
  }
  output$download_5mer_signatures <- downloadHandler(
      filename = function() { paste0('5mer_Signatures_(N=',as.numeric(input$num_5mer_sigs_to_download),').pdf') },
      content = function(file) {
          withProgress(message = 'Downloading...', value = 0.5,{
              
              k_mer <- 5
              destination_folder <- paste0("output/signatures/",as.character(k_mer),"_mer/")
              selected_N <- as.numeric(input$num_5mer_sigs_to_download)
              
              if( ! file.exists(paste0(destination_folder,'5mer_Signatures_(N=',selected_N,').tsv'))){
                  sendSweetAlert(
                      session = session,
                      title = "Error",
                      text = paste0('Signature file for N = ',selected_N,
                                    ' does not exist. Make sure that the required file is inside ./output/',k_mer,'_mer/ directory.'),
                      type = "error"
                  )
              } else {
                  No_of_core_muts <- dim(fread(paste0(destination_folder,'5mer_Signatures_(N=',selected_N,').tsv')))[1]/16
                  ggsave(file, plot = plot_sigs_5mer(selected_N), device = "pdf", width = 3*No_of_core_muts, height = 5)
                  setProgress(1)
              }
          })
      }
  )
  
  
  output$download_files_5mer <- downloadHandler(
      filename = function() {
          paste("5mer_signature_analysis_results", "zip", sep=".")
      },
      content = function(fname) {
          k_mer <- 5
          if(length(list.files(paste0('output/signatures/',k_mer,'_mer/'))) == 0){
              sendSweetAlert(
                  session = session,
                  title = "Error",
                  text = paste0('directory ./output/signatures/',k_mer,'_mer/ is empty.'),
                  type = "error"
              )
          } else {
              setwd(paste0('output/signatures/',k_mer,'_mer/'))
              fs <- list.files()
              zip(zipfile=fname, files=fs)
              setwd('../../..')
          }
      },
      contentType = "application/zip"
  )
  
  
  
  
  
  
  
  # ------------------------------------------------------------------------------------------------------------------
  # Clustering tab ---------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------------
  observe({
    changed_option <- input$clustering_options_1
    shinyjs::hide("select_motif_error_3mer_clst")
    shinyjs::hide("select_motif_error_5mer_clst")
    shinyjs::hide("error_sig_proportion_3mer_clst")
    shinyjs::hide("error_sig_proportion_5mer_clst")
  })
  
  table_5mer_for_clst <- rhandsontable(initial_table_5mer_motifs_for_clst,useTypes=T,stretchH="all",width=600)%>%
    hot_col(col = "Selected",  halign = "htCenter")%>%
    hot_col(col = "Flanking1", halign = "htCenter")%>%
    hot_col(col = "Flanking2", halign = "htCenter")%>%
    hot_col(col = "Mutation",  halign = "htCenter")%>%
    hot_col(col = "Flanking3", halign = "htCenter")%>%
    hot_col(col = "Flanking4", halign = "htCenter")
  output$select_5mer_for_clst <- renderRHandsontable({table_5mer_for_clst})
  
  table_3mer_for_clst <- rhandsontable(initial_table_3mer_motifs_for_clst,useTypes=T,stretchH="all",width=410)%>%
    hot_col(col = "Selected",  halign = "htCenter")%>%
    hot_col(col = "Flanking1", halign = "htCenter")%>%
    hot_col(col = "Mutation",  halign = "htCenter")%>%
    hot_col(col = "Flanking2", halign = "htCenter")
  output$select_3mer_for_clst <- renderRHandsontable({table_3mer_for_clst})

  observe({if(!is.null(input$select_5mer_for_clst)){table_5mer_motifs_for_clst <<- hot_to_r(input$select_5mer_for_clst)}})
  observe({if(!is.null(input$select_3mer_for_clst)){table_3mer_motifs_for_clst <<- hot_to_r(input$select_3mer_for_clst)}})
  
  
  observeEvent(input$start_clustering,{
    shinyjs::hide("select_motif_error_3mer_clst")
    shinyjs::hide("select_motif_error_5mer_clst")
    shinyjs::hide("error_sig_proportion_3mer_clst")
    shinyjs::hide("error_sig_proportion_5mer_clst")
    
    if(input$clustering_options_1 == 'Cluster based on 3-mer mutational motifs') {
        selected_indxs <- which(table_3mer_motifs_for_clst$Selected == TRUE)
        if(length(selected_indxs) == 0) {
          shinyjs::show("select_motif_error_3mer_clst")
          output$select_motif_error_3mer_clst <- renderText({paste("Select at least one 3-mer motif.")})
        } else if(!all(table_3mer_motifs_for_clst[selected_indxs,-1] != 'NA')) {
          shinyjs::show("select_motif_error_3mer_clst")
          output$select_motif_error_3mer_clst <- renderText({
            paste('The rows which are checked as "selected" must be filled completely.')})
        } else {
          shinyjs::hide("select_motif_error_3mer_clst")
          shinyjs::disable('panel_clustering_options_1')
          nuc <- c(0:3)
          names(nuc) <- c('A','C','G','T')
          mut <- c(0:5)
          names(mut) <- c('C > A','C > G','C > T','T > A','T > C','T > G')
          selected_3mer_motifs <- table_3mer_motifs_for_clst[selected_indxs,-1]
          fl1 <- as.numeric(nuc[selected_3mer_motifs$Flanking1])
          mu <- as.numeric(mut[selected_3mer_motifs$Mutation])
          fl2 <- as.numeric(nuc[selected_3mer_motifs$Flanking2])
          selected_3mer_motifs_for_clst <- 16*mu+4*fl1+fl2+1
          write.table(selected_3mer_motifs_for_clst,file = 'result/selected_3mer_motifs_for_clst.txt',row.names=F,col.names=F)
          k_mer <- 3
          motif_or_signature <- 'motif'
          if(input$count_or_proportion_3mer == 'Cluster based on "count" values of 3-mer mutational motifs') {
            count_or_proportion <- 'count'
          } else if(input$count_or_proportion_3mer == 'Cluster based on "proportions" of 3-mer mutational motifs') {
            count_or_proportion <- 'proportion'
          }
          withProgress(message = 'Preparing data for clustering', value = 0.5,{source('src/PrepareForClustering.R',local = TRUE)})
          withProgress(message = 'Clustering', value = 0,{source('src/Clustering.R',local = TRUE)})
        }

    } else if(input$clustering_options_1 == 'Cluster based on 5-mer mutational motifs') {
        selected_indxs <- which(table_5mer_motifs_for_clst$Selected == TRUE)
        if(length(selected_indxs) == 0) {
          shinyjs::show("select_motif_error_5mer_clst")
          output$select_motif_error_5mer_clst <- renderText({paste("Select at least one 5-mer motif.")})
        } else if(!all(table_5mer_motifs_for_clst[selected_indxs,-1] != 'NA')) {
            shinyjs::show("select_motif_error_5mer_clst")
            output$select_motif_error_5mer_clst <- renderText({
              paste('The rows which are checked as "selected" must be filled completely.')})
        } else {
          shinyjs::hide("select_motif_error_5mer_clst")
          shinyjs::disable('panel_clustering_options_1')
          nuc <- c(0:3)
          names(nuc) <- c('A','C','G','T')
          mut <- c(0:5)
          names(mut) <- c('C > A','C > G','C > T','T > A','T > C','T > G')
          selected_5mer_motifs <- table_5mer_motifs_for_clst[selected_indxs,-1]
          fl1 <- as.numeric(nuc[selected_5mer_motifs$Flanking1])
          fl2 <- as.numeric(nuc[selected_5mer_motifs$Flanking2])
          mu <- as.numeric(mut[selected_5mer_motifs$Mutation])
          fl3 <- as.numeric(nuc[selected_5mer_motifs$Flanking3])
          fl4 <- as.numeric(nuc[selected_5mer_motifs$Flanking4])
          selected_5mer_motifs_for_clst <- 16*(16*mu+4*fl2+fl3)+4*fl1+fl4+1
          write.table(selected_5mer_motifs_for_clst,file = 'result/selected_5mer_motifs_for_clst.txt',row.names=F,col.names=F)
          k_mer <- 5
          motif_or_signature <- 'motif'
          if(input$count_or_proportion_5mer == 'Cluster based on "count" values of 5-mer mutational motifs') {
            count_or_proportion <- 'count'
          } else if(input$count_or_proportion_5mer == 'Cluster based on "proportions" of 5-mer mutational motifs') {
            count_or_proportion <- 'proportion'
          }
          withProgress(message = 'Preparing data for clustering', value = 0.5,{source('src/PrepareForClustering.R',local = TRUE)})
          withProgress(message = 'Clustering', value = 0,{source('src/Clustering.R',local = TRUE)})
        }


        
    } else if(input$clustering_options_1 == 'Cluster based on proportions of 3-mer signatures' |
              input$clustering_options_1 == 'Cluster based on proportions of 5-mer signatures') {
      
      option <- input$clustering_options_1
      k_mer <- as.numeric(substr(option,33,33))

      files <- list.files(path=paste0('output/signatures/',as.character(k_mer),'_mer/'),pattern = "\\.*$")
      if(length(files) == 0)
      {
        shinyjs::show(paste0("error_sig_proportion_",as.character(k_mer),"mer_clst"))
        output[[paste0("error_sig_proportion_",as.character(k_mer),"mer_clst")]] <- renderText({
          paste0('Error: The output/signatures/',as.character(k_mer),'_mer/ folder is empty.')})
      } else {
        
        shinyjs::hide(paste0("error_sig_proportion_",as.character(k_mer),"mer_clst"))
        shinyjs::disable('panel_clustering_options_1')
        
        motif_or_signature <- 'signature'
        count_or_proportion <- 'proportion'
          
        withProgress(message = 'Preparing data for clustering', value = 0.5,{source('src/PrepareForClustering.R',local = TRUE)})
        withProgress(message = 'Clustering', value = 0,{source('src/Clustering.R',local = TRUE)})
      }
    }
    
    shinyjs::show('final_message_clustering')
  })
  
}

shinyApp(ui, server)


