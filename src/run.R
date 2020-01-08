

library(configr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

user_config_path <- args[1]
default_config_path <- args[2]
source_codes_path <- args[3]

user_config <- read.config(user_config_path)[['config']]

config <- read.config(default_config_path)[['default config']]

# overwrite default config
for(o in names(user_config)){
    config[[o]] <- user_config[[o]]
}


inspect_config <- function(config, option){
    if(option %in% names(config)){
        return(config[[option]])
    } else {
        return('none')
    }
}


#--------------------------------------------------------------------------------------------------------
if(inspect_config(config, 'preprocessing') == 'yes'){
    source(paste0(source_codes_path,'/CountMutations.R'))
    CountMutations(input_file_path = config[['input_file']],
                   output_dir = paste0(config[['output_dir']],'/','preprocessing'))
} else {
    # check for preprocessing files
    cat("Assuming that preprocessing has already been done\n")
}
#--------------------------------------------------------------------------------------------------------
decipher_signatures_if_required <- function(k_mer){
    if(inspect_config(config, paste0('infer_',k_mer,'mer_signatures')) == 'yes'){
        
        selected_motifs <- NULL
        if(k_mer == 5){
            if(inspect_config(config, config[['selected_3mer_motifs_for_5mer_signatures']]) == 'none'){
                stop()
            }
            selected_motifs <- gsub(' ','',strsplit(config[['selected_3mer_motifs_for_5mer_signatures']], ',')[[1]])
        }
        
        source(paste0(source_codes_path,'/DecipherSignatures.R'))
        DecipherSignatures(input_file = paste0(config[['output_dir']],'/','preprocessing','/','motif-matrix-',k_mer,'mer.csv'),
                           output_dir = paste0(config[['output_dir']],'/','infered_',k_mer,'mer_signatures'),
                           
                           k_mer = k_mer,
                           
                           selected_motifs = selected_motifs,
                           
                           N_min = as.numeric(config[[paste0('N_min_',k_mer,'mer')]]),
                           N_max = as.numeric(config[[paste0('N_max_',k_mer,'mer')]]),
                           
                           number_of_CPUs = as.numeric(config[[paste0('CPU_',k_mer,'mer')]]),
                           
                           NMF_iters = as.numeric(config[[paste0('nmf_iters_',k_mer,'mer')]]),
                           NMF_conv = as.numeric(config[[paste0('nmf_conv_',k_mer,'mer')]]),
                           NMF_total_max = as.numeric(config[[paste0('nmf_max_',k_mer,'mer')]]),
                           
                           Boot_iters = as.numeric(config[[paste0('boot_iters_',k_mer,'mer')]]),
                           Boot_conv = as.numeric(config[[paste0('boot_conv_',k_mer,'mer')]]),
                           Boot_total_max = as.numeric(config[[paste0('boot_max_',k_mer,'mer')]])
                           )
    }
}
decipher_signatures_if_required(k_mer = 3)
decipher_signatures_if_required(k_mer = 5)
#--------------------------------------------------------------------------------------------------------
cluster_based_on_signatures_if_required <- function(k_mer){
    if(inspect_config(config, paste0('cluster_samples_based_on_',k_mer,'mer_signatures')) == 'yes'){
        if(inspect_config(config, paste0('optimum_number_of_',k_mer,'mer_signatures_for_clustering')) == 'none'){
            stop(paste0('missing parameter -> "optimum_number_of_',k_mer,'mer_signatures_for_clustering"'),call. = F)
        }
        input_file <- paste0(config[['output_dir']],'/',
                             'infered_',k_mer,'mer_signatures','/',
                             'Exposures_to_',k_mer,'mer_Signatures_(N=',
                             config[[paste0('optimum_number_of_',k_mer,'mer_signatures_for_clustering')]],
                             ').tsv')
        if( ! file.exists(input_file)){
            stop(paste0('missing file -> "',input_file,'"'),call. = F)
        }
        source(paste0(source_codes_path,'/Clustering.R'))
        Clustering(input_file = input_file,
                   output_dir = paste0(config[['output_dir']],'/','cluster_samples_based_on_',k_mer,'mer_signatures'),
                   basis = 'signatures'
        )
    }
}
cluster_based_on_signatures_if_required(k_mer = 3)
cluster_based_on_signatures_if_required(k_mer = 5)
#--------------------------------------------------------------------------------------------------------
cluster_based_on_motifs_if_required <- function(k_mer){
    if(inspect_config(config, paste0('cluster_samples_based_on_',k_mer,'mer_motifs')) == 'yes'){
        if(inspect_config(config, paste0('selected_',k_mer,'mer_motifs_for_clustering')) == 'none'){
            stop(paste0('missing parameter -> "selected_',k_mer,'mer_motifs_for_clustering"'),call. = F)
        }
        input_file <- paste0(config[['output_dir']],'/','preprocessing','/',paste0('motif-matrix-',k_mer,'mer.csv'))
        if( ! file.exists(input_file)){
            stop(paste0('missing file -> "',input_file,'"'),call. = F)
        }
        source(paste0(source_codes_path,'/Clustering.R'))
        Clustering(input_file = input_file,
                   output_dir = paste0(config[['output_dir']],'/','cluster_samples_based_on_',k_mer,'mer_motifs'),
                   basis = 'motifs',
                   selected_motifs = gsub(' ','',strsplit(config[[paste0('selected_',k_mer,'mer_motifs_for_clustering')]], ',')[[1]])
        )
    }
}
cluster_based_on_motifs_if_required(k_mer = 3)
cluster_based_on_motifs_if_required(k_mer = 5)


cat('Finished\n')


