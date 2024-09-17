##########################
#Olink analysis
#August-2024
#Ziyi Li
#CARD NIH
#########################

#### PACKAGES ######################################################################################
package_list = c('OlinkAnalyze','ggplot2', 
                 'clusterProfiler','DOSE','enrichplot','org.Hs.eg.db',
                 'dplyr')
cat("INFO: Loading required packages\n      ")
cat(paste(package_list, collapse='\n      ')); cat('\n')

defaultW <- getOption("warn"); options(warn = -1)   # Temporarily disable warnings for quiet loading
if(all((lapply(package_list, require, character.only=TRUE)))) {
  cat("INFO: All packages successfully loaded\n")
} else {
  cat("ERROR: One or more packages not available. Are you running this within the container?\n")
}
options(warn = defaultW)    # Turn warnings back on# Load OlinkAnalyze

explore_npx <- read_NPX("/Users/liz36/Documents/Brain_Sample/olink/Q-13258_Qi_Extended_NPX_2024-06-06.parquet")
my_NPX_data <- read_NPX(filename = "/Users/liz36/Documents/Brain_Sample/olink/_Users_haoy4_Documents_Research_Project_Braintissue_olink_result_202406201331.csv")
explore_npx %>%Olink_dist_plot() +
  +   theme(axis.text.x = element_blank()) 

explore_npx %>% 
  filter(SampleType=='SAMPLE')%>%# For this example only plotting one panel.
  olink_dist_plot(color_g = "SampleQC") +
  theme(axis.text.x = element_blank()) # Due to the number of samples one can remove the text or rotate it

explore_npx %>% 
  filter(!str_detect(SampleID, 'CONTROL_SAMPLE'),
         Panel == 'Olink Inflammation') %>% 
  olink_qc_plot(color_g = "AssayQC")   


olink_lod(my_NPX_data, lod_method = "NCLOD")
data("npx_data2")
data(npx_data1)
olink_lod(npx_data1, lod_method = "NCLOD")



