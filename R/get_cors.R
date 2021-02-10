#Use SuperLearner objects to get median correlation between obs and predicted ic50/ic80
#Within each epitope

#import packages
library(dplyr)

#get the bnabs based on directory names in ic50ic80/docker_output
bnabs <- list.files(paste("./docker_output"))
#remove bnabs whose containers didn't finish
bnabs <- bnabs[!bnabs %in% c("dh270.1", "dh270.5", "dh270.6", "vrc38.01")]
# TODO: add a tryCatch to the function to generalize this ^


#The main function
get_cors <- function(bnabs, cor_type = "pearson", use_cv_learner = TRUE){
  
  '''
  This function takes in a list of bnab names, finds their SL objects,
  and retrieves the observed and predicted ic50 and ic80 values. Correlations
  are calculated and all results outputed in a df.

  Input: bnabs (chr vector), cor_type (str, "pearson" or "spearman")
         use_cv_learner (bool, TRUE means use th cvlearner SL object)

  Output: dataframe (cols = bnab, ic50_cors, ic80_cors)
  '''
  
  # init empty vecs for the final df
  bnab_name <- c()
  ic50_cors <- c()
  ic80_cors <- c()


  # cycle through each bnab
  for(bnab in bnabs){
    
    # get either the cvlearner or learner SL objects
    if(use_cv_learner){
      ic50_results <- readRDS(Sys.glob(paste("./docker_output/", 
                                             bnab, "/cvlearner_ic50","*", ".rds", sep = "")))
      ic80_results <- readRDS(Sys.glob(paste("./docker_output/", 
                                             bnab, "/cvlearner_ic80","*", ".rds", sep = "")))
    }
    
    else{
      
      ic50_results <- readRDS(Sys.glob(paste("./docker_output/", 
                                             bnab, "/learner_ic50","*", ".rds", sep = "")))
      ic80_results <- readRDS(Sys.glob(paste("./docker_output/", 
                                             bnab, "/learner_ic80","*", ".rds", sep = "")))
    }
    
    # get the ic50/ic80 predictions from SL object
    ic50_preds <- ic50_results$SL.predict
    ic80_preds <- ic80_results$SL.predict
    
    # get the ic50/ic80 observed vals from SL object
    ic50_obs <- ic50_results$Y
    ic80_obs <- ic80_results$Y
    
    # calculate the correlations
    ic50_cor <- cor(ic50_preds, ic50_obs, method = cor_type)
    ic80_cor <- cor(ic80_preds, ic80_obs, method = cor_type)
    
    # update the vectors with this bnabs info
    bnab_name <- c(bnab_name, bnab)
    ic50_cors <- c(ic50_cors, ic50_cor)
    ic80_cors <- c(ic80_cors, ic80_cor)
  }
  
  return(data.frame(bnab_name, ic50_cors, ic80_cors))
}


# call our function for cvlearner and learner
cvlearner_cors <- get_cors(bnabs, use_cv_learner = TRUE)
learner_cors <- get_cors(bnabs, use_cv_learner = FALSE)

# change df to tibble with epitope variable
cvlearner_tib <- cvlearner_cors %>%
  as_tibble() %>%
  mutate(
    bnab_name = toupper(bnab_name),
    epitope = case_when(
      bnab_name %in% c("PG16", "PG9", "PGDM1400",
                  "PGT145", "VRC26.08", "VRC26.25") ~ "V1V2",
      bnab_name %in% c("10-1074", "10-996", "DH270.1",
                  "DH270.5", "DH270.6", "PGT121", "PGT128",
                  "PGT135", "VRC29.03", "VRC38.01", "2G12") ~ "V3",
      bnab_name %in% c("3BNC117", "B12", "CH01",
                  "HJ16", "NIH45-46", "VRC-CH31", "VRC-PG04",
                  "VRC01", "VRC03", "VRC07") ~ "CD4bs",
      bnab_name %in% c("PGT151", "VRC34.01") ~ "Fusion peptide",
      bnab_name %in% c("35O22", "8ANC195") ~ "Subunit interface",
      bnab_name %in% c("2F5", "4E10") ~ "MPER"
    ))

# change df to tibble with epitope variable
learner_tib <- learner_cors %>%
  as_tibble() %>%
  mutate(
    bnab_name = toupper(bnab_name),
    epitope = case_when(
      bnab_name %in% c("PG16", "PG9", "PGDM1400",
                       "PGT145", "VRC26.08", "VRC26.25") ~ "V1V2",
      bnab_name %in% c("10-1074", "10-996", "DH270.1",
                       "DH270.5", "DH270.6", "PGT121", "PGT128",
                       "PGT135", "VRC29.03", "VRC38.01", "2G12") ~ "V3",
      bnab_name %in% c("3BNC117", "B12", "CH01",
                       "HJ16", "NIH45-46", "VRC-CH31", "VRC-PG04",
                       "VRC01", "VRC03", "VRC07") ~ "CD4bs",
      bnab_name %in% c("PGT151", "VRC34.01") ~ "Fusion peptide",
      bnab_name %in% c("35O22", "8ANC195") ~ "Subunit interface",
      bnab_name %in% c("2F5", "4E10") ~ "MPER"
    ))


## Create tables with median correlation by epitope for ic50 and ic80 ##

#cvlearner results, grouped by epitope
cvlearner_med_grp <- cvlearner_tib %>%
  group_by(epitope) %>% 
  summarise(ic50_median=median(ic50_cors, na.rm = TRUE), 
            ic80_median=median(ic80_cors, na.rm = TRUE))

#cvlearner results, overall
cvlearner_med <- cvlearner_tib %>%
  summarise(ic50_median=median(ic50_cors, na.rm = TRUE), 
            ic80_median=median(ic80_cors, na.rm = TRUE))

#learner results, grouped by epitope
learner_med_grp <- learner_tib %>%
  group_by(epitope) %>% 
  summarise(ic50_median=median(ic50_cors, na.rm = TRUE), 
            ic80_median=median(ic80_cors, na.rm = TRUE))

#learner results, overall
learner_med <- learner_tib %>%
  summarise(ic50_median=median(ic50_cors, na.rm = TRUE), 
            ic80_median=median(ic80_cors, na.rm = TRUE))

#write tables to csv files
write.csv(cvlearner_med_grp, file = "./R_output/cvlearner_med_cors_grp.csv")
write.csv(learner_med_grp, file = "./R_output/learner_med_cors_grp.csv")
write.csv(cvlearner_med, file = "./R_output/cvlearner_med_cors.csv")
write.csv(learner_med, file = "./R_output/learner_med_cors.csv")
