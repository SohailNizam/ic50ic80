bnabs <- list.files(paste("./docker_output"))
bnabs <- bnabs[!bnabs %in% c("dh270.1", "dh270.5", "dh270.6", "vrc38.01")]



get_cors <- function(bnabs, cor_type = "pearson", use_cv_learner = TRUE){
  bnab_name <- c()
  ic50_cors <- c()
  ic80_cors <- c()
  problematic_bnabs <- c()

  for(bnab in bnabs){
    
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
    
  
    ic50_preds <- ic50_results$SL.predict
    ic80_preds <- ic80_results$SL.predict
    
    csv_file <- read.csv(Sys.glob(paste("./docker_output/", bnab, "/slapnap", "*", ".csv", sep = "")))
    ic50_obs <- csv_file$ic50[!is.na(csv_file$ic50)]
    ic80_obs <- csv_file$ic80[!is.na(csv_file$ic80)]
    
    
    print(length(ic50_preds))
    print(length(ic50_obs))
    print("")
    print(length(ic80_preds))
    print(length(ic80_obs))
    if(length(ic50_preds) != length(ic50_obs)){
      problematic_bnabs <- c(problematic_bnabs, bnab)
      next
      }
    if(length(ic80_preds) != length(ic80_obs)){
      problematic_bnabs <- c(problematic_bnabs, bnab)
      next
      }
    
    ic50_cor <- cor(ic50_preds, ic50_obs, method = cor_type)
    ic80_cor <- cor(ic80_preds, ic80_obs, method = cor_type)
    
    bnab_name <- c(bnab_name, bnab)
    ic50_cors <- c(ic50_cors, ic50_cor)
    ic80_cors <- c(ic80_cors, ic80_cor)
    
    print(paste(bnab, "done"))
    

  }
  print(paste("Problematic bnabs", problematic_bnabs))
  return(data.frame(bnab_name, ic50_cors, ic80_cors))
}


cvlearner_cors <- get_cors(bnabs)
'''
"Problematic bnabs 2f5"      "Problematic bnabs 2g12"     "Problematic bnabs 4e10"    
 [4] "Problematic bnabs b12"      "Problematic bnabs pg16"     "Problematic bnabs pg9"     
[7] "Problematic bnabs pgt121"   "Problematic bnabs pgt128"   "Problematic bnabs pgt135"  
[10] "Problematic bnabs pgt145"   "Problematic bnabs vrc-pg04" "Problematic bnabs vrc01"
'''

learner_cors <- get_cors(bnabs)

'''
 [1] "Problematic bnabs 2f5"      "Problematic bnabs 2g12"     "Problematic bnabs 4e10"    
[4] "Problematic bnabs b12"      "Problematic bnabs pg16"     "Problematic bnabs pg9"     
[7] "Problematic bnabs pgt121"   "Problematic bnabs pgt128"   "Problematic bnabs pgt135"  
[10] "Problematic bnabs pgt145"   "Problematic bnabs vrc-pg04" "Problematic bnabs vrc01"  
'''
