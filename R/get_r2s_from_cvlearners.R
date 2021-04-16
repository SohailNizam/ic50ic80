#Functions from slapnap/code/00_utils.R to prepare a given cvlearner object
summary.myCV.SuperLearner <- function (object, obsWeights = NULL, method = NULL, opts, ...) {
  if ("env" %in% names(object)) {
    env = object$env
  }else {
    env = parent.frame()
  }
  
  is_sl <- "SuperLearner" %in% class(object)
  is_cvsl <- "CV.SuperLearner" %in% class(object)
  if(is_sl | is_cvsl){
    library.names <- object$libraryNames
    if(is_cvsl){
      V <- object$V
    }else{
      V <- length(object$validRows)
    }
    n <- length(object$SL.predict)
    if (is.null(obsWeights)) {
      obsWeights <- rep(1, length(object$Y))
    }
    
    if(is_cvsl){
      folds <- object$folds
    }else if(is_sl){
      folds <- object$validRows
    }
    
    if(is_cvsl){
      # only will use this if multiple learners selected
      SL.predict <- object$SL.predict
      # this will be only output if single learner used and opts$cvtune
      discreteSL.predict <- object$discreteSL.predict
      # only will use this if multiple learners selected
      library.predict <- object$library.predict
    }else if(is_sl){
      # in this case a single "default" learner was requested
      # so we can pull Z out from the object
      SL.predict <- object$Z[,1]
    }
    
    Y <- object$Y
    Risk.SL <- rep(NA, length = V)
    se.SL <- rep(NA, length = V)
    if(is_cvsl){
      Risk.dSL <- rep(NA, length = V)
      se.dSL <- rep(NA, length = V)
      Risk.library <- matrix(NA, nrow = length(library.names),
                             ncol = V)
      se.library <- matrix(NA, nrow = length(library.names),
                           ncol = V)
      rownames(Risk.library) <- library.names
    }
    if (!(all(Y %in% c(0,1)))) {
      for (ii in seq_len(V)) {
        Risk.SL[ii] <- mean(obsWeights[folds[[ii]]] * (Y[folds[[ii]]] -
                                                         SL.predict[folds[[ii]]])^2)
        if(is_cvsl){
          Risk.dSL[ii] <- mean(obsWeights[folds[[ii]]] * (Y[folds[[ii]]] -
                                                            discreteSL.predict[folds[[ii]]])^2)
          Risk.library[, ii] <- apply(library.predict[folds[[ii]],
                                                      , drop = FALSE], 2, function(x) mean(obsWeights[folds[[ii]]] *
                                                                                             (Y[folds[[ii]]] - x)^2))
        }
      }
      if_sl <- (Y - SL.predict)^2 - mean((Y - SL.predict)^2)
      if(is_cvsl){
        if_dsl <- (Y - discreteSL.predict)^2 - mean((Y - discreteSL.predict)^2)
        if_library <- apply(library.predict, 2, function(x){ (Y - x)^2 - mean((Y - x)^2) })
      }
      if_varY <- (Y - mean(Y))^2 - mean((Y - mean(Y))^2)
      get_log_se <- function(if_risk, if_varY, risk, varY,
                             n = length(if_risk)){
        grad <- matrix(c(1 / risk, - 1 /varY), nrow = 2)
        Sig <- cov(cbind(if_risk, if_varY))
        se_log <- t(grad) %*% Sig %*% grad
        return(se_log)
      }
      
      if(is_cvsl){
        if(length(opts$learners) > 1){
          se <- (1/sqrt(n)) * c(
            get_log_se(if_risk = if_sl, if_varY = if_varY, risk = mean(Risk.SL), varY = var(Y)),
            get_log_se(if_risk = if_dsl, if_varY = if_varY, risk = mean(Risk.dSL), varY = var(Y)),
            mapply(if1 = split(if_library, col(if_library)), risk = split(Risk.library, row(Risk.library)),
                   function(if1, risk){ get_log_se(if_risk = if1, if_varY = if_varY, risk = mean(risk), varY = var(Y))})
          )
          Table <- data.frame(Algorithm = c("Super Learner", "Discrete SL",
                                            library.names), Ave = c(1 - mean(Risk.SL)/var(Y), 1 - mean(Risk.dSL)/var(Y),
                                                                    apply(Risk.library, 1, function(x){ 1 - mean(x)/var(Y) })), log_se = se, Min = c(min(1 - Risk.SL/var(Y)),
                                                                                                                                                     min(1 - Risk.dSL/var(Y)), apply(Risk.library, 1, function(x){ min(1 - mean(x)/var(Y))})), Max = c(max(1 - Risk.SL/var(Y)),
                                                                                                                                                                                                                                                       max(1 - Risk.dSL/var(Y)), apply(Risk.library, 1, function(x){ max(1 - mean(x)/var(Y)) })))
          
        }else{
          se <- (1/sqrt(n)) * c(
            get_log_se(if_risk = if_dsl, if_varY = if_varY, risk = mean(Risk.dSL), varY = var(Y)),
            mapply(if1 = split(if_library, col(if_library)), risk = split(Risk.library, row(Risk.library)),
                   function(if1, risk){ get_log_se(if_risk = if1, if_varY = if_varY, risk = mean(risk), varY = var(Y))})
          )
          
          Table <- data.frame(Algorithm = c("Discrete SL",
                                            library.names), Ave = c(1 - mean(Risk.dSL)/var(Y),
                                                                    apply(Risk.library, 1, function(x){ 1 - mean(x)/var(Y) })), log_se = se,
                              Min = c(min(1 - Risk.dSL/var(Y)), apply(Risk.library, 1, function(x){ min(1 - mean(x)/var(Y))})),
                              Max = c(max(1 - Risk.dSL/var(Y)), apply(Risk.library, 1, function(x){ max(1 - mean(x)/var(Y)) })))
        }
      }else{
        se <- (1/sqrt(n)) * get_log_se(if_risk = if_sl, if_varY = if_varY, risk = mean(Risk.SL), varY = var(Y))
        Table <- data.frame(Algorithm = c(library.names[1]), Ave = c(1 - mean(Risk.SL)/var(Y)),
                            log_se = se,
                            Min = c(min(1 - Risk.SL/var(Y))),
                            Max = c(max(1 - Risk.SL/var(Y))))
      }
    }else {
      requireNamespace("cvAUC")
      for (ii in seq_len(V)) {
        sl_auc <- cvAUC::ci.cvAUC(predictions = SL.predict[folds[[ii]]],
                                  labels = Y[folds[[ii]]], folds = NULL)
        Risk.SL[ii] <- sl_auc$cvAUC
        se.SL[ii] <- sl_auc$se
        if(is_cvsl){
          dsl_auc <- cvAUC::ci.cvAUC(predictions = discreteSL.predict[folds[[ii]]],
                                     labels = Y[folds[[ii]]], folds = NULL)
          Risk.dSL[ii] <- dsl_auc$cvAUC
          se.dSL[ii] <- dsl_auc$se
          library_auc <- apply(library.predict[folds[[ii]], , drop = FALSE], 2, function(x){
            tmp <- cvAUC::ci.cvAUC(predictions = x, labels = Y[folds[[ii]]], folds = NULL)
            return(c(tmp$cvAUC, tmp$se))
          })
          Risk.library[,ii] <- library_auc[1,]
          se.library[,ii] <- library_auc[2,]
        }
      }
      if(is_cvsl){
        if(length(opts$learners) > 1){
          se <- c(mean(se.SL, na.rm = TRUE), mean(se.dSL, na.rm = TRUE),
                  rowMeans(se.library, na.rm = TRUE))
          Table <- data.frame(Algorithm = c("Super Learner", "Discrete SL",
                                            library.names), Ave = c(mean(Risk.SL), mean(Risk.dSL),
                                                                    apply(Risk.library, 1, mean)), se = se, Min = c(min(Risk.SL),
                                                                                                                    min(Risk.dSL), apply(Risk.library, 1, min)), Max = c(max(Risk.SL),
                                                                                                                                                                         max(Risk.dSL), apply(Risk.library, 1, max)))
        }else{
          se <- c(mean(se.dSL, na.rm = TRUE),
                  rowMeans(se.library, na.rm = TRUE))
          Table <- data.frame(Algorithm = c("Discrete SL",
                                            library.names), Ave = c(mean(Risk.dSL),
                                                                    apply(Risk.library, 1, mean)), se = se, Min = c(
                                                                      min(Risk.dSL), apply(Risk.library, 1, min)), Max = c(
                                                                        max(Risk.dSL), apply(Risk.library, 1, max)))
        }
      }else{
        se <- c(mean(se.SL, na.rm = TRUE))
        Table <- data.frame(Algorithm = c(library.names[1]),
                            Ave = c(mean(Risk.SL)), se = se,
                            Min = c(min(Risk.SL)),
                            Max = c(max(Risk.SL)))
      }
    }
    out <- list(call = object$call, method = method, V = V, Table = Table)
  }
  class(out) <- "summary.myCV.SuperLearner"
  return(out)
}
get_est_and_ci <- function(idx = 1, fit_list = NULL, Rsquared = FALSE, constant = qnorm(0.975)){
  cv_fit_table <- fit_list[[idx]]
  Mean <- cv_fit_table$Table$Ave
  if(Rsquared){
    se <- cv_fit_table$Table$log_se
    Lower <- 1 - exp( log(-Mean + 1) + constant * se)
    Upper <- 1 - exp( log(-Mean + 1) - constant * se)
  }else{
    se <- cv_fit_table$Table$se
    # put AUC CI on logit scale
    grad <- 1 / (Mean - Mean^2)
    logit_se <- sqrt(se^2 * grad^2)
    Lower <- plogis(qlogis(Mean) - constant * logit_se); Upper <- plogis(qlogis(Mean) + constant * logit_se)
  }
  return(list(est = Mean[1], ci = c(Lower[1], Upper[1])))
}

bnabs <- c("10-1074", "10-996", "2f5", "2g12", "35o22",
         "3bnc117", "4e10", "8anc195", "b12", "ch01", 
         "dh270.1", "dh270.5", "dh270.6", "hj16",
         "nih45-46", "pg16", "pg9", "pgdm1400", "pgt121", 
         "pgt128", "pgt135", "pgt145", "pgt151", "vrc-ch31",
        "vrc-pg04", "vrc01", "vrc03", "vrc07", "vrc26.08",
        "vrc26.25", "vrc29.03", "vrc34.01", "vrc38.01")



#A function to cycle through each nab and get ic50 and ic80 r2, ub, lb from a prepared cvlearner object
get_all_est_and_cis <- function(bnabs){
  
  #initialize a df with cols: method (ic50/ic80), bnab, r2,  cil, ciu, n
  cvr2_df <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(cvr2_df) <- c("method", "bnab", "r2", "cil", "ciu", "n")
  
  for(bnab in bnabs){
    
    #read in ic50 and ic80 cvlearner objects
    ic50_file <- Sys.glob(paste0("./new_docker_output/", bnab, "/cvlearner_ic50","*", ".rds"))
    if(length(ic50_file > 0)){ #proceed only if that file exists (some containers failed)
      ic50_cvlearner <- readRDS(ic50_file)
      
      #prepare each cvlearner obj using summary.myCV.SuperLearner
      opts = list(1, 2, 3) # need an arbitrary list of length > 1
      ic50_cv_summary <- summary.myCV.SuperLearner(object = ic50_cvlearner, opts = opts)
      
      #get the cvr2, cil, ciu, n and add to df
      ic50_est_and_ci <- get_est_and_ci(idx = 1, fit_list = list(ic50_cv_summary), Rsquared = TRUE)
      ic50_r2 <- ic50_est_and_ci$est
      ic50_cil <- ic50_est_and_ci$ci[1]
      ic50_ciu <- ic50_est_and_ci$ci[2]
      n_ic50 <- length(ic50_cvlearner$Y)
      ic50_row <- list("ic50", bnab, ic50_r2, ic50_cil, ic50_ciu, n_ic50)
      
      #add the new rows to the dataframe
      cvr2_df[nrow(cvr2_df)+1,] <- ic50_row
    }
    
    
    ic80_file <- Sys.glob(paste0("./new_docker_output/", bnab, "/cvlearner_ic80","*", ".rds"))
    if(length(ic80_file) > 0){ #if that file exists (some containers failed)
      ic80_cvlearner <- readRDS(ic80_file)
      
      #prepare each cvlearner obj using summary.myCV.SuperLearner
      opts = list(1, 2, 3) # need an arbitrary list of length > 1
      ic80_cv_summary <- summary.myCV.SuperLearner(object = ic80_cvlearner, opts = opts)
      
      #get the cvr2, cil, ciu, n and add to df
      ic80_est_and_ci <- get_est_and_ci(idx = 1, fit_list = list(ic80_cv_summary), Rsquared = TRUE)
      ic80_r2 <- ic80_est_and_ci$est
      ic80_cil <- ic80_est_and_ci$ci[1]
      ic80_ciu <- ic80_est_and_ci$ci[2]
      n_ic80 <- length(ic80_cvlearner$Y)
      ic80_row <- list("ic80", bnab, ic80_r2, ic80_cil, ic80_ciu, n_ic80)
      
      #add the new rows to the dataframe
      cvr2_df[nrow(cvr2_df)+1,] <- ic80_row
    }
  }
  return(cvr2_df)
}

#call function to get cvr2 dataframe
cvr2_df <- get_all_est_and_cis(bnabs)

#write to csv
write.csv(cvr2_df, file = "./R_output/new_cvr2_df.csv")

