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

ic80_10_1074 <- readRDS("/Users/sohailnizam/slapnap_redo/docker_output/cvlearner_ic80_10-1074_13Apr2021.rds")
ic50_10_1074 <- readRDS("/Users/sohailnizam/slapnap/slapnap_supplemental/ic50ic80/new_docker_output/10-1074/cvlearner_ic50_10-1074_08Apr2021.rds")

length(ic50_10_1074$Y)
length(ic80_10_1074$Y)

foo <- summary.myCV.SuperLearner(object = ic80_10_1074, opts = opts)
get_est_and_ci(idx = 1, fit_list = list(ic50_10_1074), Rsquared = TRUE)
get_est_and_ci(idx = 1, fit_list = list(foo), Rsquared = TRUE)


  
#opts list with attr 'learner', length > 1
