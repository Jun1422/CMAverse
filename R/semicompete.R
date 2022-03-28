semicompete <- function(method = "breslow"){
####################generate time,status vector, trans metrix, covs vector, grid vectors################

time = c(NA, mediator, outcome)
status = c(NA, mediator_event, event)
trans <- transMat(x = list(c(2, 3), c(3), c()), names = c("Start", "Mediator", "Outcome"))
#generate keep vector
test = c()
test = append(test,exposure)
test = append(test,mediator)
keep = append(test,basec)
grid = seq(0,total_duration,time_grid)
ref = data.frame(astar,1,t(basecval))
colnames(ref) = c(exposure, mediator, colnames(ref)[3:length(colnames(ref))])
for (k in 2:length(names(ref))){
  if (is.numeric(data[,names(ref)[k]])) ref[k] = as.numeric(ref[k])
}
#get the dataset in right ref level
relevel(data[,exposure], ref = astar)

#method
    n <- nrow(data)
    all_term <- matrix(NA,1,7)
    est <- matrix(NA,1,3)
######################Build the model use original dataset#####################################
   #prepare the dataset
    msdata = msprep(time = time, status = status,
                    data = data, trans = trans, keep = keep)
    msdata <- expand.covs(msdata, keep, append = TRUE, longnames = FALSE)
   #build the formula
    var_name = names(msdata)
    var_length = length(var_name)
    cov_length = length(keep)
    mark_number = which(var_name ==keep[cov_length])+1

    variable = var_name[mark_number:var_length]
    variable_new = variable[which(grepl(mediator, variable) == F)]
    level_exposure = length(levels(data[,exposure]))
    
    if (is.null(formula) == T) {
      formula = paste0()
      for (m in 1:length(variable_new)) {
        if (m==1) {
          formula = paste0(formula, variable_new[m])
        } else {
          formula = paste0(formula, " + ", variable_new[m])
        }
      }
      if(level_exposure >2){
      for (x in 1:(level_exposure-1)){
        formula = paste0(formula, " + ", exposure, x, ".3*",mediator,".3")
      }
      }else{
        formula = paste0(formula, " + ", exposure, ".3*", mediator, ".3")
      }
      saveformula = paste0("Surv(Tstart, Tstop, status) ~ ", formula, " + strata(trans)")
      formula = as.formula(paste0("Surv(Tstart, Tstop, status) ~ ", formula, " + strata(trans)"))
    }

    model_fit <- coxph(formula, data=msdata, method=method, model=TRUE,control = coxph.control(timefix = FALSE))

#############################Prediction#########################################################
    #set up new data (exposure = 0)
    newd0 =  rbind(ref,ref,ref)
    newd0$trans = unique(msdata$trans)
    cg_var = names(which(sapply(data[keep], is.factor)=="TRUE"))
    for (i in 1:length(cg_var)){
      newd0[,cg_var[i]] = factor(newd0[,cg_var[i]],levels =levels(data[,cg_var[i]]))
    }
    newd0$strata=1:length(as.numeric(unique(msdata$trans)))
    attr(newd0, "trans") <- trans
    class(newd0) <- c("msdata", "data.frame")
    newd0 = expand.covs(newd0, keep, longnames = FALSE)
    w1_0<-msfit(model_fit,newd0,trans=trans)
    #set up new data (exposure = 1)
    ref_1 = ref
    ref_1[,exposure] = a
    newd1 <- rbind(ref_1,ref_1,ref_1)
    newd1$trans = unique(msdata$trans)
    for (i in 1:length(cg_var)){
      newd1[,cg_var[i]] = factor(newd1[,cg_var[i]],levels =levels(data[,cg_var[i]]))
    }
    newd1$strata=1:length(as.numeric(unique(msdata$trans)))
    attr(newd1, "trans") <- trans
    class(newd1) <- c("msdata", "data.frame")
    newd1 = expand.covs(newd1, keep, longnames = FALSE)
    w1_1<-msfit(model_fit,newd1,trans=trans)
    #get cumulative hazards for transitions 0-1 and 0-2 for race=0
    for (k in 1:nboot){

      ### Data preparation
      data_b <-data[sample(1:n, size = n, replace = TRUE),]
      msbmt <- msprep(time = time, status = status,
                      data = data_b, trans = trans, keep = keep)
      msbmt <- expand.covs(msbmt, keep, append = TRUE, longnames = FALSE)
      ### Estimation
      fit <- coxph(formula, data=msbmt, method=method, model=TRUE,,control = coxph.control(timefix = FALSE))
      ### Prediction
      ###cumulative hazards for transitions 0-1 and 0-2 for Exposure = 0
      #set up new data (exposure = 0)
      newd0 =  rbind(ref,ref,ref)
      newd0$trans = unique(msbmt$trans)
      cg_var = names(which(sapply(data_b[keep], is.factor)=="TRUE"))
      for (l in 1:length(cg_var)){
        newd0[,cg_var[l]] = factor(newd0[,cg_var[l]],levels =levels(data_b[,cg_var[l]]))
      }
      newd0$strata=1:length(as.numeric(unique(msdata$trans)))
      attr(newd0, "trans") <- trans
      class(newd0) <- c("msdata", "data.frame")
      newd0 = expand.covs(newd0, keep, longnames = FALSE)
      w1_0<-msfit(fit,newd0,trans=trans)
      #set up new data
      #a <- as.numeric(unique(msbmt$trans))
      #b <- names(which(sapply(data_b[keep], is.factor)=="TRUE"))
      #c <- rbind(ref,ref)
      #c[exposure] <- unique(data_b[exposure])
      #temp <- unique(data_b[b])
      #newd_update=merge(data.frame(temp[rep(seq_len(nrow(temp)), length(a)), ]),
                        #c, all = T)

      #attr(newd_update, "trans") <- trans
      #class(newd_update) <- c("msdata", "data.frame")
      #newd_update$strata=1:length(a)
      #newd_update$trans=newd_update$strata
      #newd_update <- expand.covs(newd_update, keep, longnames = FALSE)
      #get cumulative hazards for transitions 0-1 and 0-2 for race=0
      #get  baseline cumulative hazard(covs at reference level)
      #newd_update <- newd_update[complete.cases(newd_update),]
      #newd0 <- merge(newd_update,ref, all=F)
      #w1_0<-msfit(fit,newd0,trans=trans)

      #save the time variable
      m = sapply(1:length(grid),function(i){which.min(abs(w1_0$Haz[which(w1_0$Haz[,3]==2),1]-grid[i]))})
      #save the cumulative hazard for transitions 0-1 and 0-2 as variables
      cumhaz_01_0 <-w1_0$Haz[which(w1_0$Haz[,3]==1),2][m] #alpha(s) 0->s
      cumhaz_02_0 <- w1_0$Haz[which(w1_0$Haz[,3]==2),2][m]

      ###cumulative hazards for transitions 0-1 and 0-2 for Exposure=1
      newd1 <- rbind(ref_1,ref_1,ref_1)
      newd1$trans = unique(msbmt$trans)
      #cg_var = names(which(sapply(data_b[keep], is.factor)=="TRUE"))
      for (l in 1:length(cg_var)){
        newd1[,cg_var[l]] = factor(newd1[,cg_var[l]],levels =levels(data_b[,cg_var[l]]))
      }
      newd1$strata=1:length(as.numeric(unique(msdata$trans)))
      attr(newd1, "trans") <- trans
      class(newd1) <- c("msdata", "data.frame")
      newd1 = expand.covs(newd1, keep, longnames = FALSE)
      w1_1<-msfit(fit,newd1,trans=trans)


      #newd1 <- setdiff(newd_update,newd0)
      #w1_1<-msfit(fit,newd1,trans=trans)

      cumhaz_01_1 <-w1_1$Haz[which(w1_1$Haz[,3]==1),2][m]
      cumhaz_02_1 <-w1_1$Haz[which(w1_1$Haz[,3]==2),2][m]

      #prepare cumulative hazards for transition 1-2 matrix to save grid of time values
      #to be integrated from 0 to s, where s can take all any value in the time vector

      cumhaz_12_0 <-matrix(NA, length(grid), length(grid))
      cumhaz_12_1 <-matrix(NA, length(grid), length(grid))
      alpha01_0 <- rep(NA,length(grid))
      alpha01_1 <- rep(NA,length(grid))
      lambda12uu_1 <- rep(NA,length(grid))
      lambda12uu_0 <- rep(NA,length(grid))
      TE <- rep(NA, length(grid))
      RD <- rep(NA, length(grid))
      term <- matrix(NA, length(grid), 7)

      #the integral from 0-s can be computed for all possible time levels s
      for (i in 1:length(grid)){
        # print(i)
        newd_grid <- ref
        newd_grid[mediator] <- grid[i]
        newd_grid

        newd0 =  rbind(newd_grid,newd_grid,newd_grid)
        newd0$trans = unique(msbmt$trans)
        cg_var = names(which(sapply(data_b[keep], is.factor)=="TRUE"))
        for (j in 1:length(cg_var)){
          newd0[,cg_var[j]] = factor(newd0[,cg_var[j]],levels =levels(data_b[,cg_var[j]]))
        }
        newd0$strata=1:length(as.numeric(unique(msdata$trans)))
        attr(newd0, "trans") <- trans
        class(newd0) <- c("msdata", "data.frame")
        newd0 = expand.covs(newd0, keep, longnames = FALSE)



        w1_0_time<-msfit(fit,newd0,trans=trans)
        newd_grid_1 <- ref_1
        newd_grid_1[mediator] <- grid[i]
        newd_grid_1
        newd1 <- rbind(newd_grid_1,newd_grid_1,newd_grid_1)
        newd1$trans = unique(msbmt$trans)
        cg_var = names(which(sapply(data_b[keep], is.factor)=="TRUE"))
        for (j in 1:length(cg_var)){
          newd1[,cg_var[j]] = factor(newd1[,cg_var[j]],levels =levels(data_b[,cg_var[j]]))
        }
        newd1$strata=1:length(as.numeric(unique(msdata$trans)))
        attr(newd1, "trans") <- trans
        class(newd1) <- c("msdata", "data.frame")
        newd1 = expand.covs(newd1, keep, longnames = FALSE)

        w1_1_time<-msfit(fit,newd1,trans=trans)

        s <-grid[i]

        #save cumulative hazard up to time s for transitions 0-1 and 0-2 for race=1 and race=0
        mrow <- which.min(abs(grid-s))

        lambda01s_0 <-cumhaz_01_0[mrow][[1]]#^01 A=0
        lambda01s_1 <-cumhaz_01_1[mrow][[1]]#^01 A=1
        lambda02s_0 <-cumhaz_02_0[mrow][[1]]#^02 A=0
        lambda02s_1 <-cumhaz_02_1[mrow][[1]]#^02 A=1

        #these are the terms in the disparity and residual disparity measure that depend on the
        #cumulative hazard up to time s for transitions 0-1 and 0-2 for Exposure=1 and Exposure=0
        Prg00E <- exp(-lambda01s_0-lambda02s_1)#Pr g A=1
        Pr00e <- exp(-lambda01s_0-lambda02s_0)#Pr00 A=0
        Pr00E <-exp(-lambda01s_1-lambda02s_1)#Pr00 A=1

        #save the cumulative hazard for transitions 1-2 by time=s for all time s that we set as possible right limit of the integral and for all times at which we want to set the time to surgery in the new data
        cumhaz_12_0[,i] <-(w1_0_time$Haz[which(w1_0_time$Haz[,3]==3),2])[m]#^12 A=0
        cumhaz_12_1[,i] <-(w1_1_time$Haz[which(w1_1_time$Haz[,3]==3),2])[m]#^12 A=1
        #save the cumulative hazards up to time i when time to treatment takes time i value (elements on the diagonal of the cumulative hazard matrix for transitions 1-2)
        lambda12uu_1[i] <- cumhaz_12_1[i,i]#^12 A=0
        lambda12uu_0[i] <- cumhaz_12_0[i,i]#^12 A=1
        #estimate the hazard for transitions 0-1 and 0-2 at time=s objects
        if (i==1){
          alpha01_0[i] <- cumhaz_01_0[1]
          alpha01_1[i] <- cumhaz_01_1[1]
        }else{
          alpha01_0[i] <-cumhaz_01_0[i]-cumhaz_01_0[i-1]
          alpha01_1[i] <-cumhaz_01_1[i]-cumhaz_01_1[i-1]}

        #terms that involve integration
        Prg01E <-sum(
          exp(-cumhaz_01_0[1:mrow]
              -cumhaz_02_1[1:mrow])*
            alpha01_0[1:mrow]*
            exp(-cumhaz_12_1[mrow,1:mrow]+
                  lambda12uu_1[1:mrow]))

        Pr01e <-sum(
          exp(-cumhaz_01_0[1:mrow]
              -cumhaz_02_0[1:mrow])*
            alpha01_0[1:mrow]*
            exp(-cumhaz_12_0[mrow,1:mrow]+
                  lambda12uu_0[1:mrow]))


        Pr01E <-sum(
          exp(-cumhaz_01_1[1:mrow]
              -cumhaz_02_1[1:mrow])*
            alpha01_1[1:mrow]*
            exp(-cumhaz_12_1[mrow,1:mrow]+
                  lambda12uu_1[1:mrow]))

        #estimate the effects
        (Prg00E+Prg01E)-(Pr00e+Pr01e)
        (Pr00E+Pr01E)-(Pr00e+Pr01e)

        RD[i] <- (Prg00E+Prg01E)-(Pr00e+Pr01e)
        TE[i] <- (Pr00E+Pr01E)-(Pr00e+Pr01e)

        term[i,1] = grid[i]
        term[i,2] = Prg00E
        term[i,3] = Prg01E
        term[i,4] = Pr00e
        term[i,5] = Pr01e
        term[i,6] = Pr00E
        term[i,7] = Pr01E

      }

      all_term = rbind(all_term, term)
      est = rbind(est,cbind(TE,RD,TE-RD))
    }



    all_term = as.data.frame(all_term)
    names(all_term) = c("grid","Prg00E", "Prg01E", "Pr00e", "Pr01e", "Pr00E", "Pr01E")
    all_term = all_term[-1,]

    est = as.data.frame(est)
    names(est) = c("TE","RD","SD")
    est = est[-1,]
    est$time_grid=rep(grid, nboot)
    est = est[c("time_grid","TE","RD","SD")]
 ######################Prepare output##############################################################
    out <- list()
    method_list <- list()
    variable_list <- list()

    out$data = data

    method_list$model = 'Multistate'
    method_list$total_duration = total_duration
    method_list$cutoff = time_grid
    method_list$nboot = nboot
    out$methods = method_list

    variable_list$outcome = outcome
    variable_list$exposure = exposure
    variable_list$mediator = mediator
    variable_list$basec = basec
    out$variables =variable_list

    reg_input <- list()
    data_event = events(msdata)
    reg_input$data_event = data_event
    out$reg.input = reg_input

    ref_list <-list()
    ref_list$a = a
    ref_list$astar = astar
    ref_list$ref_table = ref
    ref_list$basecval = basecval
    out$ref = ref_list

    reg_output <- list()
    reg_output$formula = paste0("coxph(",saveformula,")")
    reg_output$model_result = model_fit
    out$reg.output = reg_output


    #out$est <- est
    #out$all_terms <- all_term
    allterms_est <- all_term %>% group_by(grid) %>% summarise(Prg00E = mean(Prg00E),
                                                                  Prg01E = mean(Prg01E),
                                                                  Pr00e = mean(Pr00e),
                                                                  Pr01e = mean(Pr01e),
                                                                  Pr00E = mean(Pr00E),
                                                                  Pr01E = mean(Pr01E))
    results <- est %>% group_by(time_grid) %>% summarise(sd_te = sd(TE),
                                                         sd_rd = sd(RD),
                                                         sd_sd = sd(SD),
                                                         TE = mean(TE),
                                                         RD = mean(RD),
                                                         SD = mean(SD),
                                                         se_te = sd_te/nboot,
                                                         se_rd = sd_rd/nboot,
                                                         se_sd = sd_sd/nboot) %>% as.data.frame()
    out$results = results

    te = quantile(sort(est$TE[est$time_grid == survival_time_fortable]),c(0.025,0.975))
    rd = quantile(sort(est$RD[est$time_grid == survival_time_fortable]),c(0.025,0.975))
    sd = quantile(sort(est$SD[est$time_grid == survival_time_fortable]),c(0.025,0.975))
    ci=data.frame(cbind(TE=paste0("(",round(te[1],4),", ",round(te[2],4),")"), RD=paste0("(",round(rd[1],4),", ",round(rd[2],4),")"), SD=paste0("(",round(sd[1],4),", ",round(sd[2],4),")")))

#the point estimator for certain time point(0-max); the user should specify the time intervals, width for grid,
    #out$time_of_interest = survival_time_fortable
    summary_table = rbind(est %>%filter(time_grid == survival_time_fortable) %>%
      group_by(time_grid) %>%
      summarise(TE = round(mean(TE),4),
                RD = round(mean(RD),4),
                SD = round(mean(SD),4)) %>% select(TE,RD,SD)
    %>% as.data.frame(), ci)
    rownames(summary_table) = c('Point estimate', '95% CI')
    out$summary_table = summary_table

    point_effect = results %>% filter(time_grid == survival_time_fortable) %>% select(TE,RD,SD)
    out$effect.pe = point_effect

    point_effect_sr = results %>% filter(time_grid == survival_time_fortable) %>% select(se_te,se_rd,se_sd)
    out$effect.se = point_effect_sr

    point_ci_low = c(te[1],rd[1],sd[1])
    names(point_ci_low) = c('TE','RD','SD')
    out$effect.ci.low = point_ci_low

    point_ci_high = c(te[2],rd[2],sd[2])
    names(point_ci_high) = c('TE','RD','SD')
    out$effect.ci.high = point_ci_high

    out$range = total_duration
    out$cutoff = time_grid
    out$time_of_interest = survival_time_fortable
    return(out)
  }
