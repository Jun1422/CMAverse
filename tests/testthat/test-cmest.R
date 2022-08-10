context("cmest estimates causal effects correctly")

test_that("cmest works correctly for survival Y and survival M", {
  
  set.seed(1)
  # data simulation
  gen_srv <- function(n, lambda, gamma, beta, X){
    X = as.matrix(X)
    beta = as.matrix(beta, ncol=1)
    T = (-log(runif(n)) / (lambda * exp(X %*% beta)))^(1/gamma) #weibull distribution
    return(T)
  }
  n <- 10000
  c1 = sample(c(0,1),replace=TRUE, size=n,c(0.6, 0.4)) #binary confounder
  c2 = rnorm(n, mean = 1, sd = 1) #con confounder
  A = sample(c(0,1),replace=TRUE, size=n, c(0.7,0.3)) #binary exposure
  M = gen_srv(n=n, lambda = 0.1,gamma = 0.8, beta = c(-0.3,0.4,0.5), X=data.frame(A,c1,c2)) #time to event mediator
  Y = gen_srv(n=n, lambda = 0.07, gamma = 0.12, beta = c(0.2,0.3,0.4), X=data.frame(A,c1,c2)) #time to event outcome
  data = data.frame(id = c(1:n),A,c1,c2, M, Y)
  # indicator for event
  data$M_ind = ifelse(data$M <= data$Y, 1, 0)
  data$Y_ind = 1
  #data <- merge(data,full, by = "id")
  #modify Y distribution
  trans_matrix = transMat(x = list(c(2, 3), c(3), c()), names = c("A", "M", "Y"))
  covs = c("A","M", "c1","c2")
  pre_data = msprep(time = c(NA, "M", "Y"), status = c(NA, "M_ind", "Y_ind"),
                    data = data, trans = trans_matrix, keep = covs)
  pre_data = expand.covs(pre_data, covs, append = TRUE, longnames = FALSE)
  pre_data$A_M.3 = pre_data$A.3*pre_data$M.3
  # resample for T < S
  data_23= pre_data[which(pre_data$trans == 3),]
  data_23_tem = data.frame(id = rep(NA,dim(data_23)[1]),
                           new_y = rep(NA,dim(data_23)[1]))
  
  for(i in 1:dim(data_23)[1]){
    data_23_tem$id[i] = data_23$id[i]
    repeat {
      # do something
      time_test = gen_srv(n = 1, 
                          lambda = 0.1,  
                          gamma = 0.5,
                          beta = c(as.numeric(0.4),
                                   0,
                                   as.numeric(0.5),
                                   as.numeric(0.6),
                                   as.numeric(-0.2)), 
                          X = data_23[i, c("A.3", "M.3", "c1.3","c2.3", "A_M.3")])
      # exit if the condition is met
      if (time_test > data_23[i,"M.3"]) break
    }
    data_23_tem$new_y[i] = time_test
  }
  data_temp = merge(data, data_23_tem, by = "id", all = T)
  #modify Y and M
  data_temp$Y[which(data_temp$M_ind == 1)] = data_temp$new_y[which(data_temp$M_ind == 1)]
  data_temp$M[which(data_temp$M_ind == 0)] = data_temp$Y[which(data_temp$M_ind == 0)]
  data_final = data_temp
  data_final$Y_day = data_final$Y*30
  data_final$M_day = data_final$M*30
  data_final$Y_ind[which(data_final$Y > 24)] = 0 #censored data
  data_final$Y[which(data_final$Y> 24)] = 24
  data_final$Y_day[which(data_final$sY_day > 24*30)] = 24*30
  data_final$M_ind[which(data_final$M > 24)] = 0
  data_final$M[which(data_final$M > 24)] = 24
  data_final$M_day[which(data_final$M_day > 24*30)] = 24*30
  data_final$A = as.factor(data_final$A) #generate a factor exposure
  
  data = data_final %>% select(id,A,M,Y,M_ind,Y_ind,c1,c2)
  data_sub = data[sample(nrow(data), 5000), ]
  # results of cmest
  res_survsurv_multi <- cmest(data = data_sub, model = 'multistate',total_duration = 24, 
                              time_grid = 1,survival_time_fortable = 22, exposure = 'A',mediator = 'M', 
                              outcome = 'Y', event = "Y_ind",mediator_event = "M_ind", basec = c('c1','c2'),
                              basecval = c('c1' = '0','c2' = '0'),astar = '0',a='1',nboot=10)

  # ref results
  data_full = data[sample(nrow(data), 25000, replace = TRUE), ]
  
  res_survsurv_multi_ref <- cmest(data =   data_full, model = 'multistate',total_duration = 24, 
                              time_grid = 1,survival_time_fortable = 22, exposure = 'A',mediator = 'M', 
                              outcome = 'Y', event = "Y_ind",mediator_event = "M_ind", basec = c('c1','c2'),
                              basecval = c('c1' = '0','c2' = '0'),astar = '0',a='1',nboot=10)
  ref = unname(res_survsurv_multi_ref$effect.pe)
  # test
  expect_equal(unname(res_survsurv_multi$effect.pe),ref, tolerance = 0.1)
  
  
})


test_that("cmest works correctly for binary Y and continuous M", {
  
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 100
  C1 <- rnorm(n, mean = 1, sd = 1)
  C2 <- rbinom(n, 1, 0.6)
  pa <- expit(0.2 + 0.5*C1 + 0.1*C2)
  A <- rbinom(n, 1, pa)
  M <- rnorm(n, 1 + 2*A + 1.5*C1 + 0.8*C2, 1)
  py <- expit(-1 + 0.8*A + 0.5*M + 0.5*A*M + 0.3*C1 - 0.6*C2)
  Y <- rbinom(n, 1, py)
  data <- data.frame(A, M, Y, C1, C2)
  
  # results of cmest
  res_bincont <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                       mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                       mreg = list("linear"), yreg = "logistic",
                       astar = 0, a = 1, mval = list(1),
                       estimation = "paramfunc", inference = "delta")
  
  # reference results
  yreg <- glm(Y ~ A*M + C1 + C2, family = binomial, data = data)
  mreg <- lm(M ~ A + C1 + C2, data = data)
  thetas <- unname(coef(yreg))
  betas <- unname(coef(mreg))
  beta0 <- betas[1]
  beta1 <- betas[2]
  beta2 <- betas[3]
  beta3 <- betas[4]
  theta0 <- thetas[1]
  theta1 <- thetas[2]
  theta2 <- thetas[3]
  theta3 <- thetas[6]
  theta4 <- thetas[4]
  theta5 <- thetas[5]
  m <- 1
  a <- 1
  astar <- 0
  meanc1 <- mean(C1)
  meanc2 <- mean(C2)
  s <- sigma(mreg)
  cde_bincont <- exp((theta1+theta3*m)*(a-astar))
  pnde_bincont <- exp((theta1+theta3*(beta0+beta1*astar+beta2*meanc1+beta3*meanc2+theta2*s^2))*(a-astar)+0.5*theta3^2*s^2*(a-astar)^2)
  tnde_bincont <- exp((theta1+theta3*(beta0+beta1*a+beta2*meanc1+beta3*meanc2+theta2*s^2))*(a-astar)+0.5*theta3^2*s^2*(a-astar)^2)
  pnie_bincont <- exp((theta2*beta1+theta3*beta1*astar)*(a-astar))
  tnie_bincont <- exp((theta2*beta1+theta3*beta1*a)*(a-astar))
  te_bincont <- pnde_bincont*tnie_bincont
  cde_err_bincont <- (exp(theta1*(a-astar)+theta3*a*m)-exp(theta3*astar*m))*
    exp(theta2*m-(theta2+theta3*astar)*(beta0+beta1*astar+beta2*meanc1+beta3*meanc2)-0.5*(theta2+theta3*astar)^2*s^2)
  intref_err_bincont <- pnde_bincont - 1 - cde_err_bincont
  intmed_err_bincont <- tnie_bincont * pnde_bincont - pnde_bincont - pnie_bincont + 1
  pnie_err_bincont <- pnie_bincont - 1
  total_err_bincont <- te_bincont - 1
  cde_err_prop_bincont <- cde_err_bincont/total_err_bincont
  intmed_err_prop_bincont <- intmed_err_bincont/total_err_bincont
  intref_err_prop_bincont <- intref_err_bincont/total_err_bincont
  pnie_err_prop_bincont <- pnie_err_bincont/total_err_bincont
  pm_bincont <- (pnie_err_bincont+intmed_err_bincont)/total_err_bincont
  int_bincont <- (intref_err_bincont+intmed_err_bincont)/total_err_bincont
  pe_bincont <- (intref_err_bincont+intmed_err_bincont+pnie_err_bincont)/total_err_bincont
  ref <- c(cde_bincont, pnde_bincont, tnde_bincont, pnie_bincont, tnie_bincont, te_bincont, cde_err_bincont,
           intref_err_bincont, intmed_err_bincont, pnie_err_bincont, cde_err_prop_bincont,
           intref_err_prop_bincont, intmed_err_prop_bincont, pnie_err_prop_bincont,
           pm_bincont, int_bincont, pe_bincont)
  # test
  expect_equal(unname(res_bincont$effect.pe), ref)
  
})

test_that("cmest works correctly for continuous Y and continuous M", {
  
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 10000
  C1 <- rnorm(n, mean = 1, sd = 1)
  C2 <- rbinom(n, 1, 0.6)
  pa <- expit(0.2 + 0.5*C1 + 0.1*C2)
  A <- rbinom(n, 1, pa)
  M <- rnorm(n, 1 + 2*A + 1.5*C1 + 0.8*C2, 1)
  Y <- rnorm(n, -1 + 0.8*A + 0.5*M + 0.5*A*M + 0.3*C1 - 0.6*C2, 1)
  data <- data.frame(A, M, Y, C1, C2)
  
  # results of cmest
  res_contcont <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                        mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                        mreg = list("linear"), yreg = "linear",
                        astar = 0, a = 1, mval = list(1),
                        estimation = "paramfunc", inference = "delta")
  
  # reference results
  yreg <- lm(Y ~ A*M + C1 + C2, data = data)
  mreg <- lm(M ~ A + C1 + C2, data = data)
  thetas <- unname(coef(yreg))
  betas <- unname(coef(mreg))
  beta0 <- betas[1]
  beta1 <- betas[2]
  beta2 <- betas[3]
  beta3 <- betas[4]
  theta0 <- thetas[1]
  theta1 <- thetas[2]
  theta2 <- thetas[3]
  theta3 <- thetas[6]
  theta4 <- thetas[4]
  theta5 <- thetas[5]
  m <- 1
  a <- 1
  astar <- 0
  meanc1 <- mean(C1)
  meanc2 <- mean(C2)
  cde_contcont <- (theta1+theta3*m)*(a-astar)
  pnde_contcont <- (theta1+theta3*(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))*(a-astar)
  tnde_contcont <- (theta1+theta3*(beta0+beta1*a+beta2*meanc1+beta3*meanc2))*(a-astar)
  pnie_contcont <- (theta2*beta1+theta3*beta1*astar)*(a-astar)
  tnie_contcont <- (theta2*beta1+theta3*beta1*a)*(a-astar)
  te_contcont <- pnde_contcont+tnie_contcont
  intref_contcont <- pnde_contcont-cde_contcont
  intmed_contcont <- tnie_contcont-pnie_contcont
  cde_prop_contcont <- cde_contcont/te_contcont
  intref_prop_contcont <- intref_contcont/te_contcont
  intmed_prop_contcont <- intmed_contcont/te_contcont
  pnie_prop_contcont <- pnie_contcont/te_contcont
  pm_contcont <- (pnie_contcont+intmed_contcont)/te_contcont
  int_contcont <- (intref_contcont+intmed_contcont)/te_contcont
  pe_contcont <- (intref_contcont+intmed_contcont+pnie_contcont)/te_contcont
  ref <- c(cde_contcont, pnde_contcont, tnde_contcont, pnie_contcont, tnie_contcont, te_contcont,
           intref_contcont, intmed_contcont, cde_prop_contcont, intref_prop_contcont, intmed_prop_contcont,
           pnie_prop_contcont, pm_contcont, int_contcont, pe_contcont)
  
  # test
  expect_equal(unname(res_contcont$effect.pe), ref)
})

test_that("cmest works correctly for continuous Y and binary M", {
  
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 10000
  C1 <- rnorm(n, mean = 1, sd = 1)
  C2 <- rbinom(n, 1, 0.6)
  pa <- expit(0.2 + 0.5*C1 + 0.1*C2)
  A <- rbinom(n, 1, pa)
  pm <- expit(1 + 2*A + 1.5*C1 + 0.8*C2)
  M <- rbinom(n, 1, pm)
  Y <- rnorm(n, -1 + 0.8*A + 0.5*M + 0.5*A*M + 0.3*C1 - 0.6*C2, 1)
  data <- data.frame(A, M, Y, C1, C2)
  yreg <- lm(Y ~ A*M + C1 + C2, data = data)
  mreg <- glm(M ~ A + C1 + C2, family = binomial, data = data)
  
  # results of cmest
  res_contbin_rb_param_delta <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                                      mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                      mreg = list("logistic"), yreg = "linear",
                                      astar = 0, a = 1, mval = list(1),
                                      estimation = "paramfunc", inference = "delta")
  res_contbin_rb_param_bootstrap <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                                          mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                          mreg = list("logistic"), yreg = "linear",
                                          astar = 0, a = 1, mval = list(1),
                                          estimation = "paramfunc", inference = "bootstrap")
  res_contbin_rb_impu_bootstrap <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                         mreg = list("logistic"), yreg = "linear",
                                         astar = 0, a = 1, mval = list(1),
                                         estimation = "imputation", inference = "bootstrap")
  res_contbin_wb <- cmest(data = data, model = "wb", outcome = "Y", exposure = "A",
                          mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                          ereg = "logistic", yreg = "linear",
                          astar = 0, a = 1, mval = list(1),
                          estimation = "imputation", inference = "bootstrap")
  res_contbin_iorw <- cmest(data = data, model = "iorw", outcome = "Y", exposure = "A",
                            mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                            ereg = "logistic", yreg = "linear",
                            astar = 0, a = 1, mval = list(1),
                            estimation = "imputation", inference = "bootstrap")
  res_contbin_msm <- cmest(data = data, model = "msm", outcome = "Y", exposure = "A",
                           mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                           ereg = "logistic", yreg = "linear", mreg = list("logistic"),
                           wmnomreg = list("logistic"), wmdenomreg = list("logistic"),
                           astar = 0, a = 1, mval = list(1),
                           estimation = "imputation", inference = "bootstrap")
  res_contbin_ne <- cmest(data = data, model = "ne", outcome = "Y", exposure = "A",
                          mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                          yreg = "linear",
                          astar = 0, a = 1, mval = list(1),
                          estimation = "imputation", inference = "bootstrap")
  res_contbin_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                                mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                mreg = list("logistic"), yreg = "linear",
                                astar = 0, a = 1, mval = list(1),
                                estimation = "imputation", inference = "bootstrap")
  
  # reference results
  thetas <- unname(coef(yreg))
  betas <- unname(coef(mreg))
  beta0 <- betas[1]
  beta1 <- betas[2]
  beta2 <- betas[3]
  beta3 <- betas[4]
  theta0 <- thetas[1]
  theta1 <- thetas[2]
  theta2 <- thetas[3]
  theta3 <- thetas[6]
  theta4 <- thetas[4]
  theta5 <- thetas[5]
  m <- 1
  a <- 1
  astar <- 0
  meanc1 <- mean(C1)
  meanc2 <- mean(C2)
  cde_contbin <- (theta1+theta3*m)*(a-astar)
  pnde_contbin <- theta1*(a-astar) + theta3*(a-astar)*expit(beta0+beta1*astar+beta2*meanc1+beta3*meanc2)
  tnde_contbin <- theta1*(a-astar) + theta3*(a-astar)*expit(beta0+beta1*a+beta2*meanc1+beta3*meanc2)
  pnie_contbin <- (theta2+theta3*astar)*(expit(beta0+beta1*a+beta2*meanc1+beta3*meanc2)-expit(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))
  tnie_contbin <- (theta2+theta3*a)*(expit(beta0+beta1*a+beta2*meanc1+beta3*meanc2)-expit(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))
  te_contbin <- pnde_contbin+tnie_contbin
  intref_contbin <- pnde_contbin-cde_contbin
  intmed_contbin <- tnie_contbin-pnie_contbin
  cde_prop_contbin <- cde_contbin/te_contbin
  intref_prop_contbin <- intref_contbin/te_contbin
  intmed_prop_contbin <- intmed_contbin/te_contbin
  pnie_prop_contbin <- pnie_contbin/te_contbin
  pm_contbin<- (pnie_contbin+intmed_contbin)/te_contbin
  int_contbin <- (intref_contbin+intmed_contbin)/te_contbin
  pe_contbin <- (intref_contbin+intmed_contbin+pnie_contbin)/te_contbin
  ref <- c(cde_contbin, pnde_contbin, tnde_contbin, pnie_contbin, tnie_contbin, te_contbin,
           intref_contbin, intmed_contbin, cde_prop_contbin, intref_prop_contbin, intmed_prop_contbin,
           pnie_prop_contbin, pm_contbin, int_contbin, pe_contbin)
  
  # test
  expect_equal(unname(res_contbin_rb_param_delta$effect.pe), ref)
  expect_equal(unname(res_contbin_rb_param_bootstrap$effect.pe), ref)
  expect_equal(unname(res_contbin_rb_impu_bootstrap$effect.pe), ref, tolerance = 0.1)
  expect_equal(unname(res_contbin_wb$effect.pe), ref, tolerance = 0.1)
  expect_equal(unname(res_contbin_msm$effect.pe), ref, tolerance = 0.1)
  expect_equal(unname(res_contbin_ne$effect.pe), ref, tolerance = 0.1)
  expect_equal(unname(res_contbin_iorw$effect.pe), ref[c(6,2,5,15)], tolerance = 0.1)
  expect_equal(unname(res_contbin_gformula$effect.pe), ref, tolerance = 0.1)
  
})

test_that("cmest works correctly for multiple imputations", {
  
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 10000
  C1 <- rnorm(n, mean = 1, sd = 1)
  C2 <- rbinom(n, 1, 0.6)
  pa <- expit(0.2 + 0.5*C1 + 0.1*C2)
  A <- rbinom(n, 1, pa)
  pm <- expit(1 + 2*A + 1.5*C1 + 0.8*C2)
  M <- rbinom(n, 1, pm)
  Y <- rnorm(n, -1 + 0.8*A + 0.5*M + 0.5*A*M + 0.3*C1 - 0.6*C2, 1)
  data <- data.frame(A, M, Y, C1, C2)
  yreg <- lm(Y ~ A*M + C1 + C2, data = data)
  mreg <- glm(M ~ A + C1 + C2, family = binomial, data = data)
  
  # results of cmest
  res_contbin_rb_param_delta <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                                      mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                      mreg = list("logistic"), yreg = "linear",
                                      astar = 0, a = 1, mval = list(1),
                                      estimation = "paramfunc", inference = "delta", multimp = TRUE)
  res_contbin_rb_param_bootstrap <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                                          mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                          mreg = list("logistic"), yreg = "linear",
                                          astar = 0, a = 1, mval = list(1),
                                          estimation = "paramfunc", inference = "bootstrap", multimp = TRUE)
  res_contbin_rb_impu_bootstrap <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                         mreg = list("logistic"), yreg = "linear",
                                         astar = 0, a = 1, mval = list(1),
                                         estimation = "imputation", inference = "bootstrap", multimp = TRUE)
  res_contbin_wb <- cmest(data = data, model = "wb", outcome = "Y", exposure = "A",
                          mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                          ereg = "logistic", yreg = "linear",
                          astar = 0, a = 1, mval = list(1),
                          estimation = "imputation", inference = "bootstrap", multimp = TRUE)
  res_contbin_iorw <- cmest(data = data, model = "iorw", outcome = "Y", exposure = "A",
                            mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                            ereg = "logistic", yreg = "linear",
                            astar = 0, a = 1, mval = list(1),
                            estimation = "imputation", inference = "bootstrap", multimp = TRUE)
  res_contbin_msm <- cmest(data = data, model = "msm", outcome = "Y", exposure = "A",
                           mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                           ereg = "logistic", yreg = "linear", mreg = list("logistic"),
                           wmnomreg = list("logistic"), wmdenomreg = list("logistic"),
                           astar = 0, a = 1, mval = list(1),
                           estimation = "imputation", inference = "bootstrap", multimp = TRUE)
  res_contbin_ne <- cmest(data = data, model = "ne", outcome = "Y", exposure = "A",
                          mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                          yreg = "linear",
                          astar = 0, a = 1, mval = list(1),
                          estimation = "imputation", inference = "bootstrap", multimp = TRUE)
  res_contbin_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                                mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                mreg = list("logistic"), yreg = "linear",
                                astar = 0, a = 1, mval = list(1),
                                estimation = "imputation", inference = "bootstrap", multimp = TRUE)
  
  # reference results
  thetas <- unname(coef(yreg))
  betas <- unname(coef(mreg))
  beta0 <- betas[1]
  beta1 <- betas[2]
  beta2 <- betas[3]
  beta3 <- betas[4]
  theta0 <- thetas[1]
  theta1 <- thetas[2]
  theta2 <- thetas[3]
  theta3 <- thetas[6]
  theta4 <- thetas[4]
  theta5 <- thetas[5]
  m <- 1
  a <- 1
  astar <- 0
  meanc1 <- mean(C1)
  meanc2 <- mean(C2)
  cde_contbin <- (theta1+theta3*m)*(a-astar)
  pnde_contbin <- theta1*(a-astar) + theta3*(a-astar)*expit(beta0+beta1*astar+beta2*meanc1+beta3*meanc2)
  tnde_contbin <- theta1*(a-astar) + theta3*(a-astar)*expit(beta0+beta1*a+beta2*meanc1+beta3*meanc2)
  pnie_contbin <- (theta2+theta3*astar)*(expit(beta0+beta1*a+beta2*meanc1+beta3*meanc2)-expit(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))
  tnie_contbin <- (theta2+theta3*a)*(expit(beta0+beta1*a+beta2*meanc1+beta3*meanc2)-expit(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))
  te_contbin <- pnde_contbin+tnie_contbin
  intref_contbin <- pnde_contbin-cde_contbin
  intmed_contbin <- tnie_contbin-pnie_contbin
  cde_prop_contbin <- cde_contbin/te_contbin
  intref_prop_contbin <- intref_contbin/te_contbin
  intmed_prop_contbin <- intmed_contbin/te_contbin
  pnie_prop_contbin <- pnie_contbin/te_contbin
  pm_contbin<- (pnie_contbin+intmed_contbin)/te_contbin
  int_contbin <- (intref_contbin+intmed_contbin)/te_contbin
  pe_contbin <- (intref_contbin+intmed_contbin+pnie_contbin)/te_contbin
  ref <- c(cde_contbin, pnde_contbin, tnde_contbin, pnie_contbin, tnie_contbin, te_contbin,
           intref_contbin, intmed_contbin, cde_prop_contbin, intref_prop_contbin, intmed_prop_contbin,
           pnie_prop_contbin, pm_contbin, int_contbin, pe_contbin)
  
  # test
  expect_equal(unname(res_contbin_rb_param_delta$effect.pe), ref)
  expect_equal(unname(res_contbin_rb_param_bootstrap$effect.pe), ref)
  expect_equal(unname(res_contbin_rb_impu_bootstrap$effect.pe), ref, tolerance = 0.1)
  expect_equal(unname(res_contbin_wb$effect.pe), ref, tolerance = 0.1)
  expect_equal(unname(res_contbin_msm$effect.pe), ref, tolerance = 0.1)
  expect_equal(unname(res_contbin_ne$effect.pe), ref, tolerance = 0.1)
  expect_equal(unname(res_contbin_iorw$effect.pe), ref[c(6,2,5,15)], tolerance = 0.1)
  expect_equal(unname(res_contbin_gformula$effect.pe), ref, tolerance = 0.1)
  
})


test_that("cmest works correctly for bca", {
  
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 10000
  C1 <- rnorm(n, mean = 1, sd = 1)
  C2 <- rbinom(n, 1, 0.6)
  pa <- expit(0.2 + 0.5*C1 + 0.1*C2)
  A <- rbinom(n, 1, pa)
  pm <- expit(1 + 2*A + 1.5*C1 + 0.8*C2)
  M <- rbinom(n, 1, pm)
  Y <- rnorm(n, -1 + 0.8*A + 0.5*M + 0.5*A*M + 0.3*C1 - 0.6*C2, 1)
  data <- data.frame(A, M, Y, C1, C2)
  
  # results of cmest
  res_contbin_rb_param_bootstrap <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                                          mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                          mreg = list("logistic"), yreg = "linear",
                                          astar = 0, a = 1, mval = list(1),
                                          estimation = "paramfunc", inference = "bootstrap",
                                          boot.ci.type = "bca")
  res_contbin_rb_impu_bootstrap <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                         mreg = list("logistic"), yreg = "linear",
                                         astar = 0, a = 1, mval = list(1),
                                         estimation = "imputation", inference = "bootstrap",
                                         boot.ci.type = "bca")
  res_contbin_wb <- cmest(data = data, model = "wb", outcome = "Y", exposure = "A",
                          mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                          ereg = "logistic", yreg = "linear",
                          astar = 0, a = 1, mval = list(1),
                          estimation = "imputation", inference = "bootstrap",
                          boot.ci.type = "bca")
  res_contbin_iorw <- cmest(data = data, model = "iorw", outcome = "Y", exposure = "A",
                            mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                            ereg = "logistic", yreg = "linear",
                            astar = 0, a = 1, mval = list(1),
                            estimation = "imputation", inference = "bootstrap",
                            boot.ci.type = "bca")
  res_contbin_msm <- cmest(data = data, model = "msm", outcome = "Y", exposure = "A",
                           mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                           ereg = "logistic", yreg = "linear", mreg = list("logistic"),
                           wmnomreg = list("logistic"), wmdenomreg = list("logistic"),
                           astar = 0, a = 1, mval = list(1),
                           estimation = "imputation", inference = "bootstrap",
                           boot.ci.type = "bca")
  res_contbin_ne <- cmest(data = data, model = "ne", outcome = "Y", exposure = "A",
                          mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                          yreg = "linear",
                          astar = 0, a = 1, mval = list(1),
                          estimation = "imputation", inference = "bootstrap",
                          boot.ci.type = "bca")
  res_contbin_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                                mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                mreg = list("logistic"), yreg = "linear",
                                astar = 0, a = 1, mval = list(1),
                                estimation = "imputation", inference = "bootstrap",
                                boot.ci.type = "bca")
  
  # test
  ref <- summary(res_contbin_rb_param_bootstrap)$summarydf$Estimate
  expect_equal(summary(res_contbin_rb_impu_bootstrap)$summarydf$Estimate, ref, tolerance = 0.1)
  expect_equal(summary(res_contbin_wb)$summarydf$Estimate, ref, tolerance = 0.1)
  expect_equal(summary(res_contbin_msm)$summarydf$Estimate, ref, tolerance = 0.1)
  expect_equal(summary(res_contbin_ne)$summarydf$Estimate, ref, tolerance = 0.1)
  expect_equal(summary(res_contbin_iorw)$summarydf$Estimate, ref[c(6,2,5,15)], tolerance = 0.1)
  expect_equal(summary(res_contbin_gformula)$summarydf$Estimate, ref, tolerance = 0.1)
  expect_equal(class(print(res_contbin_rb_param_bootstrap)), "list")
  expect_equal(class(print(res_contbin_rb_impu_bootstrap)), "list")
  expect_equal(class(print(res_contbin_wb)), "list")
  expect_equal(class(print(res_contbin_msm)), "list")
  expect_equal(class(print(res_contbin_ne)), "list")
  expect_equal(class(print(res_contbin_iorw)), "list")
  expect_equal(class(print(res_contbin_gformula)), "list")
  expect_equal(class(print(summary(res_contbin_rb_impu_bootstrap))), "list")
  expect_equal(class(print(summary(res_contbin_wb))), "list")
  expect_equal(class(print(summary(res_contbin_msm))), "list")
  expect_equal(class(print(summary(res_contbin_ne))), "list")
  expect_equal(class(print(summary(res_contbin_iorw))), "list")
  expect_equal(class(print(summary(res_contbin_gformula))), "list")
  
})


test_that("cmest works correctly for survival Y and count M", {
  
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 10000
  C1 <- rnorm(n, mean = 1, sd = 1)
  C2 <- rbinom(n, 1, 0.6)
  pa <- expit(0.2 + 0.5*C1 - 0.3*C2)
  A <- rbinom(n, 1, pa)
  M <- rpois(n, exp(1 + 0.2*A + 0.2*C1 + 0.1*C2))
  Y <- rexp(n,exp(-1 - 0.3*A - 0.1*M - 0.1*A*M + 0.3*C1 - 0.6*C2))
  cen <- quantile(Y, probs = 0.01)
  delta <- as.numeric(Y > cen)
  data <- data.frame(A, M, Y, C1, C2, delta)
  ereg <- glm(A ~ C1 + C2, family = binomial(), data = data)
  yreg <- survival::survreg(survival::Surv(Y, delta) ~ A*M + C1 + C2, data = data, dist = "weibull")
  mreg <- glm(M ~ A + C1 + C2, family = poisson(), data = data)
  
  # results of cmest
  res_survcount_rb <- cmest(data = data, model = "rb", outcome = "Y", event = "delta",
                            exposure = "A", mediator = "M", basec = c("C1", "C2"), 
                            EMint = TRUE,
                            mreg = list("poisson"), yreg = "aft_weibull",
                            astar = 0, a = 1, mval = list(1),
                            estimation = "imputation", inference = "bootstrap")
  res_survcount_iorw <- cmest(data = data, model = "iorw", outcome = "Y", exposure = "A",
                              mediator = "M", basec = c("C1", "C2"), EMint = TRUE, event = "delta",
                              ereg = "logistic", yreg = "aft_weibull",
                              astar = 0, a = 1, 
                              estimation = "imputation", inference = "bootstrap")
  res_survcount_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                                  mediator = "M", basec = c("C1", "C2"), EMint = TRUE, event = "delta",
                                  mreg = list("poisson"), yreg = "aft_weibull",
                                  astar = 0, a = 1, mval = list(1),
                                  estimation = "imputation", inference = "bootstrap")
  
  # test
  expect_equal(unname(res_survcount_rb$effect.pe)[c(6,2,5,15)], 
               unname(res_survcount_iorw$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_survcount_rb$effect.pe), 
               unname(res_survcount_gformula$effect.pe), tolerance = 0.1)
  
  # results of cmest
  res_survcount_rb <- cmest(data = data, model = "rb", outcome = "Y", event = "delta",
                            exposure = "A", mediator = "M", basec = c("C1", "C2"), 
                            EMint = TRUE,
                            mreg = list("negbin"), yreg = "coxph",
                            astar = 0, a = 1, mval = list(1),
                            estimation = "imputation", inference = "bootstrap")
  res_survcount_iorw <- cmest(data = data, model = "iorw", outcome = "Y", exposure = "A",
                              mediator = "M", basec = c("C1", "C2"), EMint = TRUE, event = "delta",
                              ereg = "logistic", yreg = "coxph",
                              astar = 0, a = 1, 
                              estimation = "imputation", inference = "bootstrap")
  res_survcount_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                                  mediator = "M", basec = c("C1", "C2"), EMint = TRUE, event = "delta",
                                  mreg = list("negbin"), yreg = "coxph",
                                  astar = 0, a = 1, mval = list(1),
                                  estimation = "imputation", inference = "bootstrap")
  
  # test
  expect_equal(unname(res_survcount_rb$effect.pe)[c(6,2,5,15)], 
               unname(res_survcount_iorw$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_survcount_rb$effect.pe), 
               unname(res_survcount_gformula$effect.pe), tolerance = 0.1)
  
})


test_that("cmest works correctly for count Y and count M", {
  
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 10000
  C1 <- rnorm(n, mean = 1, sd = 1)
  C2 <- rbinom(n, 1, 0.6)
  pa <- expit(0.2 + 0.5*C1 - 0.3*C2)
  A <- rbinom(n, 1, pa)
  M <- rpois(n, exp(1 + 0.2*A + 0.2*C1 + 0.1*C2))
  Y <- rpois(n,exp(-1 - 0.3*A - 0.1*M - 0.1*A*M + 0.3*C1 - 0.6*C2))
  data <- data.frame(A, M, Y, C1, C2)
  ereg <- glm(A ~ C1 + C2, family = binomial(), data = data)
  yreg <- glm(Y ~ A*M + C1 + C2, family = quasipoisson(), data = data)
  mreg <- glm(M ~ A + C1 + C2, family = quasipoisson(), data = data)
  
  # results of cmest
  res_countcount_rb <- cmest(data = data, model = "rb", outcome = "Y", event = "delta",
                             exposure = "A", mediator = "M", basec = c("C1", "C2"), 
                             EMint = TRUE,
                             mreg = list("quasipoisson"), yreg = "quasipoisson",
                             astar = 0, a = 1, mval = list(1),
                             estimation = "imputation", inference = "bootstrap")
  res_countcount_wb <- cmest(data = data, model = "wb", outcome = "Y", event = "delta", exposure = "A",
                             mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                             ereg = "logistic", yreg = "quasipoisson",
                             astar = 0, a = 1, mval = list(1),
                             estimation = "imputation", inference = "bootstrap")
  res_countcount_iorw <- cmest(data = data, model = "iorw", outcome = "Y", exposure = "A",
                               mediator = "M", basec = c("C1", "C2"), EMint = TRUE, event = "delta",
                               ereg = "logistic", yreg = "quasipoisson",
                               astar = 0, a = 1, 
                               estimation = "imputation", inference = "bootstrap")
  res_countcount_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                                   mediator = "M", basec = c("C1", "C2"), EMint = TRUE, event = "delta",
                                   mreg = list("quasipoisson"), yreg = "quasipoisson",
                                   astar = 0, a = 1, mval = list(1),
                                   estimation = "imputation", inference = "bootstrap")
  
  # test
  expect_equal(unname(res_countcount_rb$effect.pe), 
               unname(res_countcount_wb$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_countcount_rb$effect.pe)[c(6,2,5,15)], 
               unname(res_countcount_iorw$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_countcount_rb$effect.pe), 
               unname(res_countcount_gformula$effect.pe), tolerance = 0.1)
  
})


test_that("cmest works correctly for survival Y and ordinal M", {
  
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 10000
  C1 <- rnorm(n, mean = 1, sd = 1)
  C2 <- rbinom(n, 1, 0.6)
  proba1to2 <- expit(0.4 + 0.5*C1 - 0.3*C2)
  proba2 <- expit(-0.8 + 0.5*C1 - 0.3*C2)
  proba0 <- 1 - proba1to2
  proba1 <- proba1to2 - proba2
  A <- sapply(1:n, FUN = function(x) sample(c(0,1,2),size=1,replace=TRUE,
                                            prob=c(proba0[x],
                                                   proba1[x],
                                                   proba2[x])))
  A <- as.factor(A)
  probm1to2 <- expit(0.5 + 0.2*(A == 1) + 0.1*(A == 2) + 0.2*C1 + 0.1*C2)
  probm2 <- expit(-1 + 0.2*(A == 1) + 0.1*(A == 2) + 0.2*C1 + 0.1*C2)
  probm0 <- 1 - probm1to2
  probm1 <- probm1to2 - probm2
  M <- sapply(1:n, FUN = function(x) sample(c(0,1,2),size=1,replace=TRUE,
                                            prob=c(probm0[x],
                                                   probm1[x],
                                                   probm2[x])))
  M <- as.factor(M)
  Y <- rexp(n,exp(-1 - 0.3*(A == 1) - 0.1*(A == 2) - 0.1*(M == 1) - 0.3*(M == 2) -
                    0.1*(A == 1)*(M == 1) - 0.3*(A == 1)*(M == 2) -
                    0.2*(A == 2)*(M == 1) - 0.1*(A == 2)*(M == 2) +
                    0.3*C1 - 0.6*C2))
  cen <- quantile(Y, probs = 0.01)
  delta <- as.numeric(Y > cen)
  data <- data.frame(A, M, Y, C1, C2, delta)
  ereg <- MASS::polr(A ~ C1 + C2, data = data)
  yreg <- survival::survreg(survival::Surv(Y, delta) ~ A*M + C1 + C2, data = data, dist = "exponential")
  mreg <- MASS::polr(M ~ A + C1 + C2, data = data)
  
  # results of cmest
  res_survordinal_rb <- cmest(data = data, model = "rb", outcome = "Y", event = "delta",
                              exposure = "A", mediator = "M", basec = c("C1", "C2"),
                              EMint = TRUE,
                              mreg = list("ordinal"), yreg = "aft_exp",
                              astar = 0, a = 1, mval = list(1),
                              estimation = "imputation", inference = "bootstrap")
  res_survordinal_iorw <- cmest(data = data, model = "iorw", outcome = "Y", exposure = "A",
                                mediator = "M", basec = c("C1", "C2"), EMint = TRUE, event = "delta",
                                ereg = "ordinal", yreg = "aft_exp",
                                astar = 0, a = 1,
                                estimation = "imputation", inference = "bootstrap")
  res_survordinal_msm <- cmest(data = data, model = "msm",
                               outcome = "Y", exposure = "A",
                               mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                               ereg = "ordinal", yreg = "aft_exp", mreg = list("ordinal"),
                               wmnomreg = list("ordinal"), wmdenomreg = list("ordinal"),
                               astar = 0, a = 1, mval = list(1),
                               estimation = "imputation", inference = "bootstrap")
  res_survordinal_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                                    mediator = "M", basec = c("C1", "C2"), EMint = TRUE, event = "delta",
                                    mreg = list("ordinal"), yreg = "aft_exp",
                                    astar = 0, a = 1, mval = list(1),
                                    estimation = "imputation", inference = "bootstrap")
  
  # test
  expect_equal(unname(res_survordinal_rb$effect.pe)[c(6,2,5,15)],
               unname(res_survordinal_iorw$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_survordinal_rb$effect.pe),
               unname(res_survordinal_msm$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_survordinal_rb$effect.pe),
               unname(res_survordinal_gformula$effect.pe), tolerance = 0.1)
  
})


test_that("cmest works correctly for ordinal Y and ordinal M", {
  
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 10000
  C1 <- rnorm(n, mean = 1, sd = 1)
  C2 <- rbinom(n, 1, 0.6)
  proba1to2 <- expit(0.4 + 0.5*C1 - 0.3*C2)
  proba2 <- expit(-0.8 + 0.5*C1 - 0.3*C2)
  proba0 <- 1 - proba1to2
  proba1 <- proba1to2 - proba2
  A <- sapply(1:n, FUN = function(x) sample(c(0,1,2),size=1,replace=TRUE,
                                            prob=c(proba0[x],
                                                   proba1[x],
                                                   proba2[x])))
  A <- as.factor(A)
  probm1to2 <- expit(0.5 + 0.2*(A == 1) + 0.1*(A == 2) + 0.2*C1 + 0.1*C2)
  probm2 <- expit(-1 + 0.2*(A == 1) + 0.1*(A == 2) + 0.2*C1 + 0.1*C2)
  probm0 <- 1 - probm1to2
  probm1 <- probm1to2 - probm2
  M <- sapply(1:n, FUN = function(x) sample(c(0,1,2),size=1,replace=TRUE,
                                            prob=c(probm0[x],
                                                   probm1[x],
                                                   probm2[x])))
  M <- as.factor(M)
  proby1to2 <- expit(1 - 0.3*(A == 1) - 0.1*(A == 2) - 0.1*(M == 1) - 0.3*(M == 2) -
                       0.1*(A == 1)*(M == 1) - 0.3*(A == 1)*(M == 2) -
                       0.2*(A == 2)*(M == 1) - 0.1*(A == 2)*(M == 2) +
                       0.3*C1 - 0.6*C2)
  proby2 <- expit(-0.5 - 0.3*(A == 1) - 0.1*(A == 2) - 0.1*(M == 1) - 0.3*(M == 2) -
                    0.1*(A == 1)*(M == 1) - 0.3*(A == 1)*(M == 2) -
                    0.2*(A == 2)*(M == 1) - 0.1*(A == 2)*(M == 2) +
                    0.3*C1 - 0.6*C2)
  proby0 <- 1 - proby1to2
  proby1 <- proby1to2 - proby2
  Y <- sapply(1:n, FUN = function(x) sample(c(0,1,2), size = 1, replace = TRUE,
                                            prob=c(proby0[x],
                                                   proby1[x],
                                                   proby2[x])))
  Y <- as.factor(Y)
  data <- data.frame(A, M, Y, C1, C2)
  ereg <- MASS::polr(A ~ C1 + C2, data = data)
  yreg <- MASS::polr(Y ~ A*M + C1 + C2, data = data)
  mreg <- MASS::polr(M ~ A + C1 + C2, data = data)
  
  # results of cmest
  res_ordinalordinal_rb <- cmest(data = data, model = "rb", outcome = "Y", event = "delta",
                                 exposure = "A", mediator = "M", basec = c("C1", "C2"),
                                 EMint = TRUE,
                                 mreg = list("ordinal"), yreg = "ordinal", yref = "1",
                                 astar = 0, a = 1, mval = list(1),
                                 estimation = "imputation", inference = "bootstrap")
  res_ordinalordinal_iorw <- cmest(data = data, model = "iorw", outcome = "Y", exposure = "A",
                                   mediator = "M", basec = c("C1", "C2"), EMint = TRUE, event = "delta",
                                   ereg = "ordinal", yreg = "ordinal", yref = "1",
                                   astar = 0, a = 1,
                                   estimation = "imputation", inference = "bootstrap")
  res_ordinalordinal_msm <- cmest(data = data, model = "msm",
                                  outcome = "Y", exposure = "A",
                                  mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                  ereg = "ordinal", yreg = "ordinal", mreg = list("ordinal"),
                                  wmnomreg = list("ordinal"), wmdenomreg = list("ordinal"),
                                  astar = 0, a = 1, mval = list(1), yref = "1",
                                  estimation = "imputation", inference = "bootstrap")
  res_ordinalordinal_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                                       mediator = "M", basec = c("C1", "C2"), EMint = TRUE, event = "delta",
                                       mreg = list("ordinal"), yreg = "ordinal",
                                       astar = 0, a = 1, mval = list(1), yref = "1",
                                       estimation = "imputation", inference = "bootstrap")
  
  # test
  expect_equal(unname(res_ordinalordinal_rb$effect.pe)[c(6,2,5)],
               unname(res_ordinalordinal_iorw$effect.pe)[1:3], tolerance = 0.1)
  expect_equal(unname(res_ordinalordinal_rb$effect.pe),
               unname(res_ordinalordinal_msm$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_ordinalordinal_rb$effect.pe),
               unname(res_ordinalordinal_gformula$effect.pe), tolerance = 0.1)
  
})


test_that("cmest works correctly for continuous Y and gamma M", {
  
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 10000
  C1 <- rnorm(n, mean = 1, sd = 0.1)
  C2 <- rbinom(n, 1, 0.6)
  pa <- expit(0.2 + 0.5*C1 + 0.1*C2)
  A <- rbinom(n, 1, pa)
  M <- rgamma(n, shape = 2, scale = exp(-3 + 2*A + 1.5*C1 + 0.8*C2)/2)
  Y <- rnorm(n, -3 + 0.8*A - 1.8*M + 0.5*A*M + 0.3*C1 - 0.6*C2, 0.5)
  data <- data.frame(A, M, Y, C1, C2)
  mreg <- glm(M ~ A + C1 + C2, data = data, family = Gamma("log"))
  # results of cmest
  res_gammalinear_rb <- cmest(data = data, model = "rb", outcome = "Y", 
                              exposure = "A", mediator = "M", basec = c("C1", "C2"), 
                              EMint = TRUE,
                              mreg = list(mreg), yreg = "linear",
                              astar = 0, a = 1, mval = list(1),
                              estimation = "imputation", inference = "bootstrap")
  res_gammalinear_wb <- cmest(data = data, model = "wb", outcome = "Y", 
                              exposure = "A", mediator = "M", basec = c("C1", "C2"), 
                              EMint = TRUE,
                              ereg = "logistic", yreg = "linear",
                              astar = 0, a = 1, mval = list(1),
                              estimation = "imputation", inference = "bootstrap")
  res_gammalinear_iorw <- cmest(data = data, model = "iorw", outcome = "Y", exposure = "A",
                                mediator = "M", basec = c("C1", "C2"), EMint = TRUE, 
                                ereg = "logistic", yreg = "linear",
                                astar = 0, a = 1, 
                                estimation = "imputation", inference = "bootstrap")
  res_gammalinear_ne <- cmest(data = data, model = "ne", outcome = "Y", exposure = "A",
                              mediator = "M", basec = c("C1", "C2"), EMint = TRUE, 
                              yreg = "linear",
                              astar = 0, a = 1, mval = list(1),
                              estimation = "imputation", inference = "bootstrap")
  res_gammalinear_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                                    mediator = "M", basec = c("C1", "C2"), EMint = TRUE, 
                                    mreg = list(mreg), yreg = "linear",
                                    astar = 0, a = 1, mval = list(1),
                                    estimation = "imputation", inference = "bootstrap")
  # test
  expect_equal(unname(res_gammalinear_rb$effect.pe), 
               unname(res_gammalinear_wb$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_gammalinear_rb$effect.pe)[c(6,2,5,15)], 
               unname(res_gammalinear_iorw$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_gammalinear_rb$effect.pe), 
               unname(res_gammalinear_ne$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_gammalinear_rb$effect.pe), 
               unname(res_gammalinear_gformula$effect.pe), tolerance = 0.1)
  
})


test_that("cmest works correctly for multiple mediators", {
  
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 10000
  C1 <- rnorm(n, mean = 1, sd = 0.1)
  C2 <- rbinom(n, 1, 0.6)
  pa <- expit(0.2 + 0.5*C1 + 0.1*C2)
  A <- rbinom(n, 1, pa)
  pm <- expit(1 + 2*A + 1.5*C1 + 0.8*C2)
  M1 <- rbinom(n, 1, expit(1 + 2*A + 1.5*C1 + 0.8*C2))
  M2 <- rpois(n, exp(1 + 0.2*A + 0.3*M1 + 0.2*C1 + 0.1*C2))
  Y <- rnorm(n, -3 + 0.8*A - 1.8*M1 + 0.6*M2 + 0.5*A*M1 + 0.8*A*M2 + 0.3*C1 - 0.6*C2, 0.5)
  data <- data.frame(A, M1, M2, Y, C1, C2)
  
  # results of cmest
  res_multipleM_rb <- cmest(data = data, model = "rb", outcome = "Y", 
                            exposure = "A", mediator = c("M1", "M2"), basec = c("C1", "C2"), 
                            EMint = TRUE,
                            mreg = list("logistic", "poisson"), yreg = "linear",
                            astar = 0, a = 1, mval = list(1, 5),
                            estimation = "imputation", inference = "bootstrap")
  res_multipleM_wb <- cmest(data = data, model = "wb", outcome = "Y", 
                            exposure = "A", mediator = c("M1", "M2"), basec = c("C1", "C2"), 
                            EMint = TRUE,
                            ereg = "logistic", yreg = "linear",
                            astar = 0, a = 1, mval = list(1, 5),
                            estimation = "imputation", inference = "bootstrap")
  res_multipleM_iorw <- cmest(data = data, model = "iorw", outcome = "Y", exposure = "A",
                              mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE, 
                              ereg = "logistic", yreg = "linear",
                              astar = 0, a = 1, 
                              estimation = "imputation", inference = "bootstrap")
  res_multipleM_ne <- cmest(data = data, model = "ne", outcome = "Y", exposure = "A",
                            mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE, 
                            yreg = "linear",
                            astar = 0, a = 1, mval = list(1, 5),
                            estimation = "imputation", inference = "bootstrap")
  res_multipleM_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                                  mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE, 
                                  mreg = list("logistic", "poisson"), yreg = "linear",
                                  astar = 0, a = 1, mval = list(1, 5),
                                  estimation = "imputation", inference = "bootstrap")
  # test
  expect_equal(unname(res_multipleM_rb$effect.pe), 
               unname(res_multipleM_wb$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_multipleM_rb$effect.pe)[c(6,2,5,15)], 
               unname(res_multipleM_iorw$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_multipleM_rb$effect.pe), 
               unname(res_multipleM_ne$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_multipleM_rb$effect.pe), 
               unname(res_multipleM_gformula$effect.pe), tolerance = 0.1)
  
})


test_that("cmest works correctly for binary Y and binary M with postc", {
  
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 10000
  C1 <- rnorm(n, mean = 1, sd = 0.1)
  C2 <- rbinom(n, 1, 0.6)
  pa <- expit(0.2 + 0.5*C1 + 0.1*C2)
  A <- rbinom(n, 1, pa)
  L <- rnorm(n, mean = 0.5 - A + 0.2*C1 + 0.6*C2)
  pm <- expit(1 + 2*A + 1.5*C1 + 0.8*C2 + L)
  M <- rbinom(n, 1, pm)
  py <- expit(-3 + 0.8*A - 1.8*M + 0.5*A*M + 0.3*C1 - 0.6*C2 + L)
  Y <- rbinom(n, 1, py)
  data <- data.frame(A, M, Y, C1, C2, L)
  
  # results of cmest
  expect_error(cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                     mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                     mreg = list("logistic"), yreg = "logistic",
                     astar = 0, a = 1, mval = list(1),
                     estimation = "paramfunc", inference = "delta"),
               "When postc is not empty, select model from 'msm' and 'gformula'")
  expect_error(cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                     mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                     mreg = list("logistic"), yreg = "logistic",
                     astar = 0, a = 1, mval = list(1),
                     estimation = "paramfunc", inference = "delta"),
               "When postc is not empty, select model from 'msm' and 'gformula'")
  expect_error(cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                     mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                     mreg = list("logistic"), yreg = "logistic",
                     astar = 0, a = 1, mval = list(1),
                     estimation = "paramfunc", inference = "bootstrap"),
               "When postc is not empty, select model from 'msm' and 'gformula'")
  expect_error(cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                     mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                     mreg = list("logistic"), yreg = "logistic",
                     astar = 0, a = 1, mval = list(1),
                     estimation = "imputation", inference = "bootstrap"),
               "When postc is not empty, select model from 'msm' and 'gformula'")
  expect_error(cmest(data = data, model = "wb", outcome = "Y", exposure = "A",
                     mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                     ereg = "logistic", yreg = "logistic",
                     astar = 0, a = 1, mval = list(1),
                     estimation = "imputation", inference = "bootstrap"),
               "When postc is not empty, select model from 'msm' and 'gformula'")
  expect_error(cmest(data = data, model = "iorw", outcome = "Y", exposure = "A",
                     mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                     ereg = "logistic", yreg = "logistic",
                     astar = 0, a = 1, mval = list(1),
                     estimation = "imputation", inference = "bootstrap"),
               "When postc is not empty, select model from 'msm' and 'gformula'")
  res_binbin_msm <- cmest(data = data, model = "msm", outcome = "Y", exposure = "A",
                          mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                          ereg = "logistic", yreg = "logistic", mreg = list("logistic"),
                          wmnomreg = list("logistic"), wmdenomreg = list("logistic"),
                          astar = 0, a = 1, mval = list(1),
                          estimation = "imputation", inference = "bootstrap")
  expect_error(cmest(data = data, model = "ne", outcome = "Y", exposure = "A",
                     mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                     yreg = "logistic",
                     astar = 0, a = 1, mval = list(1),
                     estimation = "imputation", inference = "bootstrap"),
               "When postc is not empty, select model from 'msm' and 'gformula'")
  res_binbin_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                               mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                               mreg = list("logistic"), yreg = "logistic", postcreg = list("linear"),
                               astar = 0, a = 1, mval = list(1),
                               estimation = "imputation", inference = "bootstrap")
  
  # test
  expect_equal(summary(res_binbin_gformula)$summarydf$Estimate,
               summary(res_binbin_msm)$summarydf$Estimate, tolerance = 0.1)
  expect_equal(class(print(res_binbin_msm)), "list")
  expect_equal(class(print(res_binbin_gformula)), "list")
  expect_equal(class(print(summary(res_binbin_msm))), "list")
  expect_equal(class(print(summary(res_binbin_gformula))), "list")
  
  res_binbin_msm <- cmest(data = data, model = "msm", outcome = "Y", exposure = "A",
                          mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                          ereg = "logistic", yreg = "logistic", mreg = list("logistic"),
                          wmnomreg = list("logistic"), wmdenomreg = list("logistic"),
                          astar = 0, a = 1, mval = list(1),
                          estimation = "imputation", inference = "bootstrap", multimp = TRUE)
  res_binbin_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                               mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                               mreg = list("logistic"), yreg = "logistic", postcreg = list("linear"),
                               astar = 0, a = 1, mval = list(1),
                               estimation = "imputation", inference = "bootstrap", multimp = TRUE)
  
  # test
  expect_equal(summary(res_binbin_gformula)$summarydf$Estimate,
               summary(res_binbin_msm)$summarydf$Estimate, tolerance = 0.1)
  expect_equal(class(print(res_binbin_msm)), "list")
  expect_equal(class(print(res_binbin_gformula)), "list")
  expect_equal(class(print(summary(res_binbin_msm))), "list")
  expect_equal(class(print(summary(res_binbin_gformula))), "list")
  
})

test_that("cmest works correctly for regressions with prior weights", {
  
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 1000000
  C1 <- rnorm(n, mean = 1, sd = 0.1)
  C2 <- rbinom(n, 1, 0.6)
  pa <- expit(0.2 + 0.5*C1 + 0.1*C2)
  A <- rbinom(n, 1, pa)
  M1 <- rbinom(n, 1, expit(0.5 + 0.2*A + 0.5*C1 - 0.8*C2))
  sum(M1)
  M2 <- rbinom(n, 1, expit(-1 + A + M1 - 0.5*C1 + 0.5*C2))
  sum(M2)
  py <- expit(-6 + 0.8*A - 1.2*M1 + 0.5*M2 + 0.3*A*M1 + 0.2*A*M2 + 0.3*C1 - 0.6*C2)
  Y <- rbinom(n, 1, py)
  yprevalence <- sum(Y)/n
  data <- data.frame(A, M1, M2, Y, C1, C2)
  case_indice <- sample(which(data$Y == 1), 2000, replace = FALSE)
  control_indice <- sample(which(data$Y == 0), 2000, replace = FALSE)
  data <- data[c(case_indice, control_indice), ]
  mreg1 <- glm(M1 ~ A + C1 + C2, data = data, family = binomial, weights = rep(2, 4000))
  mreg2 <- glm(M2 ~ A + C1 + C2, data = data, family = binomial, weights = rep(2, 4000))
  mreg1_msm <- glm(M1 ~ A, data = data, family = binomial, weights = rep(2, 4000))
  mreg2_msm <- glm(M2 ~ A, data = data, family = binomial, weights = rep(2, 4000))
  
  # yrare = TRUE
  # results of cmest
  res_binbin_rb_impu_bootstrap <- cmest(data = data, model = "rb", casecontrol = TRUE, yrare = TRUE,
                                        outcome = "Y", exposure = "A",
                                        mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                                        mreg = list(mreg1, mreg2), yreg = "logistic",
                                        astar = 0, a = 1, mval = list(1, 1),
                                        estimation = "imputation", inference = "bootstrap")
  res_binbin_wb <- cmest(data = data, model = "wb", casecontrol = TRUE, yrare = TRUE,
                         outcome = "Y", exposure = "A",
                         mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                         ereg = "logistic", yreg = "logistic",
                         astar = 0, a = 1, mval = list(1, 1),
                         estimation = "imputation", inference = "bootstrap")
  res_binbin_iorw <- cmest(data = data, model = "iorw", casecontrol = TRUE, yrare = TRUE,
                           outcome = "Y", exposure = "A",
                           mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                           ereg = "logistic", yreg = "logistic",
                           astar = 0, a = 1, mval = list(1, 1),
                           estimation = "imputation", inference = "bootstrap")
  res_binbin_msm <- cmest(data = data, model = "msm", casecontrol = TRUE, yrare = TRUE,
                          outcome = "Y", exposure = "A",
                          mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                          ereg = "logistic", yreg = "logistic", mreg = list(mreg1_msm, mreg2_msm),
                          wmnomreg = list("logistic", "logistic"), wmdenomreg = list("logistic", "logistic"),
                          astar = 0, a = 1, mval = list(1, 1),
                          estimation = "imputation", inference = "bootstrap")
  res_binbin_ne <- cmest(data = data, model = "ne", casecontrol = TRUE, yrare = TRUE,
                         outcome = "Y", exposure = "A",
                         mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                         yreg = "logistic",
                         astar = 0, a = 1, mval = list(1, 1),
                         estimation = "imputation", inference = "bootstrap")
  res_binbin_gformula <- cmest(data = data, model = "gformula", casecontrol = TRUE, yrare = TRUE,
                               outcome = "Y", exposure = "A",
                               mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                               mreg = list(mreg1, mreg2), yreg = "logistic",
                               astar = 0, a = 1, mval = list(1, 1),
                               estimation = "imputation", inference = "bootstrap")
  
  # test
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe), 
               unname(res_binbin_wb$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe)[c(6,2,5,15)], 
               unname(res_binbin_iorw$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe), 
               unname(res_binbin_ne$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe), 
               unname(res_binbin_msm$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe), 
               unname(res_binbin_gformula$effect.pe), tolerance = 0.1)
  
  # yprevalence = yprevalence
  # results of cmest
  res_binbin_rb_impu_bootstrap <- cmest(data = data, model = "rb",
                                        casecontrol = TRUE, yprevalence = yprevalence,
                                        outcome = "Y", exposure = "A",
                                        mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                                        mreg = list(mreg1, mreg2), yreg = "logistic",
                                        astar = 0, a = 1, mval = list(1, 1),
                                        estimation = "imputation", inference = "bootstrap")
  res_binbin_wb <- cmest(data = data, model = "wb", casecontrol = TRUE, yprevalence = yprevalence,
                         outcome = "Y", exposure = "A",
                         mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                         ereg = "logistic", yreg = "logistic",
                         astar = 0, a = 1, mval = list(1, 1),
                         estimation = "imputation", inference = "bootstrap")
  res_binbin_iorw <- cmest(data = data, model = "iorw", casecontrol = TRUE, yprevalence = yprevalence,
                           outcome = "Y", exposure = "A",
                           mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                           ereg = "logistic", yreg = "logistic",
                           astar = 0, a = 1, mval = list(1, 1),
                           estimation = "imputation", inference = "bootstrap")
  res_binbin_msm <- cmest(data = data, model = "msm", casecontrol = TRUE, yprevalence = yprevalence,
                          outcome = "Y", exposure = "A",
                          mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                          ereg = "logistic", yreg = "logistic", mreg = list(mreg1_msm, mreg2_msm),
                          wmnomreg = list("logistic", "logistic"), wmdenomreg = list("logistic", "logistic"),
                          astar = 0, a = 1, mval = list(1, 1),
                          estimation = "imputation", inference = "bootstrap")
  res_binbin_ne <- cmest(data = data, model = "ne", casecontrol = TRUE, yprevalence = yprevalence,
                         outcome = "Y", exposure = "A",
                         mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                         yreg = "logistic",
                         astar = 0, a = 1, mval = list(1, 1),
                         estimation = "imputation", inference = "bootstrap")
  res_binbin_gformula <- cmest(data = data, model = "gformula", casecontrol = TRUE, yprevalence = yprevalence,
                               outcome = "Y", exposure = "A",
                               mediator = c("M1", "M2"), basec = c("C1", "C2"), EMint = TRUE,
                               mreg = list(mreg1, mreg2), yreg = "logistic",
                               astar = 0, a = 1, mval = list(1, 1),
                               estimation = "imputation", inference = "bootstrap")
  
  # test
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe), 
               unname(res_binbin_wb$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe)[c(6,2,5,15)], 
               unname(res_binbin_iorw$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe), 
               unname(res_binbin_ne$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe), 
               unname(res_binbin_msm$effect.pe), tolerance = 0.1)
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe), 
               unname(res_binbin_gformula$effect.pe), tolerance = 0.1)
  
})

test_that("cmest works correctly for binary Y and binary M", {
  
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 10000
  C1 <- rnorm(n, mean = 1, sd = 0.1)
  C2 <- rbinom(n, 1, 0.6)
  pa <- expit(0.2 + 0.5*C1 + 0.1*C2)
  A <- rbinom(n, 1, pa)
  pm <- expit(1 + 2*A + 1.5*C1 + 0.8*C2)
  M <- rbinom(n, 1, pm)
  py <- expit(-3 + 0.8*A - 1.8*M + 0.5*A*M + 0.3*C1 - 0.6*C2)
  Y <- rbinom(n, 1, py)
  data <- data.frame(A, M, Y, C1, C2)
  yreg <- glm(Y ~ A*M + C1 + C2, family = binomial(), data = data)
  mreg <- glm(M ~ A + C1 + C2, family = binomial(), data = data)
  
  # results of cmest
  res_binbin_rb_param_delta <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                                     mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                     mreg = list("logistic"), yreg = "logistic",
                                     astar = 0, a = 1, mval = list(1),
                                     estimation = "para", inference = "delt")
  res_binbin_rb_param_bootstrap <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                         mreg = list("logistic"), yreg = "logistic",
                                         astar = 0, a = 1, mval = list(1),
                                         estimation = "paramfunc", inference = "bootstrap")
  res_binbin_rb_impu_bootstrap <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                                        mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                        mreg = list("logistic"), yreg = "logistic",
                                        astar = 0, a = 1, mval = list(1),
                                        estimation = "impu", inference = "boot")
  res_binbin_wb <- cmest(data = data, model = "wb", outcome = "Y", exposure = "A",
                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                         ereg = "logistic", yreg = "logistic",
                         astar = 0, a = 1, mval = list(1),
                         estimation = "imputation", inference = "bootstrap")
  res_binbin_iorw <- cmest(data = data, model = "iorw", outcome = "Y", exposure = "A",
                           mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                           ereg = "logistic", yreg = "logistic",
                           astar = 0, a = 1, mval = list(1),
                           estimation = "imputation", inference = "bootstrap")
  res_binbin_msm <- cmest(data = data, model = "msm", outcome = "Y", exposure = "A",
                          mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                          ereg = "logistic", yreg = "logistic", mreg = list("logistic"),
                          wmdenomreg = list("logistic"), wmnomreg = list("logistic"),
                          astar = 0, a = 1, mval = list(1),
                          estimation = "imputation", inference = "bootstrap")
  res_binbin_ne <- cmest(data = data, model = "ne", outcome = "Y", exposure = "A",
                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                         yreg = "logistic",
                         astar = 0, a = 1, mval = list(1),
                         estimation = "imputation", inference = "bootstrap")
  res_binbin_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                               mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                               mreg = list("logistic"), yreg = "logistic",
                               astar = 0, a = 1, mval = list(1),
                               estimation = "imputation", inference = "bootstrap")
  
  # reference results
  thetas <- unname(coef(yreg))
  betas <- unname(coef(mreg))
  beta0 <- betas[1]
  beta1 <- betas[2]
  beta2 <- betas[3]
  beta3 <- betas[4]
  theta0 <- thetas[1]
  theta1 <- thetas[2]
  theta2 <- thetas[3]
  theta3 <- thetas[6]
  theta4 <- thetas[4]
  theta5 <- thetas[5]
  m <- 1
  a <- 1
  astar <- 0
  meanc1 <- mean(C1)
  meanc2 <- mean(C2)
  cde_binbin <- exp((theta1+theta3*m)*(a-astar))
  pnde_binbin <- (exp(theta1*a)*(1+exp(theta2+theta3*a+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))/
    (exp(theta1*astar)*(1+exp(theta2+theta3*astar+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))
  tnde_binbin <- (exp(theta1*a)*(1+exp(theta2+theta3*a+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))/
    (exp(theta1*astar)*(1+exp(theta2+theta3*astar+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))
  pnie_binbin <- ((1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*astar+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))/
    ((1+exp(beta0+beta1*a+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*astar+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))
  tnie_binbin <- ((1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*a+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))/
    ((1+exp(beta0+beta1*a+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*a+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))
  te_binbin <- pnde_binbin*tnie_binbin
  cde_err_binbin <- exp(theta2*m)*(exp(theta1*(a-astar)+theta3*a*m)-exp(theta3*astar*m))*
    (1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))/
    (1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2+theta2+theta3*astar))
  intref_err_binbin <- pnde_binbin - 1 - cde_err_binbin
  intmed_err_binbin <- tnie_binbin * pnde_binbin - pnde_binbin - pnie_binbin + 1
  pnie_err_binbin <- pnie_binbin - 1
  total_err_binbin <- te_binbin - 1
  cde_err_prop_binbin <- cde_err_binbin/total_err_binbin
  intmed_err_prop_binbin <- intmed_err_binbin/total_err_binbin
  intref_err_prop_binbin <- intref_err_binbin/total_err_binbin
  pnie_err_prop_binbin <- pnie_err_binbin/total_err_binbin
  pm_binbin <- (pnie_err_binbin+intmed_err_binbin)/total_err_binbin
  int_binbin <- (intref_err_binbin+intmed_err_binbin)/total_err_binbin
  pe_binbin <- (intref_err_binbin+intmed_err_binbin+pnie_err_binbin)/total_err_binbin
  
  ref <- c(cde_binbin, pnde_binbin, tnde_binbin, pnie_binbin, tnie_binbin, te_binbin, cde_err_binbin,
           intref_err_binbin, intmed_err_binbin, pnie_err_binbin, cde_err_prop_binbin,
           intref_err_prop_binbin, intmed_err_prop_binbin, pnie_err_prop_binbin,
           pm_binbin, int_binbin, pe_binbin)
  # test
  expect_equal(summary(res_binbin_rb_param_delta)$summarydf$Estimate, ref)
  expect_equal(class(ggcmest(res_binbin_rb_param_delta)), c("gg", "ggplot"))
  expect_equal(summary(res_binbin_rb_param_bootstrap)$summarydf$Estimate, ref)
  expect_equal(summary(res_binbin_rb_impu_bootstrap)$summarydf$Estimate, ref, tolerance = 0.1)
  expect_equal(summary(res_binbin_wb)$summarydf$Estimate, ref, tolerance = 0.1)
  expect_equal(summary(res_binbin_msm)$summarydf$Estimate, ref, tolerance = 0.1)
  expect_equal(summary(res_binbin_ne)$summarydf$Estimate, ref, tolerance = 0.1)
  expect_equal(summary(res_binbin_iorw)$summarydf$Estimate, ref[c(6,2,5,15)], tolerance = 0.1)
  expect_equal(summary(res_binbin_gformula)$summarydf$Estimate, ref, tolerance = 0.1)
  expect_equal(class(print(res_binbin_rb_param_delta)), "list")
  expect_equal(class(print(res_binbin_rb_param_bootstrap)), "list")
  expect_equal(class(print(res_binbin_rb_impu_bootstrap)), "list")
  expect_equal(class(print(res_binbin_wb)), "list")
  expect_equal(class(print(res_binbin_msm)), "list")
  expect_equal(class(print(res_binbin_ne)), "list")
  expect_equal(class(print(res_binbin_iorw)), "list")
  expect_equal(class(print(res_binbin_gformula)), "list")
  expect_equal(class(print(summary(res_binbin_rb_param_delta))), "list")
  expect_equal(class(print(summary(res_binbin_rb_param_bootstrap))), "list")
  expect_equal(class(print(summary(res_binbin_rb_impu_bootstrap))), "list")
  expect_equal(class(print(summary(res_binbin_wb))), "list")
  expect_equal(class(print(summary(res_binbin_msm))), "list")
  expect_equal(class(print(summary(res_binbin_ne))), "list")
  expect_equal(class(print(summary(res_binbin_iorw))), "list")
  expect_equal(class(print(summary(res_binbin_gformula))), "list")
  
})


test_that("cmest works correctly for binary Y and binary M in a case control study", {
  
  set.seed(1)
  # data simulation
  expit <- function(x) exp(x)/(1+exp(x))
  n <- 1000000
  C1 <- rnorm(n, mean = 1, sd = 0.1)
  C2 <- rbinom(n, 1, 0.6)
  pa <- expit(0.2 + 0.5*C1 + 0.1*C2)
  A <- rbinom(n, 1, pa)
  M <- rbinom(n, 1, expit(1 + 2*A + 1.5*C1 + 0.8*C2))
  py <- expit(-3 + 0.8*A - 1.8*M + 0.5*A*M + 0.3*C1 - 0.6*C2)
  Y <- rbinom(n, 1, py)
  yprevalence <- sum(Y)/n
  data <- data.frame(A, M, Y, C1, C2)
  case_indice <- sample(which(data$Y == 1), 2000, replace = FALSE)
  control_indice <- sample(which(data$Y == 0), 2000, replace = FALSE)
  data <- data[c(case_indice, control_indice), ]
  
  # yrare = TRUE
  # results of cmest
  res_binbin_rb_param_delta <- cmest(data = data, model = "rb", casecontrol = TRUE, yrare = TRUE,
                                     outcome = "Y", exposure = "A",
                                     mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                     mreg = list("logistic"), yreg = "loglinear",
                                     astar = 0, a = 1, mval = list(1),
                                     estimation = "paramfunc", inference = "delta")
  res_binbin_rb_param_bootstrap <- cmest(data = data, model = "rb", casecontrol = TRUE, yrare = TRUE,
                                         outcome = "Y", exposure = "A",
                                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                         mreg = list("logistic"), yreg = "loglinear",
                                         astar = 0, a = 1, mval = list(1),
                                         estimation = "paramfunc", inference = "bootstrap")
  res_binbin_rb_impu_bootstrap <- cmest(data = data, model = "rb", casecontrol = TRUE, yrare = TRUE,
                                        outcome = "Y", exposure = "A",
                                        mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                        mreg = list("logistic"), yreg = "logistic",
                                        astar = 0, a = 1, mval = list(1),
                                        estimation = "imputation", inference = "bootstrap")
  res_binbin_wb <- cmest(data = data, model = "wb", casecontrol = TRUE, yrare = TRUE,
                         outcome = "Y", exposure = "A",
                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                         ereg = "logistic", yreg = "logistic",
                         astar = 0, a = 1, mval = list(1),
                         estimation = "imputation", inference = "bootstrap")
  res_binbin_iorw <- cmest(data = data, model = "iorw", casecontrol = TRUE, yrare = TRUE,
                           outcome = "Y", exposure = "A",
                           mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                           ereg = "logistic", yreg = "logistic",
                           astar = 0, a = 1, mval = list(1),
                           estimation = "imputation", inference = "bootstrap")
  res_binbin_msm <- cmest(data = data, model = "msm", casecontrol = TRUE, yrare = TRUE,
                          outcome = "Y", exposure = "A",
                          mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                          ereg = "logistic", yreg = "logistic", mreg = list("logistic"),
                          wmnomreg = list("logistic"), wmdenomreg = list("logistic"),
                          astar = 0, a = 1, mval = list(1),
                          estimation = "imputation", inference = "bootstrap")
  res_binbin_ne <- cmest(data = data, model = "ne", casecontrol = TRUE, yrare = TRUE,
                         outcome = "Y", exposure = "A",
                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                         yreg = "logistic",
                         astar = 0, a = 1, mval = list(1),
                         estimation = "imputation", inference = "bootstrap")
  res_binbin_gformula <- cmest(data = data, model = "gformula", casecontrol = TRUE, yrare = TRUE,
                               outcome = "Y", exposure = "A",
                               mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                               mreg = list("logistic"), yreg = "logistic",
                               astar = 0, a = 1, mval = list(1),
                               estimation = "imputation", inference = "bootstrap")
  
  # reference results
  thetas <- unname(coef(res_binbin_rb_param_delta$reg.output$yreg))
  betas <- unname(coef(res_binbin_rb_param_delta$reg.output$mreg[[1]]))
  beta0 <- betas[1]
  beta1 <- betas[2]
  beta2 <- betas[3]
  beta3 <- betas[4]
  theta0 <- thetas[1]
  theta1 <- thetas[2]
  theta2 <- thetas[3]
  theta3 <- thetas[6]
  theta4 <- thetas[4]
  theta5 <- thetas[5]
  m <- 1
  a <- 1
  astar <- 0
  meanc1 <- mean(data$C1)
  meanc2 <- mean(data$C2)
  cde_binbin <- exp((theta1+theta3*m)*(a-astar))
  pnde_binbin <- (exp(theta1*a)*(1+exp(theta2+theta3*a+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))/
    (exp(theta1*astar)*(1+exp(theta2+theta3*astar+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))
  tnde_binbin <- (exp(theta1*a)*(1+exp(theta2+theta3*a+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))/
    (exp(theta1*astar)*(1+exp(theta2+theta3*astar+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))
  pnie_binbin <- ((1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*astar+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))/
    ((1+exp(beta0+beta1*a+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*astar+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))
  tnie_binbin <- ((1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*a+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))/
    ((1+exp(beta0+beta1*a+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*a+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))
  te_binbin <- pnde_binbin*tnie_binbin
  cde_err_binbin <- exp(theta2*m)*(exp(theta1*(a-astar)+theta3*a*m)-exp(theta3*astar*m))*
    (1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))/
    (1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2+theta2+theta3*astar))
  intref_err_binbin <- pnde_binbin - 1 - cde_err_binbin
  intmed_err_binbin <- tnie_binbin * pnde_binbin - pnde_binbin - pnie_binbin + 1
  pnie_err_binbin <- pnie_binbin - 1
  total_err_binbin <- te_binbin - 1
  cde_err_prop_binbin <- cde_err_binbin/total_err_binbin
  intmed_err_prop_binbin <- intmed_err_binbin/total_err_binbin
  intref_err_prop_binbin <- intref_err_binbin/total_err_binbin
  pnie_err_prop_binbin <- pnie_err_binbin/total_err_binbin
  pm_binbin <- (pnie_err_binbin+intmed_err_binbin)/total_err_binbin
  int_binbin <- (intref_err_binbin+intmed_err_binbin)/total_err_binbin
  pe_binbin <- (intref_err_binbin+intmed_err_binbin+pnie_err_binbin)/total_err_binbin
  
  ref <- c(cde_binbin, pnde_binbin, tnde_binbin, pnie_binbin, tnie_binbin, te_binbin, cde_err_binbin,
           intref_err_binbin, intmed_err_binbin, pnie_err_binbin, cde_err_prop_binbin,
           intref_err_prop_binbin, intmed_err_prop_binbin, pnie_err_prop_binbin,
           pm_binbin, int_binbin, pe_binbin)
  # test
  expect_equal(unname(res_binbin_rb_param_delta$effect.pe), ref)
  expect_equal(unname(res_binbin_rb_param_bootstrap$effect.pe), ref)
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe), ref, tolerance = 0.1)
  expect_equal(unname(res_binbin_wb$effect.pe), ref, tolerance = 0.1)
  expect_equal(unname(res_binbin_msm$effect.pe), ref, tolerance = 0.1)
  expect_equal(unname(res_binbin_ne$effect.pe), ref, tolerance = 0.1)
  expect_equal(unname(res_binbin_iorw$effect.pe), ref[c(6,2,5,15)], tolerance = 0.1)
  expect_equal(unname(res_binbin_gformula$effect.pe), ref, tolerance = 0.1)
  
  # yprevalence = yprevalence
  # results of cmest
  res_binbin_rb_param_delta <- cmest(data = data, model = "rb",
                                     casecontrol = TRUE, yprevalence = yprevalence,
                                     outcome = "Y", exposure = "A",
                                     mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                     mreg = list("logistic"), yreg = "loglinear",
                                     astar = 0, a = 1, mval = list(1),
                                     estimation = "paramfunc", inference = "delta")
  res_binbin_rb_param_delta_multinom <- cmest(data = data, model = "rb",
                                              casecontrol = TRUE, yprevalence = yprevalence,
                                              outcome = "Y", exposure = "A",
                                              mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                              mreg = list("multinomial"), yreg = "loglinear",
                                              astar = 0, a = 1, mval = list(1),
                                              estimation = "paramfunc", inference = "delta", 
                                              multimp = TRUE, m = 1)
  res_binbin_rb_param_bootstrap <- cmest(data = data, model = "rb",
                                         casecontrol = TRUE, yprevalence = yprevalence,
                                         outcome = "Y", exposure = "A",
                                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                         mreg = list("logistic"), yreg = "loglinear",
                                         astar = 0, a = 1, mval = list(1),
                                         estimation = "paramfunc", inference = "bootstrap")
  res_binbin_rb_impu_bootstrap <- cmest(data = data, model = "rb",
                                        casecontrol = TRUE, yprevalence = yprevalence,
                                        outcome = "Y", exposure = "A",
                                        mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                                        mreg = list("logistic"), yreg = "logistic",
                                        astar = 0, a = 1, mval = list(1),
                                        estimation = "imputation", inference = "bootstrap")
  res_binbin_wb <- cmest(data = data, model = "wb", casecontrol = TRUE, yprevalence = yprevalence,
                         outcome = "Y", exposure = "A",
                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                         ereg = "logistic", yreg = "logistic",
                         astar = 0, a = 1, mval = list(1),
                         estimation = "imputation", inference = "bootstrap")
  res_binbin_iorw <- cmest(data = data, model = "iorw", casecontrol = TRUE, yprevalence = yprevalence,
                           outcome = "Y", exposure = "A",
                           mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                           ereg = "logistic", yreg = "logistic",
                           astar = 0, a = 1, mval = list(1),
                           estimation = "imputation", inference = "bootstrap")
  res_binbin_msm <- cmest(data = data, model = "msm", casecontrol = TRUE, yprevalence = yprevalence,
                          outcome = "Y", exposure = "A",
                          mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                          ereg = "logistic", yreg = "logistic", mreg = list("logistic"),
                          wmnomreg = list("logistic"), wmdenomreg = list("logistic"),
                          astar = 0, a = 1, mval = list(1),
                          estimation = "imputation", inference = "bootstrap")
  res_binbin_ne <- cmest(data = data, model = "ne", casecontrol = TRUE, yprevalence = yprevalence,
                         outcome = "Y", exposure = "A",
                         mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                         yreg = "logistic",
                         astar = 0, a = 1, mval = list(1),
                         estimation = "imputation", inference = "bootstrap")
  res_binbin_gformula <- cmest(data = data, model = "gformula", casecontrol = TRUE, yprevalence = yprevalence,
                               outcome = "Y", exposure = "A",
                               mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                               mreg = list("logistic"), yreg = "logistic",
                               astar = 0, a = 1, mval = list(1),
                               estimation = "imputation", inference = "bootstrap")
  
  # reference results
  thetas <- unname(coef(res_binbin_rb_param_delta$reg.output$yreg))
  betas <- unname(coef(res_binbin_rb_param_delta$reg.output$mreg[[1]]))
  beta0 <- betas[1]
  beta1 <- betas[2]
  beta2 <- betas[3]
  beta3 <- betas[4]
  theta0 <- thetas[1]
  theta1 <- thetas[2]
  theta2 <- thetas[3]
  theta3 <- thetas[6]
  theta4 <- thetas[4]
  theta5 <- thetas[5]
  m <- 1
  a <- 1
  astar <- 0
  meanc1 <- mean(data$C1)
  meanc2 <- mean(data$C2)
  cde_binbin <- exp((theta1+theta3*m)*(a-astar))
  pnde_binbin <- (exp(theta1*a)*(1+exp(theta2+theta3*a+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))/
    (exp(theta1*astar)*(1+exp(theta2+theta3*astar+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))
  tnde_binbin <- (exp(theta1*a)*(1+exp(theta2+theta3*a+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))/
    (exp(theta1*astar)*(1+exp(theta2+theta3*astar+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))
  pnie_binbin <- ((1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*astar+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))/
    ((1+exp(beta0+beta1*a+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*astar+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))
  tnie_binbin <- ((1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*a+beta0+beta1*a+beta2*meanc1+beta3*meanc2)))/
    ((1+exp(beta0+beta1*a+beta2*meanc1+beta3*meanc2))*(1+exp(theta2+theta3*a+beta0+beta1*astar+beta2*meanc1+beta3*meanc2)))
  te_binbin <- pnde_binbin*tnie_binbin
  cde_err_binbin <- exp(theta2*m)*(exp(theta1*(a-astar)+theta3*a*m)-exp(theta3*astar*m))*
    (1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2))/
    (1+exp(beta0+beta1*astar+beta2*meanc1+beta3*meanc2+theta2+theta3*astar))
  intref_err_binbin <- pnde_binbin - 1 - cde_err_binbin
  intmed_err_binbin <- tnie_binbin * pnde_binbin - pnde_binbin - pnie_binbin + 1
  pnie_err_binbin <- pnie_binbin - 1
  total_err_binbin <- te_binbin - 1
  cde_err_prop_binbin <- cde_err_binbin/total_err_binbin
  intmed_err_prop_binbin <- intmed_err_binbin/total_err_binbin
  intref_err_prop_binbin <- intref_err_binbin/total_err_binbin
  pnie_err_prop_binbin <- pnie_err_binbin/total_err_binbin
  pm_binbin <- (pnie_err_binbin+intmed_err_binbin)/total_err_binbin
  int_binbin <- (intref_err_binbin+intmed_err_binbin)/total_err_binbin
  pe_binbin <- (intref_err_binbin+intmed_err_binbin+pnie_err_binbin)/total_err_binbin
  
  ref <- c(cde_binbin, pnde_binbin, tnde_binbin, pnie_binbin, tnie_binbin, te_binbin, cde_err_binbin,
           intref_err_binbin, intmed_err_binbin, pnie_err_binbin, cde_err_prop_binbin,
           intref_err_prop_binbin, intmed_err_prop_binbin, pnie_err_prop_binbin,
           pm_binbin, int_binbin, pe_binbin)
  # test
  expect_equal(unname(res_binbin_rb_param_delta$effect.pe), ref)
  expect_equal(unname(res_binbin_rb_param_delta_multinom$effect.pe), ref, tolerance = 0.001)
  expect_equal(class(print(res_binbin_rb_param_delta)), "list")
  expect_equal(class(print(res_binbin_rb_param_delta_multinom)), "list")
  expect_equal(unname(res_binbin_rb_param_bootstrap$effect.pe), ref)
  expect_equal(unname(res_binbin_rb_impu_bootstrap$effect.pe), ref, tolerance = 0.1)
  expect_equal(unname(res_binbin_wb$effect.pe), ref, tolerance = 0.1)
  expect_equal(unname(res_binbin_msm$effect.pe), ref, tolerance = 0.1)
  expect_equal(unname(res_binbin_ne$effect.pe), ref, tolerance = 0.1)
  expect_equal(unname(res_binbin_iorw$effect.pe), ref[c(6,2,5,15)], tolerance = 0.1)
  expect_equal(unname(res_binbin_gformula$effect.pe), ref, tolerance = 0.1)
  expect_equal(class(print(summary(res_binbin_rb_param_delta_multinom))), "list")
  
})