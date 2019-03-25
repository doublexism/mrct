## functions
## Non-parametric bootstrap
## simulation
library(plyr)
library(dplyr)
library(purrr)
library(survival)
library(ggplot2)
library(grid)
library(cowplot)
library(gsDesign)
library(doParallel)
## survival cut-off
survCutOff <- function(data, cutoff){
  data$event <- as.numeric(data$event == 1 & data$stop <= cutoff)
  data$stop <- if_else(data$stop < cutoff, data$stop, cutoff)
  data$tte <- data$stop - data$start
  return(data)
}

## generate parameter setting 
paramGenerate <- function(num_pat, num_region,prop_region,type = "mean",num_arm = 2,ratio = NULL, num_param = NULL,name_param = NULL, par1_arm1 =NULL, par1_arm2=NULL, more_par=NULL){
  ## The function is to generate parameters setting for simulation. 
  ## The value of function is a dataframe with 6 colums (Region, Arm, Param, Value of param, Proportion of region and Num of patients per arm and region)
  
  # num_pat: number of patients in the study
  # num_region: number of regions
  # prop_region: proportion of each region
  # num_arm: number of treatment arms
  # ratio: the randomization ratio
  # num_param: number of parameters; 2 for continous(Mean and SD), 1 for binary (P), 
  # 6 for survival (enrollment rate, shape parameter for event, scale parameter for event, shape parameter of censor, scale parameter for censor and max followup period)
  # name_parameter: optional, name of the parameters
  # par1_arm1: a vector containing the value of parameter1 in the first arm across regions
  # par1_arm2: a vector containing the value of parameter1 in the second arm across regions
  # more_par: a list of vectors that contain additional parameters.
  
  if (is.null(ratio)){
    ratio <- rep(1,num_arm)
  }
  params <- unlist(c(par1_arm1, par1_arm2, more_par))
  param_dat <- expand.grid(1:num_region, 1:num_arm, 1:num_param)
  param_dat$value <- params
  quantum <- kronecker(ratio, prop_region)/sum(ratio)
  param_dat$proportion <- rep(quantum,num_param)
  param_dat$num <- round(num_pat*param_dat$proportion)
  colnames(param_dat) <- c("Region", "Arm", "Param","Value","Prop","Num")
  param_dat$Param <- paste0("Arg-",param_dat$Param)
  param_dat$type <- type
  return(param_dat)
}

sim.continuous <- function(par, dist = "normal"){
  ## internal function to simulate continous outome
  ## return a data.frame for given parameter subset
  
  # par, paramter
  # dist, distributions
  num <- par$Num[1]
  mu_sigma <- par$Value
  if (dist == "normal"){
    val <- rnorm(num,mean =   mu_sigma[1], sd = mu_sigma[2])
  } else {
    stop("Distribution not yet supported") 
  }
  return(data.frame(value = val))
} 

sim.binary <- function(par){
  ## internal function to simulate binary outome
  ## return a data.frame for given parameter subset
  
  # par, paramter
  num <- par$Num[1]
  p <- par$Value
  val <- rbinom(num, 1, p)
  return(data.frame(value = val))
}

sim.survival <- function(par){
  ## internal function to simulate surviavl data
  # six parameters in total for weibull survival:
  # enrollment rate
  # median parameter for the event
  # median parameter for the censoring
  # maximum length of follow-up
  num <- par$Num[1]
  par_event <- par$Value[2]
  par_censor <- par$Value[3]
  par_enr <- par$Value[1]
  par_fol <- par$Value[4]
  par_events <- par$Value[5]
  enroll_time <- round(seq(0, num/par_enr, length.out = num)*30)
  event_free_time <- round(rexp(num, par_event)*30) + enroll_time
  censor_time <- round(rexp(num, par_censor)*30) + enroll_time
  EOS <- pmin(event_free_time, censor_time)
  event <- as.numeric(EOS == event_free_time)
  return(data.frame(start = enroll_time, stop = EOS, tte = EOS - enroll_time,event = event, weight = 1))
}

simulate <- function(params, sim_func){
  ## the main function to simulate dataset
  # params: the parameter generated
  # sim_func: the internal simulation function to use
  value <- params %>% split(.,list(.$Region,.$Arm)) %>% map_dfr(sim_func)
  dat <- distinct(params[c("Region","Arm","Num")])
  dat <- params[rep(1:nrow(dat),dat$Num),c("Region","Arm")]
  dat <- bind_cols(dat, value)
  if (params$type[1] == "survival"){
    par_event <- tail(params$Value,1)
    par_fol <- tail(params$Value[params$Param == "Arg-4"],1)
    if (par_fol == 0){
    study_time <- sort(dat$stop[dat$event == 1])[par_event]
    dat <- survCutOff(dat, study_time)
    } else {
      dat <- survCutOff(dat, par_fol)
    }
  }

  dat <- structure(dat, class = c("sims","data.frame"))
  return(dat)
}
# bootstrap
# non-parametric bootstrap
# multinomial sampling prob
npProbability <- function(data){
  ## the non-parametric bootstrap probability for each subject in the dataset
  return(rep(1/nrow(data), nrow(data)))
}

bootWeighting <- function(num_boot, data,sample_size, probability = npProbability){
  ## generate bootstrap result
  ## returns a matrix (num_row * num_boot) representing the bootstrap result
  # num_boot: the number of bootstrap resampling
  # data: the dataframe to resample
  # sample_size: the size of resampling
  # probability: the probability of resampling for each row in the data
  prob <- probability(data)
  mat <- rmultinom(num_boot, sample_size, prob)
  return(mat)
}
# weighted.mean.difference
boot_t_ci <- function(dat, weight = 1, control = 2, treatment =1){
  ## calculate bootstrap for t statistics
  ## return a dataframe  
  index <- dat$Arm == treatment
  if (length(weight) <= 1){
    weight_trt <- matrix(rep(weight, sum(index)))
    weight_ctrl <- matrix(rep(weight, sum(!index)))
  } else{
    weight_trt <- weight[index,]
    weight_ctrl <- weight[!index,]
  }
  N_trt <- colSums(weight_ctrl)
  N_ctrl <- colSums(weight_ctrl)
  Mean_trt <- c(dat$value[index] %*% weight_trt)/N_trt 
  Mean_ctrl <- c(dat$value[!index] %*% weight_ctrl)/N_ctrl
  Var <- (colSums((dat$value[index] * weight_trt - Mean_trt* weight_trt)**2)+ colSums((dat$value[!index]*weight_ctrl - Mean_ctrl*weight_ctrl)**2))/(N_trt+N_ctrl-2)
  se_diff <- (sqrt(Var/N_trt + Var/N_ctrl))
  df <- N_trt + N_ctrl - 2
  CI_low <- c(Mean_trt -Mean_ctrl + se_diff * qt(0.025, df))
  CI_high <- c(Mean_trt -Mean_ctrl+ se_diff * qt(0.975, df))
  return(data.frame(mean_diff = Mean_trt -Mean_ctrl,low_CI = CI_low, high_CI = CI_high))
}

bootDiffMean <- function(dat, weight = 1, control = 2, treatment =1){
  ## bootstrap for the mean difference for the continous variable
  index <- 2*(dat$Arm == treatment) - 1
  weight[index == 1,] <- weight[index == 1,]%*% diag(1/colSums(weight[index==1,]))
  weight[index == -1,] <- weight[index == -1,] %*% diag(1/colSums(weight[index==-1,]))
  weight <- weight * index
  boot_diff <- c(dat$value %*% weight)
  return(boot_diff)
}

bootHR <- function(dat, weight, control, treatment){
  weight.list <- split(weight,col(weight))
  safe_cox <- safely(coxph, otherwise = NULL)
  result <- map(weight.list, 
                    ~coxph(Surv(tte, event) ~ Arm, data = dat,weights = .,subset = . > 0)$coefficients)
  result <- unlist(result)
  return(result)
}
# calculate 95% bootstrap CI overlap
bootCover <- function(LCI, NLCI, upper = FALSE){
  ## calculate the coverage rate based on bootstrap CI for Local and non-Local
  # LCI: bootstrap CI for Local population
  # NLCI: bootstrap CI for the Non-Local population
  if (upper == FALSE){
    bootCover <- max((min(LCI[2],NLCI[2]) - max(LCI[1],NLCI[1]))/(LCI[2] - LCI[1]),0, na.rm = TRUE)
  } else {
    bootCover <- max((LCI[2] - max(LCI[1],NLCI[1]))/(LCI[2] - LCI[1]),0, na.rm = TRUE)
  }
  return(bootCover)
}

## calculate BCICR for a simulated dataset
BCICR <- function(dat, local_index, treatment = 1, control=2, num_boot = 1000,statistics = "mean", method = "np", probability = npProbability){
  ## calculate BCICR for a simulated dataset
  ## returns a list containing the BCICR , local CI and non-local CI
  # dat: the dataset
  # local_index: the index of local population
  # treatment: the index of treatment arm
  # control: the index of comparator
  # num_boot: the number of bootstrap resampling
  # statistics: the statistics to resample
  # method: the method of resampling
  # probability: the function to calculate resamping probability
  dat_local <- dat[dat$Region == local_index,]
  dat_nonlocal <- dat[dat$Region != local_index,]
  if (statistics == "mean"){
    global_test <- t.test(dat$value[dat$Arm == treatment], dat$value[dat$Arm == control],alternative = "greater")[c("statistic","p.value")]
    pass <- (global_test$statistic > 0 & global_test$p.value < 0.025)
  } else if (statistics == "proportion"){

    n1 <- sum(dat$Arm == treatment)
    n2 <- sum(dat$Arm == control)
    p1 <- sum(dat$value[dat$Arm == treatment])/n1
    p2 <- sum(dat$value[dat$Arm == control])/n2
    z <- (p1-p2)/sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)
    pass <- z > 1.96
  } else if (statistics == "loghr"){
    chisq <- survdiff(Surv(tte, event)~Arm, data = dat)$chisq
    pass <- chisq > qchisq(0.975,1)
  }
  if (pass == FALSE){
    return(list(BCICR = NULL, BCICR_m = NULL, local_CI = NULL, nonlocal_CI = NULL))
  }
  if (statistics %in% c("mean", "proportion")){
    if (method == "np"){
      local_weight <- bootWeighting(num_boot,data = dat_local, sample_size = nrow(dat_local),probability = probability)
      local_diff <- bootDiffMean(dat_local, weight = local_weight,control = control, treatment = treatment)
      local_CI <- quantile(local_diff,probs = c(0.025, 0.975),na.rm = TRUE)
      nonlocal_weight <- bootWeighting(num_boot,data = dat_nonlocal, sample_size = nrow(dat_local),probability = probability)
      nonlocal_diff <- bootDiffMean(dat_nonlocal, weight = nonlocal_weight,control = control, treatment = treatment)
      nonlocal_CI <- quantile(nonlocal_diff,probs = c(0.025, 0.975),na.rm=TRUE)
    }
  } else if (statistics == "loghr"){
    local_weight <- bootWeighting(num_boot,data = dat_local, sample_size = nrow(dat_local),probability = probability)
    local_diff <- bootHR(dat_local, weight = local_weight,control = control, treatment = treatment)
    local_CI <- quantile(local_diff,probs = c(0.025, 0.975),na.rm = TRUE)
    nonlocal_weight <- bootWeighting(num_boot,data = dat_nonlocal, sample_size = nrow(dat_local),probability = probability)
    nonlocal_diff <- bootHR(dat_nonlocal, weight = nonlocal_weight,control = control, treatment = treatment)
    nonlocal_CI <- quantile(nonlocal_diff,probs = c(0.025, 0.975),na.rm=TRUE)
    print(nonlocal_CI)
    print(local_CI)
  }
  bcicr <- bootCover(local_CI, nonlocal_CI)
  bcicr_m <- bootCover(local_CI, nonlocal_CI, upper = TRUE)
  return(list(BCICR = bcicr, BCICR_m = bcicr_m, local_CI = local_CI, nonlocal_CI = nonlocal_CI))
}

LPP <- function(bcicr, cutoff){
  ## return the LPP for BCICR under certain cutoff
  # bcicr: a vector of bcicr
  # cutoff: cutoff for bcicr
  return(sum(bcicr > cutoff)/length(bcicr))
}

BCICR_sim <- function(params,cutoff,local, num_sim = 1000, num_boot=1000,probability = npProbability, statistics = "mean",sim_func = sim.continuous,treatment = 1, control=2,seed = 0){
  ## the main simulation function for a fixed parameters
  ## returns a list of BCICR, cutoff and LPP for each simulated dataset
  # params: the generated parameters
  # cutoff: a vector of cutoff values to explore
  # local: the index of local population
  # num_sim: number of simulated datasets
  # num_boot: number of bootstrap simulation
  # probability: the function to calculate the resampling probability
  # statistics: the statistics to calculate
  # sim_func: the simulation function
  # treatment: the index of treatment arm
  # control :the index of control arm
  # seed: random seed
  set.seed(seed)
  num_group <- nrow(params)/length(unique(params$Param))
  num_per_group <- params$Num[1:num_group]
  start_data <- cumsum(c(1, num_per_group[-num_group]*num_sim))
  end_data <- start_data + num_per_group-1
  index <- map2(start_data, end_data, seq,by = 1) %>% unlist()
  if (statistics == "loghr"){
    dat <- rerun(num_sim, simulate(params = params,sim_func = sim_func))
    bcicr_result <- map(dat, BCICR, local_index = local,treatment = treatment, control = control,num_boot = num_boot,statistics = statistics,probability = probability)
  } else {
    params$Num <- params$Num * num_sim
    dat <- simulate(params = params,sim_func = sim_func)
    bcicr_result <- map(1:num_sim, ~BCICR(dat = dat[index + (.-1)*rep(num_per_group,num_per_group),], local_index = local,treatment = treatment, control = control,num_boot = num_boot,statistics = statistics,probability = probability))
  }
  bcicr <- getListElement(bcicr_result, "BCICR")
  bcicr_m <- getListElement(bcicr_result, "BCICR_m")
  LPP <- map_dbl(cutoff, ~LPP(bcicr, .))
  LPP_m <- map_dbl(cutoff, ~LPP(bcicr_m,.))
  return(list(BCICR = bcicr,BCICR_m = bcicr_m,cutoff=cutoff,LPP = LPP, LPP_m = LPP_m))
}


#sim scenarios
scenarios_sim <- function(prop_local, sample_size, delta,num_sim = 1000, num_boot = 1000,parallel = FALSE,
                          type = "continuous", mean_diff = 1,mean_variance = 2,unequal_var = FALSE, unequal_eff = FALSE, NL_var = NULL,NL_delta = NULL,
                          treatment_rate = 0.3,control_rate = 0.1,
                          expect_hr = 0.7,enr = 30,median_tte_control = 20,censor_tte = 100,follow_up=0, seed = 20190117){
  ## perform simulation scenarios for different scenarios
  # prop_local: a vector of local population proportions
  # sample_size: a vector of the number of patients, if sample_size  = NULL, it will be calculated by power.
  # delta: the difference of effect size in percentage 0 to 1
  # unequal_var: the variance between multiple non-locol population are not equal
  # unequal_effect: the effects in non-local population are not equal
  if (parallel == TRUE){
    cl <- makeCluster(detectCores(logical = FALSE))
    registerDoParallel(cl)
  }
  scenarios <- expand.grid(sample_size = sample_size, prop_local = prop_local, delta=delta)
  if (type == "continuous"){
    
    NL_num <- ifelse(any(unequal_eff, unequal_var), max(length(NL_var), length(NL_delta)),1)
    if (unequal_var == FALSE){
      NL_var <- rep(mean_variance,NL_num)
    } 
    if (unequal_eff == FALSE){
      NL_delta <- rep(mean_diff,NL_num)
    }
    prop_region <-  map(prop_local,~c(rep((1 - .)/NL_num,NL_num),.)) %>% do.call(rbind,.)
    if (sample_size == 0){
      mean_var <-  c(prop_region %*% c(NL_var,mean_variance))**2
      sample_size <- ceiling((qnorm(0.025)+qnorm(0.1))**2*2*mean_var/mean_diff**2)*2
      scenarios$sample_size <- rep(sample_size, length(delta))
    }
    NL_delta <- NL_num*mean_diff/sum(NL_delta) * NL_delta
    params.list <- pmap(scenarios, ~paramGenerate(num_pat=..1,num_region = NL_num+1,num_param = 2,prop_region = c(rep((1 - ..2)/NL_num,NL_num),..2),
                                                  par1_arm1 = c(NL_delta,0+mean_diff*..3),par1_arm2 = rep(0,NL_num+1), more_par = list(par2_arm1 = c(NL_var,mean_variance), par2_arm2 = c(NL_var,mean_variance))))
    sim <- llply(params.list, 
                 BCICR_sim,
                 cutoff = seq(0,1,0.01),
                 local =NL_num+1,
                 num_sim = num_sim,
                 num_boot = num_boot,
                 sim_func = sim.continuous,
                 statistics = "mean",
                 seed = seed,
                 .parallel = parallel,
                 .paropts = list(.packages = c("plyr","dplyr","purrr","gsDesign","survival"),
                                 .export = as.vector(lsf.str(envir = .GlobalEnv)) ))
    scenarios$sample_size <- map_dbl(params.list, ~sum(.$Num)/length(unique(.$Param)))
    scenarios$variable_type <- type
    scenarios$unequal_variance <- unequal_var
    scenarios$unequal_effect <- unequal_eff
  } else if (type == "binary"){
    if (sample_size == 0){
      mean_diff = treatment_rate-control_rate
      sample_size <- ceiling((qnorm(0.025)+qnorm(0.1))**2*(treatment_rate*(1 - treatment_rate)+control_rate*(1-control_rate))/mean_diff**2)*2
      scenarios$sample_size <- sample_size
    }
    params.list <- pmap(scenarios, ~paramGenerate(num_pat=..1,num_region =2,num_param = 1,prop_region = c(1 - ..2,..2),
                                                  par1_arm1 = c(treatment_rate,(treatment_rate -control_rate)*..3+control_rate),
                                                  par1_arm2 = rep(control_rate,2)))
    sim <- llply(params.list,
                 BCICR_sim,
                 cutoff = seq(0,1,0.01),
                 local =2,
                 num_sim = num_sim,
                 num_boot = num_boot,
                 sim_func = sim.binary,
                 statistics = "proportion",
                 seed = seed,
                 .parallel = parallel,
                 .paropts = list(.packages = c("plyr","dplyr","purrr","gsDesign","survival"),
                                 .export = as.vector(lsf.str(envir = .GlobalEnv))
                                 ))
    scenarios$sample_size <- map_dbl(params.list, ~sum(.$Num))
    scenarios$variable_type <- type
    scenarios$treatment_rate <- treatment_rate
    scenarios$control_rate <- control_rate
   } else if (type == "survival"){
    if (sample_size == 0){
      design <- gsSurv(k=2, timing = 0.67,test.type = 4,sfu = sfLDOF,sfl = sfHSD, sflpar = -10,lambdaC = log(2)/median_tte_control, hr = expect_hr, hr0=1,eta = log(2)/censor_tte,gamma = enr,minfup = median_tte_control)
      sample_size <- round(c(tail(design$eNC,1)) + c(tail(design$eNE,1)))
      scenarios$sample_size <- sample_size
      events <- round(c(tail(design$eDC,1)) + c(tail(design$eDE,1)))
    }
      params.list <- pmap(scenarios, ~paramGenerate(num_pat=..1,
                                                    num_region =2,
                                                    num_param = 5,
                                                    type = "survival",
                                                    prop_region = c(1 - ..2,..2),
                                                    par1_arm1 = c(enr*c((1-..2),..2)),
                                                    par1_arm2 = c(enr*c((1-..2),..2)),
                                                    more_par = list(
                                                    par2_arm1 = c(log(2)*expect_hr/median_tte_control,log(2)*(expect_hr**..3)/median_tte_control),
                                                    par2_arm2 = c(log(2)/median_tte_control,log(2)/median_tte_control),
                                                    par3_arm1 = c(log(2)/censor_tte, log(2)/censor_tte),
                                                    par3_arm2 = c(log(2)/censor_tte, log(2)/censor_tte),
                                                    par4_arm1 = c(follow_up, follow_up),
                                                    par4_arm2 = c(follow_up, follow_up),
                                                    par5_arm1 = c(events,events),
                                                    par5_arm2 = c(events,events))
                                                    ))
      sim <- llply(params.list, 
                   BCICR_sim,
                   cutoff = seq(0,1,0.01),
                   local =2,
                   num_sim = num_sim,
                   num_boot = num_boot,
                   sim_func = sim.survival,
                   statistics = "loghr",
                   seed = seed,
                   .parallel = parallel,
                   .paropts = list(.packages = c("plyr","dplyr","purrr","gsDesign","survival"),
                                   .export = as.vector(lsf.str(envir = .GlobalEnv))
                   ))
      scenarios$variable_type <- type
      scenarios$control_median <- median_tte_control
      scenarios$hr <- expect_hr
      scenarios$censor_median <- censor_tte
      scenarios$follow_up <- median_tte_control
  }
  if (parallel == TRUE){
    stopCluster(cl)
  }
  return(list(scenarios = scenarios,sim_result = sim))
}


# simulated data diagnostics
summary.sims <- function(sims,...){
  # summary the simulated dataset
  # sims: simulated dataset
  # ...: summary statistics quoted, such as "mean", "sd",...
  statistics <- list(...)
  dfs <- sims %>% split(list(.$Region,.$Arm)) 
  value <- map(dfs,~`[[`(.,"value"))
  check <- map(dfs, ~distinct(.[c("Region","Arm")])) %>% bind_rows()
  stats <- list()
  for (s in statistics){
    .func <- match.fun(s)
    stats[[s]] <- map_dbl(value, .func)
  }
  check <- bind_cols(check, as.data.frame(stats))
  return(check)
}
# get element from a list
getListElement <- function(.l, name, simplify = TRUE){
  element <- map(.l,`[[`,name)
  if (simplify){
    element <- unlist(element)
  }
  return(element)
}
# summarizing and ploting 
# summarizing the outputs
# summarizing the outputs
simResult <- function(...){
  results <- list(...)
  scenarios <- getListElement(results,name = "scenarios",simplify = FALSE) %>% bind_rows()
  scenarios <- scenarios[rep(1:nrow(scenarios), each = 101),]
  result.list <- getListElement(results,name = "sim_result",simplify = FALSE) %>% 
    do.call(c,.)
  LPP <- result.list %>% 
    getListElement(., name = "LPP")
  LPP_m <- result.list %>% 
    getListElement(., name = "LPP_m")
  cutoff <- result.list %>% 
    getListElement(., name = "cutoff")
  LPP_dat <- data.frame(cutoff = cutoff, LPP = LPP,LPP_m = LPP_m)
  data <- bind_cols(scenarios, LPP_dat)
  # get raw bcicr and bcicr_m data
  data_bcicr <-getListElement(result.list, name = "BCICR", simplify = FALSE)
  data_bcicr_m <- getListElement(result.list, name = "BCICR_m", simplify = FALSE)
  return(list(data=data, bcicr = data_bcicr, bcicr_m = data_bcicr_m))
}

oddsratio <- function(control,treatment=NULL,or=NULL){
  if (is.null(or)){
    return(treatment*(1-control)/(control*(1-treatment)))
  } else {
    return(or*control/((or-1)*control + 1))
  }
}


## line plot function
linePlot <- function(data, cutoff_grp = c(0.4,0.6,0.8,0.95), filename = NULL,suffix = "",
                     x="delta", y = "LCP",device = "png", type = "continuous", 
                     variance=NULL, effect=NULL,var_val = 2, effect_val = 1, treatment=NULL, control=NULL,expect_hr = NULL,enr = NULL,median_tte = NULL, censor_tte = NULL, 
                     height = 12, width = 12){
  if (type == "continuous"){
    data_plot <- data[data$cutoff %in% cutoff_grp & data$variance == variance & data$effect == effect & data$variable_type==type,]
  } else if (type == "binary"){
    data_plot <- data[data$cutoff %in% cutoff_grp & data$treatment_rate == treatment & data$control_rate == control & data$variable_type==type,]
    
  } else if (type == "survival"){
    data_plot <- data[data$cutoff %in% cutoff_grp & data$treatment_rate == treatment & data$control_rate == control & data$variable_type==type,]
  }
  if (is.null(filename)){
    filename <- sprintf("figs/%s %s %s-eff %0.2f %s-var %0.1f %s.%s",y, x,effect, effect_val,variance,var_val,suffix,device)
  }
  if (x == "delta"){
    p <- ggplot(mapping = aes_string(x = x, y = y,group = "cutoff"), data = data_plot) +
      geom_line(aes_string(color = "cutoff"),size = 0.5) +
      geom_point(size = 1)+
      facet_grid(prop_local~.)+
      scale_color_brewer(type = "div",name = "Cut-off") +
      scale_x_continuous(limits = c(0,1)) +
      scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0,1))+
      xlab("Delta")+
      theme_gray(base_size = 10)
  } else {
    p <- ggplot(mapping = aes_string(x = x, y = y,group = "cutoff"), data = data_plot) +
      geom_line(aes_string(color = "cutoff"),size = 0.5) +
      geom_point(size = 1)+
      facet_grid(delta~.)+
      scale_color_brewer(type = "div",name = "Cut-off") +
      scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0,1))+
      xlab("Local Proportion")+
      theme_gray(base_size = 10)
  }
  label <- sprintf("Parameter: %s, %s Effect = %0.2f and %s Variance = %0.1f", type, effect,effect_val ,variance , var_val)
  title <- ggdraw() + draw_label(label,fontface = "bold",size = 10)
  p <- plot_grid(title, p,ncol=1,rel_heights = c(0.04,1))
  ggplot2::ggsave(filename, plot = p, device = device, width = width, height = height,units = "cm")
  print(p)
  return(data_plot)
}
# ROC curve: 1-LPP of delta=0 as x-axis and LPP  of delta > 0 as y-axis
# 