## functions
## Non-parametric bootstrap
## simulation
library(dplyr)
library(purrr)
library(survival)
library(ggplot2)
## generate parameter setting 
paramGenerate <- function(num_pat, num_region,prop_region,num_arm = 2,ratio = NULL, num_param = NULL,name_param = NULL, par1_arm1 =NULL, par1_arm2=NULL, more_par=NULL){
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
  # shape parameter for the event
  # scale parameter for the event
  # shape parameter for the censoring
  # scale parameter for the censoring
  # maximum length of follow-up
  num <- par$Num[1]
  par_event <- par$Value[2:3]
  par_censor <- par$Value[4:5]
  par_enr <- par$Value[1]
  par_fol <- par$Value[6]
  enroll_time <- seq(1, num/par_enr, length.out = num)
  event_free_time <- rweibull(num, par_event[1], par_event[2]) + enroll_time
  censor_time <- rweibull(num, par_censor[1], par_censor[2]) + enroll_time
  EOS <- pmin(event_free_time, censor_time, par_fol)
  event <- as.numeric(EOS == event_free_time)
  return(data.frame(start_time = enroll_time, end_time = EOS, event = event))
}

simulate <- function(params, sim_func){
  ## the main function to simulate dataset
  # params: the parameter generated
  # sim_func: the internal simulation function to use
  value <- params %>% split(.,list(.$Region,.$Arm)) %>% map_dfr(sim_func)
  dat <- distinct(params[c("Region","Arm","Num")])
  dat <- params[rep(1:nrow(dat),dat$Num),c("Region","Arm")]
  dat <- bind_cols(dat, value)
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
  if (statistics == "mean" | statistics == "proportion"){
    if (method == "np"){
      local_weight <- bootWeighting(num_boot,data = dat_local, sample_size = nrow(dat_local),probability = probability)
      local_diff <- bootDiffMean(dat_local, weight = local_weight,control = control, treatment = treatment)
      local_CI <- quantile(local_diff,probs = c(0.025, 0.975),na.rm = TRUE)
      nonlocal_weight <- bootWeighting(num_boot,data = dat_nonlocal, sample_size = nrow(dat_local),probability = probability)
      nonlocal_diff <- bootDiffMean(dat_nonlocal, weight = nonlocal_weight,control = control, treatment = treatment)
      nonlocal_CI <- quantile(nonlocal_diff,probs = c(0.025, 0.975),na.rm=TRUE)
    }
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
  params$Num <- params$Num * num_sim
  dat <- simulate(params = params,sim_func = sim_func)
  bcicr_result <- map(1:num_sim, ~BCICR(dat = dat[index + (.-1)*rep(num_per_group,num_per_group),], local_index = local,treatment = treatment, control = control,num_boot = num_boot,statistics = statistics,probability = probability))
  bcicr <- getListElement(bcicr_result, "BCICR")
  bcicr_m <- getListElement(bcicr_result, "BCICR_m")
  LPP <- map_dbl(cutoff, ~LPP(bcicr, .))
  LPP_m <- map_dbl(cutoff, ~LPP(bcicr_m,.))
  return(list(BCICR = bcicr,BCICR_m = bcicr_m,cutoff=cutoff,LPP = LPP, LPP_m = LPP_m))
}


#sim scenarios
scenarios_sim <- function(prop_local, sample_size, delta,num_sim = 1000, num_boot = 1000,
                          type = "continuous", mean_diff = 1,mean_variance = 2,unequal_var = FALSE, unequal_eff = FALSE, NL_var = NULL,NL_delta = NULL,treatment_rate = 0.3,control_rate = 0.1,seed = 20190117){
  ## perform simulation scenarios for different scenarios
  # prop_local: a vector of local population proportions
  # sample_size: a vector of the number of patients, if sample_size  = NULL, it will be calculated by power.
  # delta: the difference of effect size in percentage 0 to 1
  # unequal_var: the variance between multiple non-locol population are not equal
  # unequal_effect: the effects in non-local population are not equal
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
      mean_var <-  c(prop_region %*% c(NL_var,mean_variance))**2*2
      sample_size <- ceiling((qnorm(0.025)+qnorm(0.1))**2*2*mean_var/mean_diff**2)
      scenarios$sample_size <- rep(sample_size, length(delta))
    }
    NL_delta <- NL_num*mean_diff/sum(NL_delta) * NL_delta
    params.list <- pmap(scenarios, ~paramGenerate(num_pat=..1,num_region = NL_num+1,num_param = 2,prop_region = c(rep((1 - ..2)/NL_num,NL_num),..2),
                                                  par1_arm1 = c(NL_delta,0+mean_diff*..3),par1_arm2 = rep(0,NL_num+1), more_par = list(par2_arm1 = c(NL_var,mean_variance), par2_arm2 = c(NL_var,mean_variance))))
    sim <- map(params.list, ~BCICR_sim(params = .,
                                       cutoff = seq(0,1,0.01),
                                       local =NL_num+1,
                                       num_sim = num_sim,
                                       num_boot = num_boot,
                                       sim_func = sim.continuous,
                                       seed = seed))
    scenarios$sample_size <- map_dbl(params.list, ~sum(.$Num)/length(unique(.$Param)))
    scenarios$variable_type <- type
    scenarios$unequal_variance <- unequal_var
    scenarios$unequal_effect <- unequal_eff
  } else if (type == "binary"){
    if (sample_size == 0){
      mean_diff = treatment_rate-control_rate
      sample_size <- ceiling((qnorm(0.025)+qnorm(0.1))**2*(treatment_rate*(1 - treatment_rate)+control_rate*(1-control_rate))/mean_diff**2)
      scenarios$sample_size <- sample_size
    }
    params.list <- pmap(scenarios, ~paramGenerate(num_pat=..1,num_region =2,num_param = 1,prop_region = c(1 - ..2,..2),
                                                  par1_arm1 = c(treatment_rate,(treatment_rate -control_rate)*..3+control_rate),
                                                  par1_arm2 = rep(control_rate,2)))
    sim <- map(params.list, ~BCICR_sim(params = .,
                                       cutoff = seq(0,1,0.01),
                                       local =2,
                                       num_sim = num_sim,
                                       num_boot = num_boot,
                                       sim_func = sim.binary,
                                       seed = seed))
    scenarios$sample_size <- map_dbl(params.list, ~sum(.$Num))
    scenarios$variable_type <- type
    scenarios$treatment_rate <- treatment_rate
    scenarios$control_rate <- control_rate
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
# ROC curve: 1-LPP of delta=0 as x-axis and LPP  of delta > 0 as y-axis
# 