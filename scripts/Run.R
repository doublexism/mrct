# Run simulation
timestamp()
# 1. continuous
result1 <- scenarios_sim(prop_local = c(0.05, 0.10,0.15,0.20,0.25), # proportion of local population
                         parallel = TRUE,
                         sample_size = 0,                    # sample size of the trial
                         delta = c(0,0.25,0.5,0.75,1),       # ratio of local effect to non-local effect
                         num_sim = 1000,                       # number of simulated datasets
                         num_boot = 1000,                      # number of bootstrap resampling
                         type = "continuous",                  # continuous enpoints
                         unequal_var = FALSE,                  # equal variance between non-local population
                         unequal_eff = TRUE,                   # equal effect size between non-local population
                         NL_delta = c(2,0,2,0))                # specifing the effect size of the 4 non-local populations
timestamp()
result2 <- scenarios_sim(prop_local = c(0.05, 0.10,0.15,0.2, 0.25),
                         parallel = TRUE,
                         sample_size = 0,
                         delta = c(0,0.25,0.5,0.75,1),       # ratio of local effect to non-local effect
                         num_sim = 1000,
                         num_boot = 1000,
                         type = "continuous",
                         unequal_var = TRUE,
                         unequal_eff = FALSE,
                         NL_var = c(2,3,3,2))
timestamp()
result3 <- scenarios_sim(prop_local = c(0.05, 0.10,0.15,0.20,0.25),
                         parallel = TRUE,
                         sample_size = 0,
                         delta = c(0,0.25,0.5,0.75,1),       # ratio of local effect to non-local effect
                         num_sim = 1000,
                         num_boot = 1000,
                         type = "continuous",
                         unequal_var = TRUE,
                         unequal_eff = TRUE,
                         NL_var = c(2,3,3,2),
                         NL_delta = c(2,0,2,0))
timestamp()
result4 <- scenarios_sim(prop_local = c(0.05, 0.10,0.15,0.20,0.25), # proportion of local population
                         parallel = TRUE,
                         sample_size = 0,                    # sample size of the trial
                         delta = c(0,0.25,0.5,0.75,1),       # ratio of local effect to non-local effect
                         num_sim = 1000,                       # number of simulated datasets
                         num_boot = 1000,                      # number of bootstrap resampling
                         type = "continuous",                  # continuous enpoints
                         unequal_var = FALSE,                  # equal variance between non-local population
                         unequal_eff = FALSE)                   # equal effect size between non-local population)  
timestamp()
# 2. binary
result5 <- scenarios_sim(prop_local = c(0.05, 0.10,0.15,0.20,0.25), # proportion of local population
                         parallel = TRUE,
                         sample_size = 0,                    # sample size of the trial
                         delta = c(0,0.25,0.5,0.75,1),       # ratio of local effect to non-local effect
                         num_sim = 1000,                       # number of simulated datasets
                         num_boot = 1000,                      # number of bootstrap resampling
                         type = "binary",                  # continuous enpoints
                         treatment_rate =  0.40,                  # equal variance between non-local population
                         control_rate = 0.25)                   # equal effect size between non-local population)
timestamp()
result6 <- scenarios_sim(prop_local = c(0.05, 0.10,0.15,0.20,0.25), # proportion of local population
                         parallel = TRUE,
                         sample_size = 0,                    # sample size of the trial
                         delta = c(0,0.25,0.5,0.75,1),       # ratio of local effect to non-local effect
                         num_sim = 1000,                       # number of simulated datasets
                         num_boot = 1000,                      # number of bootstrap resampling
                         type = "binary",                  # continuous enpoints
                         treatment_rate =  0.18,                  # equal variance between non-local population
                         control_rate = 0.10)                   # equal effect size between non-local population)
timestamp()
# 3. survival -- HR_global = 0.6
result7 <- scenarios_sim(prop_local = c(0.05,0.10,0.15,0.20,0.25), # proportion of local population
                          parallel = TRUE,
                          sample_size = 0,                    # sample size of the trial
                          delta = c(0,0.25,0.5,0.75,1),       # ratio of local effect to non-local effect
                          num_sim = 1000,                       # number of simulated datasets
                          num_boot = 1000,                      # number of bootstrap resampling
                          type = "survival",                  # continuous enpoints
                          expect_hr = 0.7,
                          enr = 30,
                          median_tte_control =  20,                 
                          censor_tte =  100,
                          follow_up = 0)     


result <- simResult(result1,result2,result3,result4,result5,result6, result7)
saveRDS(result, sprintf("data/result_sum_%s.RDS",Sys.Date()))