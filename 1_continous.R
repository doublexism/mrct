result1 <- scenarios_sim(prop_local = c(0.05, 0.10,0.15,0.20), # proportion of local population
                         sample_size = 400,                    # sample size of the trial
                         delta = c(0,0.2,0.4,0.6,0.8,1),       # ratio of local effect to non-local effect
                         num_sim = 1000,                       # number of simulated datasets
                         num_boot = 1000,                      # number of bootstrap resampling
                         type = "continuous",                  # continuous enpoints
                         unequal_var = FALSE,                  # equal variance between non-local population
                         unequal_eff = TRUE,                   # equal effect size between non-local population
                         NL_delta = c(1,0,1,0))                # specifing the effect size of the 4 non-local populations
result2 <- scenarios_sim(prop_local = c(0.05, 0.10,0.15,0.20),
                         sample_size = 400,
                         delta = c(0,0.2,0.4,0.6,0.8,1),
                         num_sim = 1000,
                         num_boot = 1000,
                         type = "continuous",
                         unequal_var = TRUE,
                         unequal_eff = FALSE,
                         NL_var = c(1,3,3,1))
result3 <- scenarios_sim(prop_local = c(0.05, 0.10,0.15,0.20),
                         sample_size = 400,
                         delta = c(0,0.2,0.4,0.6,0.8,1),
                         num_sim = 1000,
                         num_boot = 1000,
                         type = "continuous",
                         unequal_var = TRUE,
                         unequal_eff = TRUE,
                         NL_var = c(1,3,3,1),
                         NL_delta = c(1,0,1,0))

## merger results of the 3 scenarios
result_sum <- simResult(result1, result2,result3) 

