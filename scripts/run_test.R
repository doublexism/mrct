## test run

# 1. continuous
result1t <- scenarios_sim(prop_local = c(0.05, 0.10,0.15,0.20,0.25), # proportion of local population
                         sample_size = 0,                    # sample size of the trial
                         delta = c(0,0.25,0.5,0.75,1),       # ratio of local effect to non-local effect
                         num_sim = 1000,                       # number of simulated datasets
                         num_boot = 1000,                      # number of bootstrap resampling
                         type = "continuous",                  # continuous enpoints
                         mean_diff = 0.5,
                         mean_variance = 2,
                         unequal_var = FALSE,                  # equal variance between non-local population
                         unequal_eff = FALSE)                # specifing the effect size of the 4 non-local populations


result2t <- scenarios_sim(prop_local = c(0.05, 0.10,0.15,0.20,0.25), # proportion of local population
                          sample_size = 0,                    # sample size of the trial
                          delta = c(0,0.25,0.5,0.75,1),       # ratio of local effect to non-local effect
                          num_sim = 100,                       # number of simulated datasets
                          num_boot = 100,                      # number of bootstrap resampling
                          type = "continuous",                  # continuous enpoints
                          mean_diff = 0.75,
                          mean_variance = 2,
                          unequal_var = FALSE,                  # equal variance between non-local population
                          unequal_eff = FALSE)                # specifing the effect size of the 4 non-local populations

result3t <- scenarios_sim(prop_local = c(0.05, 0.10,0.15,0.20,0.25), # proportion of local population
                          sample_size = 0,                    # sample size of the trial
                          delta = c(0,0.25,0.5,0.75,1),       # ratio of local effect to non-local effect
                          num_sim = 1000,                       # number of simulated datasets
                          num_boot = 1000,                      # number of bootstrap resampling
                          type = "continuous",                  # continuous enpoints
                          mean_diff = 1,
                          mean_variance = 2,
                          unequal_var = FALSE,                  # equal variance between non-local population
                          unequal_eff = FALSE)                # specifing the effect size of the 4 non-local populations

result3t <- scenarios_sim(prop_local = c(0.05, 0.10,0.15,0.20,0.25), # proportion of local population
                          sample_size = 0,                    # sample size of the trial
                          delta = c(0,0.25,0.5,0.75,1),       # ratio of local effect to non-local effect
                          num_sim = 1000,                       # number of simulated datasets
                          num_boot = 1000,                      # number of bootstrap resampling
                          type = "continuous",                  # continuous enpoints
                          mean_diff = 1,
                          mean_variance = 4,
                          unequal_var = FALSE,                  # equal variance between non-local population
                          unequal_eff = FALSE)                # specifing the effect size of the 4 non-local populations
result <- simResult(result1t, result2t, result3t)
saveRDS(result, sprintf("data/result_sum_%s.RDS",Sys.Date()))
