resultb1 <- scenarios_sim(prop_local = c(0.05, 0.10,0.15,0.20), # proportion of local population
                         sample_size = 0,                    # sample size of the trial
                         delta = c(0.5, 1),       # ratio of local effect to non-local effect
                         num_sim = 100,                       # number of simulated datasets
                         num_boot = 100,                      # number of bootstrap resampling
                         type = "binary",
                         treatment_rate = 0.40,
                         control_rate = 0.25)                # specifing the effect size of the 4 non-local populations

dat <- simulate(params = params,sim_func = sim.binary)
summary(dat, "mean")

weight <- bootWeighting(1000, dat[dat$Region !=4,], 20)
diff <- bootDiffMean(dat[dat$Region !=4,], weight  = weight, control = 2, treatment = 1)

w <- weight[,1:25]
v <- dat[dat$Region != 4,]$value
i <- 2*(dat[dat$Region != 4,]$Arm ==1)-1

w[i == 1,] <- w[i == 1,]%*% diag(1/colSums(w[i==1,]))
w[i == -1,] <- w[i == -1,] %*% diag(1/colSums(w[i==-1,]))
w <- w * i

boot_diff <- c(v %*% w)
