params <- paramGenerate(num_pat = 400, num_region = 4, prop_region = c(0.35,0.30,0.30,0.05),
                        num_arm = 2,
                        num_param = 1,name_param = c("ORR"),
                        par1_arm1 = c(0.25,0.25,0.25,0.15),par1_arm2 = c(0.1,0.1,0.1,0.1))

result4 <- scenarios_sim(prop_local = c(0.5),sample_size = 400,delta = c(0),num_sim = 100,num_boot = 1000,type = "binary",treatment_rate = 0.4,control_rate = 0.1)

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
