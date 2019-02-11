# 
params <- paramGenerate(num_pat = 600, 
                        num_region = 2, 
                        prop_region = c(0.75,0.25),
                        num_arm = 2,
                        num_param = 6,
                        name_param = c("enr_rate","a_event","b_event","a_censor","b_censor","follow-up time"),
                        par1_arm1 = round(50*c(0.75,0.25)/2),
                        par1_arm2 = round(50*c(0.75,0.25)/2),
                        more_par = list(
                        par2_arm1 = c(1,1),
                        par2_arm2 = c(1,1),
                        par3_arm1 = c(10, 10),
                        par3_arm2 = c(14.3/log(2), 13.3/log(2)),
                        par4_arm1 = c(1, 1),
                        par4_arm2 = c(1, 1),
                        par5_arm1 = c(100,100),
                        par5_arm2 = c(100,100),
                        par6_arm1 = c(60,60),
                        par6_arm1 = c(60,60)))
dat <- simulate(params, sim_func = sim.survival)                        
coxph(formula = Surv(end_time - start_time, event)~Arm,data = dat) %>% summary()
