## Extended analysis of antiviral effectiveness and proportion undergoing prompt Tx ----

xfun::pkg_attach(c("deSolve", "tidyverse", "lubridate", "here", "janitor", "readxl", "lhs", "sensitivity"))

filepath <- ""
nsim <- 100    # dummy var so 'code/0_params.R' runs

scenario_tab <- read.csv(here(filepath,'inputs/scenarios.csv'), row.names = 1)
epi_scenarios <- names(scenario_tab)

# only need baseline scenarios for different season severities
scenario_base <- str_subset(epi_scenarios, "baseline")
scenario_base_short <- unique(str_split_i(scenario_base, pattern = "_", i = 3))

# number of sub-stages per day for each infectious compartment
nstage <- 1

source(here(filepath, 'code/0_functions.R'))
source(here(filepath, 'code/0_params.R'))


# Fix previously sampled parameters at mean values -------------

pars_fixed <- 
  data.frame(vals = (maxvals + minvals)/2,
             pars = c(paste0(parbyage, 1:m), parnames[parnames != parbyage])) %>% spread(pars, vals) %>% 
  pivot_longer(cols = contains(parbyage), values_to = parbyage) %>%
  mutate(agegrp = as.numeric(str_match(name, pattern = "[0-9]"))) %>% select(-name) %>% 
  mutate(p_prescribed1 = p_seekcare * p_d1 * p_prescribed,
         p_prescribed2 = p_seekcare * p_d2 * p_prescribed,
         
         p_a1 = p_prescribed1 * p_adherence,
         p_a2 = p_prescribed2 * p_adherence
  )

# adjust p_a2 for ODE input (need % AFTER d1 care-seekers taken out)
pars_fixed <- pars_fixed %>% mutate(p_a2_adj = p_a2/(1 - p_a1) )


## Values for Av effects against transmission vs hospitalization vs proportion on tx (p_av) -------

# resolution
res <- 0.025

# vary av effectiveness & uptake parameters
av_trans_vec <- seq(0, 1, by = res)
av_hosp_vec  <- seq(0, 1, by = res)
p_av_vec     <- seq(0, 1, by = res)

total_sims <- length(av_trans_vec) * length(av_hosp_vec) + length(av_trans_vec) * length(av_hosp_vec) + length(av_trans_vec) * length(p_av_vec)

# Combine into a data.frame of all desired combinations
pars_sim <- 
  expand.grid(par1 = av_trans_vec,
             par2 = p_av_vec) %>%
  mutate(combo = "av_trans x p_av") %>%
  # repeat for sensitivity analysis with av_hosp = 0.9 (equivalent to a 10% reduction in risk)
  rbind(
    expand.grid(par1 = av_trans_vec,
                par2 = p_av_vec) %>%
      mutate(combo = "av_trans x p_av sensitivity")
  ) %>%
  rbind(
    expand.grid(par1 = av_trans_vec,
               par2 = av_hosp_vec) %>%
      mutate(combo = "av_trans x av_hosp")
    ) 


## Run model ------------------

times <- seq(from = 1, to = 365) 

out_spectrum <- list()

# run different severity scenarios
for (e in 1:length(scenario_base)) {
  epi_scenario <- scenario_base[e] 
  print(epi_scenario)
  
  # get scenario parameters (now just Re)
  Reff <- scenario_tab['Re', epi_scenario] 
  
  # adjust symptomatic proportions for vaccination
  p_sym_corr  <- p_sym * (1 - vc * ve)
  
  # approximate beta (fixed at baseline value without antivirals, so redundant sometimes)
  beta_vals <- calculate_beta(Reff = Reff, 
                              contacts = contact_mat,
                              p_immune = prop_immune,
                              asym_adj = asym_adj,
                              av_trans_adj = 1)  # no difference in transmission
  
  # adjust beta for ODE input
  beta0 <- beta_vals * contact_mat 
  
  # simulate ODEs
  tmp_out <- list()

  for(j in 1:nrow(pars_sim)){ 
      
    # determine which av effectiveness / uptake parameters are fixed for each sim
    
    if (pars_sim$combo[j] == "av_trans x p_av"){
      
      p_a1 <- p_a2 <- rep(pars_sim$par2[j]/2, m)
      p_a2_adj <- p_a2/(1 - p_a1)
      av_hosp  <- 0.8 # fixed hospitalization risk reduction of 20%
      
    } else if (pars_sim$combo[j] == "av_trans x p_av sensitivity"){
      
      p_a1 <- p_a2 <- rep(pars_sim$par2[j]/2, m)
      p_a2_adj <- p_a2/(1 - p_a1)
      av_hosp  <- 0.9 # fixed hospitalization risk reduction of 10%
      
    } else {
      p_a1 <- pars_fixed$p_a1  # fixed % on tx after 1d
      p_a2 <- pars_fixed$p_a2  # fixed % on tx after 2d
      
      p_a2_adj <- pars_fixed$p_a2_adj
      av_hosp  <- pars_sim$par2[j]
    }
    
     # see where we are in computation
     if (ceiling(j/500) == j/500) print(paste0("sim ", j, "/", total_sims))
      
      tmp_params <- 
        list(m = m,                                      # number of age groups
             beta0 = as.matrix(beta0),                   # transmission (no AVs)
             adj_asym = asym_adj,                        # asymptomatic transmission adjustment
             adj_avi  = pars_sim$par1[j],                # AV transmission adjustment
             
             p_sym  = p_sym_corr,                        # proportion asymptomatic
             p_a1   = p_a1,                              # proportion starting AVs after 1d
             p_a2   = p_a2,                              # proportion starting AVs after 2d
             
             sigma   = 1/ dur_latent,                    # rate E --> A or Ip
             gamma_p = 1 / dur_presym_inf,               # rate Ip --> Is
             gamma_a = 1 / dur_asym,                     # rate A --> R
             gamma_s = nstage * (1 / 1),                 # rate Is_i --> Is_i+1 or R
             gamma_a1 = nstage * (1 / 1),                # rate Ia1_i --> Ia1_i+1 or R
             gamma_a2 = nstage * (1 / 1),                # rate Ia2  --> R
             
             p_h = hosp_risk,                            # risk of hospitalization
             r_h = av_hosp                               # antiviral hospitalization risk adjustment
        )
      
      tmp_out[[j]] <- as.data.frame(ode(y=state, times=times, func=flu.model.1stage, parms=tmp_params, method="lsodes")) %>% 
                          mutate(scenario = epi_scenario,
                                 av_trans = pars_sim$par1[j],
                                 av_hosp = av_hosp, 
                                 p_av = mean(p_a1 + p_a2),
                                 combination = pars_sim$combo[j]) %>% 
                          # only save what's necessary
                          filter(time == max(time)) %>%
                          combine_out() %>%
                          select(scenario, av_trans, av_hosp, p_av, combination, 
                                 H_all, count_Ia1_all, count_Ia2_all)

    }
  out_spectrum[[e]] <- bind_rows(tmp_out) 
}

out_spectrum <- bind_rows(out_spectrum) %>% separate(scenario, c("X", "severity", "scenario"), sep = "_") 

save(out_spectrum, file = here(filepath, paste0("out/spectrum_", res, ".RData")))



