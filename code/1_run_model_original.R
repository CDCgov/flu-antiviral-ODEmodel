## Set up ------------

xfun::pkg_attach(c("deSolve", "tidyverse", "here", "lhs"))

filepath <- ""

set.seed(99)

# number of simulations
nsim <- 100    

# number of sub-stages per day for each infectious compartment
nstage <- 1

source(here(filepath, 'code/0_functions.R'))
source(here(filepath, 'code/0_params.R'))

# scenario information
scenario_tab <- read.csv(here(filepath,'inputs/scenarios.csv'), row.names = 1)
epi_scenarios <- names(scenario_tab)

scenario_short <- unique(str_split_i(epi_scenarios, pattern = "_", i = 3))
scenario_labels <- c("Baseline", "No transmission impact", "Small transmission impact", "Large transmission impact")

# simulation times
times <- seq(from = 1, to = 365) 


## Simulate model for different scenarios ---------------

out <- list()

for (e in 1:length(epi_scenarios)) {
  epi_scenario <- epi_scenarios[e] 
  print(epi_scenario)
  
  # get scenario parameters
  av_trans_adj         <- scenario_tab['av_transmissibility_adj', epi_scenario]      # effect of av on transmissibility
  av_hosp_risk_multiplier  <- scenario_tab['av_hosp_risk_multiplier', epi_scenario]  # effect of av on CHR
  Reff                     <- scenario_tab['Re', epi_scenario]                       # Re for low vs high severity
  
  # adjust symptomatic proportions for vaccination
  p_sym_corr  <- p_sym * (1 - vc * ve)
  
  # approximate beta (fixed at baseline value without antivirals, so redundant sometimes)
  beta_vals <- calculate_beta(Reff = Reff, 
                              contacts = contact_mat,
                              p_immune = prop_immune,
                              asym_adj = asym_adj,
                              av_trans_adj = 1) # no difference in transmission
  
  # adjust beta for ODE input
  beta0 <- beta_vals * contact_mat 

  # adjust p_a2 for ODE input (need % AFTER d1 proportion taken out)
  pars_sim_adj <- pars_sim %>% mutate(p_a2_adj = p_a2/(1 - p_a1) )
  
  # run simulations
  tmp_out <- list()

  for (n in 1:nsim) {
    tmp_params <- 
      list(m = m,                                         # number of age groups
           beta0 = as.matrix(beta0),                      # transmission
           adj_asym = asym_adj,                           # asymptomatic transmission adjustment
           adj_avi = av_trans_adj,                        # AV transmission adjustment
           p_sym  = p_sym_corr,                           # proportion symptomatic
           p_a1   = subset(pars_sim_adj, sim == n)$p_a1,      # proportion starting AVs after 1d
           p_a2   = subset(pars_sim_adj, sim == n)$p_a2_adj,  # proportion starting AVs after 2d
           sigma = 1/ dur_latent,                         # rate E --> A or Ip
           gamma_p = 1 / dur_presym_inf,                  # rate Ip --> Is
           gamma_a = 1 / dur_asym,                        # rate A --> R
           gamma_s = nstage * (1 / 1),                    # rate Is_i --> Is_i+1 or Iai
           gamma_a1 = nstage * (1 / 1),                   # rate Ia1_i --> Ia1_i+1 or R
           gamma_a2 = nstage * (1 / 1),                   # rate Ia2  --> R
           p_h = hosp_risk,                               # risk of hospitalization
           r_h = av_hosp_risk_multiplier                  # antiviral hospitalization risk adjustment
      )
    
    tmp_out[[n]] <- as.data.frame(ode(y=state, times=times, func=flu.model.1stage, parms=tmp_params, method="lsodes")) %>% 
                      mutate(scenario = epi_scenario, sim = n)
  }
  out[[e]] <- bind_rows(tmp_out) 
}

out <- bind_rows(out)

save(out, file = here(filepath, paste0("out/original_", nsim, ".RData")) )


## End --------------

