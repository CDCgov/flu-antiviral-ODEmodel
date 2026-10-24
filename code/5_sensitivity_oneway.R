## Set up ------------

xfun::pkg_attach(c("deSolve", "tidyverse", "here", "lhs",  "scico"))

filepath <- ""
figpath <- "figs"

set.seed(99)
nsim <- 100

# number of stages
nstage <- 1

source(here(filepath, 'code/0_functions.R'))
source(here(filepath, 'code/0_params.R'))

# scenario information
scenario_tab <- read.csv(here(filepath,'inputs/scenarios.csv'), row.names = 1)

# only need one set of antiviral effectiveness scenarios
# since varying Re independently and calculating baseline within each sim
epi_scenarios <- names(scenario_tab)[2:4] 

scenario_short <- unique(str_split_i(epi_scenarios, pattern = "_", i = 3))
scenario_labels <- c("No impact", "Small impact", "Large impact")

# simulation times
times <- seq(from = 1, to = 365) 

# Fix % undergoing prompt Tx at mean values for each age group
# (and adjust p_a2 for ODE input (need % AFTER d1 proportion taken out))
pars_sim_mean <- pars_sim %>% 
                  mutate(p_a2_adj = p_a2/(1 - p_a1) ) %>%
                  group_by(agegrp) %>% 
                  summarize(p_a1 = mean(p_a1),
                            p_a2_adj = mean(p_a2_adj))


## Vary select parameters ---------------------

# parameters to include in analysis
parnames <- c("p_asym", "asym_adj", "Re", "hosp_risk", "prop_immune")

# number of parameters
npar_sim <- length(parnames) 

# define parameter ranges
minvals <- c(0.2, 0.3, 1.10, 1/400, 0.05)
maxvals <- c(0.8, 0.9, 1.25, 1/20,  0.55)

# reference values
baseline <- apply(rbind(minvals, maxvals), 2, median)

names(minvals) <- names(maxvals) <- names(baseline) <- parnames

# number of levels to test within each range
n_levels <- 11

# list of all parameter combinations
pars_levels <- 
  lapply(parnames, function(p){
  varied_vals <- seq(minvals[p], maxvals[p], length.out = n_levels)

  block <- as.data.frame(matrix(baseline, nrow = n_levels, ncol = npar_sim,
                                byrow = TRUE, dimnames = list(NULL, parnames)))
  block[[p]] <- varied_vals
  block$varied_par <- p             
  block$level      <- varied_vals 
  block
})

pars_levels <- bind_rows(pars_levels) %>% mutate(sim = 1:nrow(.))

niter <- nrow(pars_levels)


## Simulate model ---------------

tmp_out <- list()

k <- 1

for (n in 1:niter) {
    
  pars_tmp <- pars_levels[n,]
  
  # reset initial conditions based on sampled proportion immune
  state[paste0("S", 1:m)] <- N_byage - 1000 - (N_byage) * pars_tmp$prop_immune
  state[paste0("R", 1:m)] <- (N_byage) * pars_tmp$prop_immune
  
  Reff <- pars_tmp$Re   
  
  # Don't correct symptomatic proportion by VC/ VE (since varying anyway; easier interpretation of values)
  p_sym_corr <- 1-rep(pars_tmp$p_asym, m)

  # calc beta 
  beta_vals <- calculate_beta(Reff = Reff, 
                              contacts = contact_mat,
                              p_immune = pars_tmp$prop_immune,
                              asym_adj = pars_tmp$asym_adj,
                              av_trans_adj = 1) # no difference in transmission
 
  # Adjust beta 
  beta0 <- beta_vals * contact_mat 
  
  for (e in 1:length(epi_scenarios)) {
  
    epi_scenario <- epi_scenarios[e] 
  
    # get scenario parameters
    av_trans_adj         <- scenario_tab['av_transmissibility_adj', epi_scenario]      # effect of av on transmissibility
    av_hosp_risk_multiplier  <- scenario_tab['av_hosp_risk_multiplier', epi_scenario]  # effect of av on CHR

    # ODE input parameters
    tmp_params <- 
      list(m = m,                                     # number age groups
          beta0 = as.matrix(beta0),                   # transmission (no AVs)
          adj_asym = pars_tmp$asym_adj,               # asymptomatic transmission adjustment
          adj_avi  = av_trans_adj,                    # AV transmission adjustment
          
          p_sym = p_sym_corr,                         # proportion asymptomatic
          p_a1   = pars_sim_mean$p_a1,                # proportion AVs ~ 1d
          p_a2   = pars_sim_mean$p_a2_adj,            # proportion AVs ~ 2d
          sigma = 1/ dur_latent,                         # rate E --> A or Ip
          gamma_p = 1 / dur_presym_inf,                  # rate Ip --> Is
          gamma_a = 1 / dur_asym,                        # rate A --> R
          gamma_s = nstage * (1 / 1),                    # rate Is_i --> Is_i+1 or Iai
          gamma_a1 = nstage * (1 / 1),                   # rate Ia1_i --> Ia1_i+1 or R
          gamma_a2 = nstage * (1 / 1),                   # rate Ia2  --> R
          p_h = rep(pars_tmp$hosp_risk, m),              # risk of hospitalization
          r_h = av_hosp_risk_multiplier                  # antiviral hospitalization risk adjustment
      )
    
    # with antivirals
    tmp <- as.data.frame(ode(y=state, times=times, func=flu.model.1stage, parms=tmp_params, method="lsodes")) 
    
    # without antivirals
    tmp_params$r_h <- tmp_params$adj_avi <- 1
    
    tmp0 <- as.data.frame(ode(y=state, times=times, func=flu.model.1stage, parms=tmp_params, method="lsodes")) 
    
    
    # Derive # prevented hospitalizations
    tmp_out[[k]] <- tmp %>% filter(time == max(time)) %>%
                      rbind(tmp0 %>% filter(time == max(time))) %>% 
                      combine_out() %>%
                      select(H_all) %>%
                      summarize(out_abs = diff(H_all),
                                out_perc = (H_all[2]-H_all[1])/H_all[2] * 100) %>%
                      mutate(sim = n, scenario = epi_scenario)
  
    k <- k+1
  }
}

out_levels <- bind_rows(tmp_out)

input_levels <- pars_levels %>% left_join(out_levels) %>% 
                  gather(metric, out, contains("out")) %>%
                  separate(scenario, c("X", "severity", "scenario"), sep = "_") %>%
                  mutate(scenario_names = rename_scenarios(scenario, scenario_short, scenario_labels),
                         metric = ifelse(metric == "out_perc", "Percent", "Absolute"))


## Tornado plot showing range of outputs ----------

scenario_cols <- scico(7, palette = 'lapaz')[c(3, 4, 6)]

# Store baseline values
baseline_out <- input_levels %>%
  group_by(varied_par, scenario_names, metric) %>%
  slice(median(1:n_levels)) %>%   
  group_by(scenario_names, metric)  %>%
  distinct(out = round(out, 4))

baseline_out <- baseline_out %>% tidyr::expand_grid(varied_par = rename_pars(parnames))

# Get min/max outputs for each varied parameter
tornado <- input_levels %>%
  group_by(varied_par, scenario_names, metric) %>%
  summarise(ymin = out[which.min(level)],
            ymax = out[which.max(level)],
          .groups = "keep") %>%
  mutate(width = abs(ymax - ymin)) %>%
  arrange(scenario_names, metric, width) %>%
  group_by(scenario_names, metric) %>%
  mutate(varied_par = rename_pars(varied_par)) %>%
  mutate(varied_par = factor(varied_par, levels = varied_par))

tornado <- input_levels %>%
  group_by(varied_par, scenario_names, metric) %>%
  summarise(ymin = min(out),
            ymax = max(out),
          .groups = "keep") %>%
  mutate(width = abs(ymax - ymin)) %>%
  arrange(scenario_names, metric, width) %>%
  group_by(scenario_names, metric) %>%
  mutate(varied_par = rename_pars(varied_par)) %>%
  mutate(varied_par = factor(varied_par, levels = varied_par))

# Create the plot
tornado %>% 
  ggplot(aes(y = varied_par, color = scenario_names, group = scenario_names)) +
  geom_linerange(aes(xmin = ymin, xmax = ymax), linewidth = 5, position = position_dodge(width = 0.6)) +
  geom_point(data = baseline_out, aes(x = out), shape = 73, size = 4.5, position = position_dodge(width = 0.6)) +
  geom_point(data = baseline_out, aes(x = out), color= "black", position = position_dodge(width = 0.6)) +
  facet_wrap(~ metric, scales = "free_x") +
  scale_color_manual(values = scenario_cols, breaks = rev(scenario_labels)) + 
  labs(y = NULL, x = "Prevented hospitalizations", color = "Antiviral effect\non transmission") +
  scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  scale_y_discrete(labels = label_wrap_gen(15)) +
  get_theme(panel.spacing.x = unit(1.2, "lines"))
ggsave(file = here(figpath, "sensitivity_oneway.pdf"), width = 9, height = 6)


## End ------------

