## Set up ------------

xfun::pkg_attach(c("deSolve", "tidyverse", "here", "lhs",  "sensitivity"))

filepath <- ""

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
scenario_labels <- c("None", "Small", "Large")

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

# parameters to include in LHS sampling
parnames <- c("p_asym", "asym_adj", "Re", "hosp_risk", "prop_immune")

# number of parameters
npar_sim <- length(parnames) 

# define parameter ranges
minvals <- c(0.2, 0.2, 1.05, 1/400, 0.1)
maxvals <- c(0.8, 0.8, 1.30, 1/50,  0.5)

# draw LHS sample
parsamples <- maximinLHS(nsim, npar_sim)

# transform LHS samples to parameter space
parsamples <- t(t(parsamples) * (maxvals - minvals) + minvals)

colnames(parsamples) <- parnames

pars_prcc <- data.frame(parsamples) %>% mutate(sim = 1:nsim) 


## Simulate model ---------------

tmp_out <- list()

k <- 1

for (n in 1:nsim) {
    
  if(n/10 == round(n/10)) print(paste0("simulation ", n, "/", nsim))
  
  pars_tmp <- subset(pars_prcc, sim == n)
  
  # reset initial conditions based on sampled proportion immune
  state[paste0("S", 1:m)] <- N_byage - 1000 - (N_byage) * pars_tmp$prop_immune
  state[paste0("R", 1:m)] <- (N_byage) * pars_tmp$prop_immune
  
  Reff  <- pars_tmp$Re   
  
  p_sym <- 1-rep(pars_tmp$p_asym, m)

  # Don't correct symptomatic proportion by VC/ VE (since varying anyway; easier interpretation of values)
  p_sym_corr  <- p_sym  #* (1 - vc * ve)

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
      list(m = m,                                      # number age groups
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
                      #summarize(H_all = diff(H_all)) %>% 
                      summarize(out_abs = diff(H_all),
                                out_perc = (H_all[2]-H_all[1])/H_all[2] * 100) %>%
                      mutate(sim = n, scenario = epi_scenario)
  
    k <- k+1
  }
}

out_prcc <- bind_rows(tmp_out)

input_prcc <- pars_prcc %>% left_join(out_prcc) %>% 
              gather(metric, out, contains("out")) %>%
              separate(scenario, c("X", "severity", "scenario"), sep = "_") %>%
              mutate(scenario_names = rename_scenarios(scenario, scenario_short, scenario_labels),
                     metric = ifelse(metric == "out_perc", "Percent", "Absolute"))


## Plotting function ----------
rename_pars <- function(parnames){
  parnames[parnames == "Reff"] <- "Re"
  parnames[parnames == "hosp_risk"]  <- "Risk of hospitalization"
  parnames[parnames == "p_asym"] <- "Proportion asymptomatic"
  parnames[parnames == "asym_adj"] <- "Asymptomatic transmission scaling"
  parnames[parnames == "prop_immune"] <- "Proportion immune"
  
  return(factor(parnames, 
                levels = c("Re", 
                           "Risk of hospitalization", 
                           "Proportion immune",
                           "Proportion asymptomatic",
                           "Asymptomatic transmission scaling")))
}

## Monotonicity check ----------

# Ranked outputs
input_prcc %>% 
  gather(par, val, all_of(parnames)) %>%
  group_by(par, metric, scenario) %>%
  mutate(rank_x = rank(val), rank_y = rank(out)) %>%
  mutate(parnames = rename_pars(par)) %>%
  ggplot(aes(x = rank_x, y = rank_y, color = scenario_names)) +
  geom_point(alpha = 0.5) + 
  facet_grid(metric~ parnames, scales = "free", labeller = label_wrap_gen(width = 15)) + 
  geom_smooth(method = "loess", se = FALSE) + 
  get_theme() + labs(x = "Input rank", y = "Output rank") +
  scale_color_brewer("Antiviral\nimpact\non\ntransmission", palette = "Dark2")

# Raw outputs
input_prcc %>% 
  gather(par, val, all_of(parnames)) %>%
    mutate(parnames = rename_pars(par)) %>%
  ggplot(aes(x = val, y = out, color = scenario_names)) +
  geom_point(alpha = 0.5) + 
  facet_grid(metric~ parnames, scales = "free", labeller = label_wrap_gen(width = 15)) +
  geom_smooth(method = "loess", se = FALSE) + 
  get_theme() + scale_y_log10() + labs(x = "Input", y = "Output") +
  scale_color_brewer("Antiviral\nimpact\non\ntransmission", palette = "Dark2")


## Calculate PRC coefficients -------------

prc_total <- input_prcc %>% group_by(scenario_names, metric) %>%
              do(pcc(X = .[parnames], y = .$out, rank = TRUE, nboot = 1000, conf = 0.95)$PRCC) %>% 
              mutate(parnames = parnames)

src_total <- input_prcc %>% group_by(scenario_names, metric) %>%
              do(src(X = .[parnames], y = .$out, rank = FALSE, nboot = 1000, conf = 0.95)$SRC) %>% 
              mutate(parnames = parnames)


figpath <- "figs"

prc_total %>%
  mutate(parnames = rename_pars(parnames)) %>%
  ggplot(aes(x = parnames, y = original, color = scenario_names)) + 
  geom_rect(data = NULL, aes(ymin = -1, ymax = -0.5, xmin = -Inf, xmax = Inf), color = NA, alpha = 0.05, fill = "grey80") +
  geom_rect(data = NULL, aes(ymin = +1, ymax = +0.5, xmin = -Inf, xmax = Inf), color = NA, alpha = 0.05, fill = "grey80") +
  geom_point(position = position_dodge(width = 0.2)) + 
  geom_errorbar(aes(ymin = `min. c.i.`, ymax = `max. c.i.`), position = position_dodge(width = 0.2), width = 0.2) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  get_theme(legend.position = "top") + ylim( c(-1, 1)) +
  labs(x = "Parameter", y = "Partial rank correlation coefficient") +
  scale_x_discrete(labels = label_wrap_gen(15)) + 
  facet_wrap(~ metric, ncol = 1) +
  scale_color_brewer("Antiviral impact on transmission", palette = "Dark2")
#ggsave(file = here(figpath, "prcc_prevented.png"), width = 10, height = 5, dpi = 600)


src_total %>%
  mutate(parnames = rename_pars(parnames)) %>%
  ggplot(aes(x = parnames, y = original, color = scenario_names)) + 
  geom_rect(data = NULL, aes(ymin = -1, ymax = -0.5, xmin = -Inf, xmax = Inf), color = NA, alpha = 0.05, fill = "grey80") +
  geom_rect(data = NULL, aes(ymin = +1, ymax = +0.5, xmin = -Inf, xmax = Inf), color = NA, alpha = 0.05, fill = "grey80") +
  geom_point(position = position_dodge(width = 0.2)) + 
  geom_errorbar(aes(ymin = `min. c.i.`, ymax = `max. c.i.`), position = position_dodge(width = 0.2), width = 0.2) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  get_theme(legend.position = "top") + ylim( c(-1, 1)) +
  labs(x = "Parameter", y = "Standardized regression coefficient") +
  scale_x_discrete(labels = label_wrap_gen(15)) +
  facet_wrap(~ metric, ncol = 1) +
  scale_color_brewer("Antiviral impact on transmission", palette = "Dark2")
#ggsave(file = here(figpath, "src_prevented.png"), width = 10, height = 5, dpi = 600)


## End ------------

