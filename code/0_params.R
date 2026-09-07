## Pop size and age groups  ------------------

# 2024 population data

age_dt <- data.frame(
  agegrp = 1:5,
  pop = c(18599314, 54533406, 143650318, 62148032, 61179918)
)

N_byage <- age_dt$pop
N <- sum(N_byage)

m <- length(N_byage)


## Contact matrix ----------

# see '0_contact_matrix.R' for source code
load(here(filepath, "inputs/contact_matrix.RData"))

# correct for reciprocity (cijNi = cjiNj) 
totalC       <- contact_mat * N_byage                    
contact_mat  <- (totalC + t(totalC)) / (2 * N_byage) 


## Natural history parameters -------------

asym_adj <- 0.6 

# latent, pre-symptomatic, symptomatic periods
dur_latent = 1.5
dur_presym_inf = 0.5
dur_inf = 3
dur_asym = dur_inf + dur_presym_inf 

# proportion (a)symptomatic
p_asym = c(0.2, 0.45, 0.6, 0.5, 0.45)   
p_sym = 1-p_asym

# risk of hospitalization 
hosp_risk <- 1/c(143.44, 364.71, 178.16, 94.3, 11)                                


## Vaccination coverage/effectiveness --------------

vc <- c(0.57, 0.48, 0.31, 0.42, 0.63)  # 24/25 percentages
ve <- c(0.50, 0.45, 0.40, 0.45, 0.40)  # calibrated 

# Translate reported VC among 6m-4yrs to VC among 0-4yrs
N0 <- 3615598    # 0-6m population

vc_6m4y <- vc[1] * (N_byage[1] - N0)        # total children 6m-4y vaccinated
vc[1]   <- round(vc_6m4y / N_byage[1], 2)   # coverage among all children 0-4y


## Simulated parameters: Care seeking, prescribing, treatment adherence -------

# parameters to include in LHS sampling
parnames <- c("p_seekcare", "p_d1", "p_d2", "p_adherence", "p_prescribed")
# which ones are varied by age
parbyage <- "p_seekcare"

# number of parameters, accounting for age differences
npar_sim <- length(parnames[parnames %in% parbyage]) * m + length(parnames[!(parnames %in% parbyage)])

# draw LHS sample
parsamples <- maximinLHS(nsim, npar_sim)

# define max/min vals for parameters
maxvals <- c(c(0.35, 0.35, 0.50, 0.55, 0.65), # any care-seeking
             0.25, # care seeking within 24hrs, among all care-seekers
             0.25, # care seeking between 24-48hrs, among all care-seekers
             0.85, # adherence
             0.35) # prescribing

minvals <- c(c(0.25, 0.25, 0.35, 0.40, 0.45), # any care-seeking
             0.20, # care seeking within 24hrs, among all care-seekers
             0.20, # care seeking between 24-48hrs, among all care-seekers
             0.75, # adherence
             0.15) # prescribed

# transform LHS samples to parameter space
parsamples <- t(t(parsamples) * (maxvals - minvals) + minvals)

colnames(parsamples) <- c(as.vector(t(outer(parbyage, 1:m, paste0))), parnames[!(parnames %in% parbyage)])


# derive other variables
pars_sim <- 
  data.frame(parsamples) %>% gather(var, val, contains(parbyage)) %>%
  mutate(agegrp = as.numeric(str_match(var, pattern = "[0-9]"))) %>%
  rename(!!parbyage := val) %>% 
  select(-var) %>%
  mutate(sim = rep(1:nsim, times = m)) %>%
  # calculate overall uptake and adherence
  mutate(# care seeking at 1d/2d
         p_care_1d = p_seekcare * p_d1,
         p_care_2d = p_seekcare * p_d2,
         
         # prescribing at 1d/2d
         p_prescribed1 = p_seekcare * p_d1 * p_prescribed,
         p_prescribed2 = p_seekcare * p_d2 * p_prescribed,
         
         # prescribing and adherence at 1d/2d
         p_a1 = p_prescribed1 * p_adherence,
         p_a2 = p_prescribed2 * p_adherence
)


## Initial conditions ----------------

prop_immune <- 0.3 

state <- c(S   = N_byage - 1000 - (N_byage) * prop_immune, 
           E   = rep(0, m),
           A   = rep(500, m), 
           Ip  = rep(500, m), 
           Is_d1  = rep(0, m),
           Is_d2  = rep(0, m),
           Is_d3  = rep(0, m),
           Ia1_d1 = rep(0, m),  
           Ia1_d2 = rep(0, m),  
           Ia2    = rep(0, m), 
           R      = (N_byage) * prop_immune,
           H      = rep(0, m),
           
           count_Ia1    = rep(0, m), 
           count_Ia2    = rep(0, m), 
           count_Is_udx = rep(0, m), 
           count_A      = rep(0, m)
)


## End ---------
