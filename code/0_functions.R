## Plotting theme --------

get_theme <- function(txt = 12, ...) {
  theme_classic() + theme(axis.text = element_text(size = txt),
                          axis.title = element_text(size = txt),
                          legend.text = element_text(size = txt),
                          legend.title = element_text(size = txt),
                          strip.text = element_text(size = txt),
                          panel.grid.major = element_line(),
                          panel.border = element_rect(fill = NA), 
                          ...)
}


## Rename scenarios ---------

rename_scenarios <- function(scenario, scenario_short, scenario_labels){
  scenario[scenario == scenario_short[1]] <- scenario_labels[1]
  scenario[scenario == scenario_short[2]]  <- scenario_labels[2]
  scenario[scenario == scenario_short[3]] <- scenario_labels[3]
  scenario[scenario == scenario_short[4]] <- scenario_labels[4]
  return(factor(scenario, levels = scenario_labels))
}


## Summarize output --------

get_summary <- function(dat, cols) {
  dat %>%
    summarise(across(.cols = {{ cols }},
                     .fns  = list(mean = ~mean(., na.rm = TRUE),
                                  lcl  = ~quantile(., 0.025, na.rm = TRUE),
                                  ucl  = ~quantile(., 0.975, na.rm = TRUE)),
                     .names = "{.col}_{.fn}"))
}


## Calculate beta -----------

# Only used to anchor beta in baseline scenario (without antivirals)

calculate_beta <- function(Reff, contacts, m = 5,
                           p_immune = 0,
                           asym_adj = 0.6,
                           p_a1, p_a2, av_trans_adj = 1
) {
  if(!(exists('contacts'))) stop("You first need to create the contact matrix by running 1_contact_matrix.R")
  
  gamma_a1 <- 1 / 1   # rate of leaving Ia1_d1 and Ia1_d2
  gamma_a2 <- 1 / 1   # rate of leaving Ia2
  
  K <- matrix(0, nrow = m, ncol = m)
  
  for (i in 1:m) {
    for (j in 1:m) {
      # ~ prob in compartment * duration of compartment
      A_contr <-  (1 - p_sym_corr[i]) * dur_asym
      Ip_contr <- p_sym_corr[i]  * dur_presym_inf
      
      # Is_d1 + Is_d2 + Is_d3
      Is_contr <- p_sym_corr[i]  * (1 * 1 + (1 - p_a1[i]) * 1 + (1 - p_a1[i] - p_a2[i]) * 1)   
      
      # Ia1_d1 + Ia1_d2 + Ia2
      Iav_contr <- p_sym_corr[i] * (p_a1[i] * (1/gamma_a1) + p_a1[i] * (1/gamma_a1) + p_a2[i] * (1/gamma_a2) ) 
      
      # ~ transmission weighting
      K[i, j] <- as.numeric(contacts[i, j] * N_byage[i]) * (Ip_contr + Is_contr + asym_adj * A_contr + av_trans_adj *  Iav_contr )
    }
  }
  
  eigenvalues <- eigen(K)$values
  Reff / (max(Re(eigenvalues)) * (1 - p_immune) )
}


## Combine output ----------------------

combine_out <- function(dat){
  tmp <- 
    dat %>% 
    mutate(# Compartments over time
      S_all = rowSums(across(matches(c("S[0-9]"), ignore.case = F))),
      E_all = rowSums(across(matches(c("E[0-9]"), ignore.case = F))),
      A_all = rowSums(across(matches(c("A[0-9]"), ignore.case = F))),
      Ip_all = rowSums(across(matches(c("Ip[0-9]"), ignore.case = F))),
      Is_all = rowSums(across(matches(c("Is_d[0-9]"), ignore.case = F))),
      Ia1_all = rowSums(across(matches(c("Ia1_d[0-9]"), ignore.case = F))),
      Ia2_all = rowSums(across(matches(c("Ia2[0-9]"), ignore.case = F))),
      R_all = rowSums(across(matches(c("R[0-9]"), ignore.case = F))),
      H_all = rowSums(across(matches(c("H[0-9]"), ignore.case = F))),
      
      # Counting cumulative infections as people transition to R
      count_A_all = rowSums(across(contains("count_A"))),
      count_Is_all = rowSums(across(contains("count_Is_udx"))),
      count_Ia1_all = rowSums(across(contains("count_Ia1"))),
      count_Ia2_all = rowSums(across(contains("count_Ia2"))),
      
      # Correct variables with too many columns summed
      A_all = A_all - count_A_all,
      Ia2_all = Ia2_all - count_Ia2_all,
      
      I_all = A_all + Ip_all + Is_all + Ia1_all + Ia2_all
    )
  tmp
}



# Model equations (one infectious sub-stage per day)  --------------------

flu.model.1stage <- function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    
    # Set-up states (m age groups within each)
    dS  <- dE <- dA <- dIp <- dR <- dH <- numeric(m)
    dIs_d1 <- dIs_d2 <- dIs_d3 <- numeric(m)
    dIa1_d1 <- dIa1_d2 <- dIa2 <- numeric(m)
    
    S     <- state[paste0("S", 1:m)]     # susceptible
    E     <- state[paste0("E", 1:m)]     # exposed
    A     <- state[paste0("A", 1:m)]     # infected, asymptomatic, undx
    Ip    <- state[paste0("Ip", 1:m)]    # infected, presymptomatic, undx
    Is_d1 <- state[paste0("Is_d1", 1:m)] # infected, symptomatic, undx, d1
    Is_d2 <- state[paste0("Is_d2", 1:m)] # infected, symptomatic, undx, d2
    Is_d3 <- state[paste0("Is_d3", 1:m)] # infected, symptomatic, undx, d3
    Ia1_d1 <- state[paste0("Ia1_d1", 1:m)]   # infected, on antivirals after 1d, d1
    Ia1_d2 <- state[paste0("Ia1_d2", 1:m)]   # infected, on antivirals, after 1d, d2
    Ia2   <- state[paste0("Ia2", 1:m)]   # infected, on antivirals, after 2d, d1
    R     <- state[paste0("R", 1:m)]     # recovered, not hospitalized
    H     <- state[paste0("H", 1:m)]     # hospitalized, hospitalized
    
    # Add seasonality 
    beta <- beta0 * (1 + .1*sin(2*pi*( (t - 45)/365))) # calibrated phase
    
    # Equations
    for (i in 1:m) {
      Is <- Is_d1 + Is_d2 + Is_d3 
      Ia <- Ia1_d1 + Ia1_d2 + Ia2
      
      dS[i]     <- -sum(S[i]* beta[i,]* (Is + Ip + adj_avi * Ia + adj_asym *A) )
      dE[i]     <-  sum(S[i]* beta[i,]* (Is + Ip + adj_avi * Ia + adj_asym *A) ) - sigma *E[i]
      dA[i]     <-  (1 - p_sym[i])* sigma *E[i] - gamma_a* A[i]                                  
      dIp[i]    <-  p_sym[i]* sigma *E[i] - gamma_p* Ip[i]                                  
      dIs_d1[i] <-  gamma_p* Ip[i] - gamma_s* Is_d1[i]  
      dIs_d2[i] <-  (1 - p_a1[i])* gamma_s* Is_d1[i] - gamma_s* Is_d2[i] 
      dIs_d3[i] <-  (1 - p_a2[i])* gamma_s* Is_d2[i] - gamma_s* Is_d3[i]  
      dIa1_d1[i] <-  p_a1[i]* gamma_s* Is_d1[i] - gamma_a1* Ia1_d1[i]
      dIa1_d2[i] <-  gamma_a1* Ia1_d1[i] - gamma_a1* Ia1_d2[i]
      dIa2[i]    <-  p_a2[i]* gamma_s* Is_d2[i] - gamma_a2* Ia2[i]
      dR[i]     <-  (1 - p_h[i]* r_h)* (gamma_a1* Ia1_d2[i] + gamma_a2* Ia2[i]) + (1 - p_h[i])* gamma_s* Is_d3[i] + gamma_a* A[i]
      dH[i]     <-  p_h[i]* (r_h* (gamma_a1* Ia1_d2[i] + gamma_a2* Ia2[i]) + gamma_s* Is_d3[i])
    }
    
    # Track transitions from infectious compartments
    count_Ia1 = gamma_a1*Ia1_d2
    count_Ia2 = gamma_a2*Ia2
    count_Is_udx = gamma_s* Is_d3
    count_A = gamma_a* A
    
    return(list(unlist(c(dS, dE, dA, dIp, dIs_d1, dIs_d2, dIs_d3, dIa1_d1, dIa1_d2, dIa2, dR, dH,
                         count_Ia1, count_Ia2, count_Is_udx, count_A
    ))))
  })
}


# End ----------

