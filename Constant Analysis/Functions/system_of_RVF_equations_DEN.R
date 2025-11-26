#--------------------------SYSTEM OF RVF EQUATIONS------------------------------
# not used within analysis - contains infection
system_of_RVF_equations_density_temp <- function(t, State, Pars) {
  
  with(as.list(c(State, Pars)), { # Variables from State & Pars brought into function
    
    # temperature_t <- temperature(t) # Converts temps to Kelvin for each day
    temperature_t <- constanttemp + 273.15 # Converts temps to Kelvin for each day when there is constant temperature
    
    # Mortality rates of species & life stage are dependent on temperature:
    
    mu_C_O <- mortality(temperature_t)[1] # Mortality rate for Culex eggs
    mu_C_L <- mortality(temperature_t)[2] # Mortality rate for Culex larvae
    mu_C_P <- mortality(temperature_t)[3] # Mortality rate for Culex pupae
    mu_C_E <- mortality(temperature_t)[4] # Mortality rate for additional Culex pupal (delta in paper)
    mu_C <- mortality(temperature_t)[5]   # Mortality rate for Culex adult
    
    mu_A_Oi <- mortality(temperature_t)[6] # Mortality rate for immature Aedes eggs
    mu_A_Om <- mortality(temperature_t)[7] # Mortality rate for mature Aedes eggs
    mu_A_L <- mortality(temperature_t)[8] # Mortality rate for Aedes larvae
    mu_A_P <- mortality(temperature_t)[9] # Mortality rate for Aedes pupae
    mu_A_E <- mortality(temperature_t)[10] # Mortality rate for additional Aedes pupal pupal (delta in paper)
    mu_A <- mortality(temperature_t)[11]  # Mortality rate for Aedes adult
    
    # Culex developmental rates (egg, larva & pupa) dependent on temperature:
    
    theta_C_O <- 2 * Schoolfield(temperature_t, Culex_Embryogenesis) #assumed to be half the time for the first instar
    theta_C_L <- Schoolfield(temperature_t, Culex_Larval_development)
    theta_C_P <- Schoolfield(temperature_t, Culex_Pupal_development)
    # based on the regression of Reisen, W. K., Milby, M. M., Presser, S. B., & Hardy, J. L. (1992). Ecology of mosquitoes and St. Louis encephalitis virus in the Los Angeles Basin of California, 1987-1990. Journal of medical entomology, 29(4), 582–98. table 3 figure 5 for Culex quinquefasciatus. Maybe you could use Culex tartalis or their mean (See their discussion about lab/field dicrepancy on page 588

    # Culex developmental rates (adults) with infinite blood meal availability:
    
    theta_C_A1_no_lim <- 0.0173*(temperature_t - 273.16 - 9.60) # eq. 12 GL2018 supplementary - from Fischer no matter Kelvin of celsius as we consider differences
    theta_C_A2_no_lim <- theta_C_A1_no_lim/reduction_C

    # Aedes dev rates (egg, larva & pupa) more complicated
    # Number of hatching eggs at time t-k will be null if:
             # k is less than the minimum desiccation period
             # Eggs were submerged before achieving min desiccation period
             
    Td <- 6 # Minimum desiccation period in days

    delta_sp <- 0.197  # proportion of spontaneous hatching without flooding
    
    # 1/rate = actual time for eggs to hatch, the bigger of:
    theta_A_O_inv <- pmax(Td, 1/Schoolfield(temperature_t, Aedes_Embryogenesis)) # eq. 15 GL2018 supplementary
    
    # Aedes developmental rates (egg, larva & pupa) dependent on temperature:
    
    theta_A_O <- 1/theta_A_O_inv
    theta_A_L <- Schoolfield(temperature_t, Aedes_Larval_development)
    theta_A_P <- Schoolfield(temperature_t, Aedes_Pupal_development)
    
    # Aedes developmental rates (adults) with infinite blood meal availability:
    
    theta_A_A1_no_lim <- 0.0173*(temperature_t - 273.16 - 9.60) # eq. 20 A1 GL2018 supplementary - from Fischer no matter Kelvin of celsius as we consider differences
    theta_A_A2_no_lim <- theta_A_A1_no_lim/reduction_A

    a_Culex <- biting_rate(temperature_t)[[1]]
    a_Aedes <- biting_rate(temperature_t)[[2]]
    
    epsilon_C <- EIP(temperature_t)[[1]] # Culex extrinsic incubation period
    epsilon_A <- EIP(temperature_t)[[2]] # Aedes extrinsic incubation period

    ov_C <- oviposition_rate_Culex(t, Pars)[[1]] # Oviposition rate Culex
    b_C <- oviposition_rate_Culex(t, Pars)[[2]] # Eggs per batch Culex
    K_C <- oviposition_rate_Culex(t, Pars)[[3]] # Carrying capacity eggs in pond
    
    ov_A <- oviposition_rate_Aedes(t, Pars)[[1]] # Oviposition rate Aedes
    b_A <- oviposition_rate_Aedes(t, Pars)[[2]] # Eggs per batch Aedes
    K_A <- oviposition_rate_Aedes(t, Pars)[[3]] # Carrying capacity eggs in pond
    tau_A_O2 <- oviposition_rate_Aedes(t, Pars)[[4]] # Hatching rate Aedes
    
    # Total livestock = sum of susceptible, exposed, infected and recovered:    
    N_L1 <- livestocktotal
    
    # Total Culex = only adults female C1 and C2 are biting
    N_Culex <- (A1_c + A2_c + I_A2_c + E2_c)
    
    # Infected Culex = proportion of females to find host species x infected
    Inf_Culex <- prop_find*(I_A2_c)
	  Inf_Culex2 <- prop_find*(I_A2_c)
		
	  # Vector to host ratio for Culex:
    m_C_L1 <- prop_find*(N_Culex/N_L1)

    # Total Aedes = only adults female A1 and A2 are biting
    N_Aedes <- (A1_a + A2_a + I_A1_a + I_A2_a + E2_a)
    
    Inf_Aedes <- prop_find*(I_A1_a + I_A2_a)
		Inf_Aedes1 <- prop_find*(I_A1_a)
		Inf_Aedes2 <- prop_find*(I_A2_a)
		
		# If sum of eggs is less than zero, automatically set it to zero
		# Then, if sum of eggs exceeds carrying capacity, set to carrying capacity:
    Egg_Aedes <- min(max(O_a1 + O_a2 + I_O_a1 + I_O_a2, 0), K_A)
    
    # Vector to host ratio for Aedes:    
    m_A_L1 <- prop_find*(N_Aedes/N_L1)
    
    # These are the absolute developmental rates for adults: Take into account the vector to host ratio
    theta_C_A1 <- theta_C_A1_no_lim/(1 + (m_A_L1 + m_C_L1)/q_divided)
		theta_C_A2 <- theta_C_A2_no_lim/(1 + (m_A_L1 + m_C_L1)/q_divided)
		theta_A_A1 <- theta_A_A1_no_lim/(1 + (m_A_L1 + m_C_L1)/q_divided)
		theta_A_A2 <- theta_A_A2_no_lim/(1 + (m_A_L1 + m_C_L1)/q_divided)
    
		
    phi_L1_C <- 1 # proportion of meal, if multiple host need to use Sota model
    a_L1_C <- a_Culex*phi_L1_C # biting rate since phi = 1
    a_L1_C_dens <- a_L1_C/(1 + (m_A_L1 + m_C_L1)/q_divided)
    
	  # This is force of infection from livestock to Culex:    
	  lambda_L1_C_dens <- beta_L1_C*a_L1_C_dens*I_L1/N_L1
    lambda_L1_C_dens2 <- lambda_L1_C_dens/reduction_C

    # Eggs laid per batch for Culex
    b_C_dens <- b_C/(1 + (m_A_L1 + m_C_L1)/q_divided)

    phi_L1_A <- 1 # proportion of meal, if multiple host need to use Sota model                             
    a_L1_A <- a_Aedes*phi_L1_A # biting rate since phi = 1                   
    a_L1_A_dens <- a_L1_A/(1 + (m_A_L1 + m_C_L1)/q_divided)                
    
    # This is force of infection from livestock to Aedes:    
    lambda_L1_A_dens <- beta_L1_A*a_L1_A_dens*I_L1/N_L1           
    lambda_L1_A_dens2 <- lambda_L1_A_dens/reduction_A
    
    # Eggs laid per batch for Aedes   
    b_A_dens <- b_A/(1 + (m_A_L1 + m_C_L1)/q_divided)
    
    lambda_A_L1 <- beta_A_L1*a_L1_A_dens*(Inf_Aedes1 + 
                                            Inf_Aedes2/reduction_A)/N_L1
    lambda_C_L1 <- beta_C_L1*a_L1_C_dens*(Inf_Culex2/reduction_C)/N_L1    

    # DIFFERENTIAL EQUATIONS ---------------------------------------------------
    
    # Number entering compartment per day = positive
    # Number exiting compartment per day = negative
    
    # AEDES SUSCEPTIBLE
    
    # Immature eggs:
    dO_a1 <- b_A_dens*(1 - Egg_Aedes/K_A)*ov_A*(F_a + (1 - q_A)*(I_F_a+E_F_a)) - mu_A_Oi*O_a1 - theta_A_O*O_a1
    # Mature eggs:
    dO_a2 <- (1 - delta_sp)*theta_A_O*O_a1 - mu_A_Om*O_a2 - tau_A_O2*O_a2
    # Larvae:
    dL_a <- tau_A_O2*O_a2 + theta_A_O*O_a1*delta_sp - mu_A_L*L_a - theta_A_L*L_a
    # Pupae:
    dP_a <- theta_A_L*L_a - mu_A_P*P_a - theta_A_P*P_a
    # Nulliparous adults:
    dA1_a <- theta_A_P*(mu_A_E/2)*P_a - mu_A*A1_a - theta_A_A1*A1_a
    # Flyers:
    dF_a <- theta_A_A1*A1_a + theta_A_A2*A2_a - ov_A*F_a - mu_A*F_a - lambda_L1_A_dens*A1_a - lambda_L1_A_dens2*A2_a
    # Parous adults:
    dA2_a <- ov_A*F_a - mu_A*A2_a - theta_A_A2*A2_a  
    
    #--------------------------------------------------------------------------#
    
    # AEDES EXPOSED
    
    # Flyers:
    dE_F_a <- theta_A_A2*E2_a + lambda_L1_A_dens*A1_a + lambda_L1_A_dens2*A2_a - ov_A*E_F_a - mu_A*E_F_a - epsilon_A*E_F_a 
    # Parous adults:
	  dE2_a <- ov_A*E_F_a - epsilon_A*E2_a - mu_A*E2_a - mu_RVF_A*E2_a - theta_A_A2*E2_a
    
	  #--------------------------------------------------------------------------#
	  
	  # AEDES INFECTIOUS
	  
	  # Immature eggs:    
		dI_O_a1 <- b_A_dens*(1 - Egg_Aedes/K_A)*ov_A*q_A*(I_F_a + E_F_a) - mu_A_Oi*I_O_a1 - theta_A_O*I_O_a1
		# Mature eggs:
    dI_O_a2 <- (1 - delta_sp)*theta_A_O*I_O_a1 - mu_A_Om*I_O_a2 - tau_A_O2*I_O_a2
    # Larvae:
    dI_L_a <- tau_A_O2*I_O_a2 + theta_A_O*I_O_a1*delta_sp - mu_A_L*I_L_a - theta_A_L*I_L_a
    # Pupae:
    dI_P_a <- theta_A_L*I_L_a - mu_A_P*I_P_a - theta_A_P*I_P_a
    # Nulliparous adults:
    dI_A1_a <- theta_A_P*(mu_A_E/2)*I_P_a - mu_A*I_A1_a - theta_A_A1*I_A1_a
    # Flyers:
    dI_F_a <- theta_A_A1*I_A1_a + theta_A_A2*I_A2_a + epsilon_A*E_F_a - ov_A*I_F_a - mu_A*I_F_a 
    # Parous adults:
    dI_A2_a <- ov_A*I_F_a + epsilon_A*E2_a - mu_A*I_A2_a - theta_A_A2*I_A2_a
        
    #--------------------------------------------------------------------------#
    
    # CULEX SUSCEPTIBLE:
    
    # Eggs:
    dO_c <- b_C_dens*(1 - min(O_c/K_C, 1))*ov_C*(F_c + (1 - 0.21)*I_F_c) - mu_C_O*O_c - theta_C_O*O_c #think missing F_C_exposed in first term
    # Larvae:
    dL_c <- theta_C_O*O_c - mu_C_L*L_c - theta_C_L*L_c
    # Pupae:
    dP_c <- theta_C_L*L_c - mu_C_P*P_c - theta_C_P*P_c
    # Nulliparous adults:
    dA1_c <- theta_C_P*(mu_C_E/2)*P_c - mu_C*A1_c - theta_C_A1*A1_c
    # Flyers:
    dF_c <- theta_C_A1*A1_c + theta_C_A2*A2_c - lambda_L1_C_dens*A1_c - lambda_L1_C_dens2*A2_c - ov_C*F_c - mu_C*F_c
    # Parous adults:
    dA2_c <- ov_C*F_c - mu_C*A2_c - theta_C_A2*A2_c
    
    #--------------------------------------------------------------------------#
    
    # CULEX EXPOSED:
    
    # Flyers:    
    dE_F_c <- lambda_L1_C_dens*A1_c + lambda_L1_C_dens2*A2_c - ov_C*E_F_c - mu_C*E_F_c - mu_RVF_C* - epsilon_C*E_F_c
    # Parous adults:
    dE2_c <- ov_C*E_F_c - epsilon_C*E2_c - mu_C*E2_c - mu_RVF_C*E2_c - theta_C_A2*E2_c
		
    #--------------------------------------------------------------------------#
    
    # CULEX INFECTIOUS:
    
    # Flyers:
		dI_F_c <- theta_C_A2*I_A2_c + epsilon_C*E_F_c - ov_C*I_F_c - mu_C*I_F_c - mu_RVF_C*I_F_c
		# Parous adults:
    dI_A2_c <- ov_C*I_F_c + epsilon_C*E2_c - mu_C*I_A2_c - theta_C_A2*I_A2_c - mu_RVF_C*I_A2_c

    #--------------------------------------------------------------------------#
    
    # LIVESTOCK SUSCEPTIBLE:
    
    dS_L1 <- b_L1*N_L1 - mu_L1*S_L1 - lambda_A_L1*S_L1 - lambda_C_L1*S_L1
    
    # LIVESTOCK EXPOSED:
    
    dE_L1 <- lambda_A_L1*S_L1 + lambda_C_L1*S_L1 - epsilon_L1*E_L1 - mu_L1*E_L1

    # LIVESTOCK INFECTIOUS:
    
    dI_L1 <- epsilon_L1*E_L1 - gamma_L1*I_L1 - mu_L1*I_L1 - mu_RVF_L1*I_L1
    
    # LIVESTOCK RECOVERED:
    
    dR_L1 <- gamma_L1*I_L1 - mu_L1*R_L1
    
    #--------------------------------------------------------------------------#
    
    return(list(c(dO_c, dL_c,
                  dP_c, dA1_c,
                  dF_c, dA2_c,
                  dE2_c, dI_F_c,
                  dI_A2_c, dO_a1,
                  dO_a2, dL_a,
                  dP_a, dA1_a,
                  dF_a, dA2_a,
                  dE2_a, dI_O_a1,
                  dI_O_a2, dI_L_a,
                  dI_P_a, dI_A1_a,
                  dI_F_a, dI_A2_a,
                  dS_L1, dE_L1,
                  dI_L1, dR_L1,
                  dE_F_c, dE_F_a)))
    
        }
		
	  )
    
}

#--------------------------SYSTEM OF RVF EQUATIONS - Culex Only------------------------------

system_of_RVF_equations_density_temp_Culex <- function(t, State, Pars) {
  
  with(as.list(c(State, Pars)), { # Variables from State & Pars brought into function
    
    # temperature_t <- temperature(t) # Converts temps to Kelvin for each day
    temperature_t <- constanttemp + 273.15 # Converts temps to Kelvin for each day when there is constant temperature
    
    # Mortality rates of species & life stage are dependent on temperature:
    
    mu_C_O <- mortality(temperature_t)[1] # Mortality rate for Culex eggs
    mu_C_L <- mortality(temperature_t)[2] # Mortality rate for Culex larvae
    mu_C_P <- mortality(temperature_t)[3] # Mortality rate for Culex pupae
    mu_C_E <- mortality(temperature_t)[4] # Mortality rate for additional Culex pupal (delta in paper)
    mu_C <- mortality(temperature_t)[5]   # Mortality rate for Culex adult
    
    # Culex developmental rates (egg, larva & pupa) dependent on temperature:
    
    theta_C_O <- 2 * Schoolfield(temperature_t, Culex_Embryogenesis) #assumed to be half the time for the first instar
    theta_C_L <- Schoolfield(temperature_t, Culex_Larval_development)
    theta_C_P <- Schoolfield(temperature_t, Culex_Pupal_development)
    # based on the regression of Reisen, W. K., Milby, M. M., Presser, S. B., & Hardy, J. L. (1992). Ecology of mosquitoes and St. Louis encephalitis virus in the Los Angeles Basin of California, 1987-1990. Journal of medical entomology, 29(4), 582–98. table 3 figure 5 for Culex quinquefasciatus. Maybe you could use Culex tartalis or their mean (See their discussion about lab/field dicrepancy on page 588
    
    # Culex developmental rates (adults) with infinite blood meal availability:
    
    theta_C_A1_no_lim <- 0.0173*(temperature_t - 273.16 - 9.60) # eq. 12 GL2018 supplementary - from Fischer no matter Kelvin of celsius as we consider differences
    theta_C_A2_no_lim <- theta_C_A1_no_lim/reduction_C
    
    a_Culex <- biting_rate(temperature_t)[[1]]

    epsilon_C <- EIP(temperature_t)[[1]] # Culex extrinsic incubation period

    ov_C <- oviposition_rate_Culex(t, Pars)[[1]] # Oviposition rate Culex
    b_C <- oviposition_rate_Culex(t, Pars)[[2]] # Eggs per batch Culex
    K_C <- oviposition_rate_Culex(t, Pars)[[3]] # Carrying capacity eggs in pond
    
    # Total livestock = sum of susceptible, exposed, infected and recovered:    
    N_L1 <- livestocktotal
    
    # Total Culex = only adults female C1 and C2 are biting
    N_Culex <- (A1_c + A2_c)
    
    # Vector to host ratio for Culex:
    m_C_L1 <- prop_find*(N_Culex/N_L1)
    
    # Total Aedes = only adults female A1 and A2 are biting
    N_Aedes <- 0
    
    # These are the absolute developmental rates for adults: Take into account the vector to host ratio
    theta_C_A1 <- theta_C_A1_no_lim/(1 + ( m_C_L1)/q_divided)
    theta_C_A2 <- theta_C_A2_no_lim/(1 + ( m_C_L1)/q_divided)

    phi_L1_C <- 1 # proportion of meal, if multiple host need to use Sota model
    a_L1_C <- a_Culex*phi_L1_C # biting rate since phi = 1
    a_L1_C_dens <- a_L1_C/(1 + ( m_C_L1)/q_divided)
    
    # Eggs laid per batch for Culex
    b_C_dens <- b_C/(1 + ( m_C_L1)/q_divided)
    
    # DIFFERENTIAL EQUATIONS ---------------------------------------------------
    
    # Number entering compartment per day = positive
    # Number exiting compartment per day = negative
    
    # CULEX SUSCEPTIBLE:
    
    # Eggs:
    dO_c <- b_C_dens*(1 - min(O_c/K_C, 1))*ov_C*(F_c) - mu_C_O*O_c - theta_C_O*O_c
    # Larvae:
    dL_c <- theta_C_O*O_c - mu_C_L*L_c - theta_C_L*L_c
    # Pupae:
    dP_c <- theta_C_L*L_c - mu_C_P*P_c - theta_C_P*P_c
    # Nulliparous adults:
    dA1_c <- theta_C_P*(mu_C_E/2)*P_c - mu_C*A1_c - theta_C_A1*A1_c
    # Flyers:
    dF_c <- theta_C_A1*A1_c + theta_C_A2*A2_c - ov_C*F_c - mu_C*F_c
    # Parous adults:
    dA2_c <- ov_C*F_c - mu_C*A2_c - theta_C_A2*A2_c
    

    
    #--------------------------------------------------------------------------#
    
    return(list(c(dO_c, dL_c,
                  dP_c, dA1_c,
                  dF_c, dA2_c )))
    
  }
  
  )
  
}

#--------------------------SYSTEM OF RVF EQUATIONS - Aedes Only------------------------------

system_of_RVF_equations_density_temp_Aedes <- function(t, State, Pars) {
  
  with(as.list(c(State, Pars)), { # Variables from State & Pars brought into function
    
    # temperature_t <- temperature(t) # Converts temps to Kelvin for each day
    temperature_t <- constanttemp + 273.15 # Converts temps to Kelvin for each day when there is constant temperature
    
    # Mortality rates of species & life stage are dependent on temperature:
    
    mu_A_Oi <- mortality(temperature_t)[6] # Mortality rate for immature Aedes eggs
    mu_A_Om <- mortality(temperature_t)[7] # Mortality rate for mature Aedes eggs
    mu_A_L <- mortality(temperature_t)[8] # Mortality rate for Aedes larvae
    mu_A_P <- mortality(temperature_t)[9] # Mortality rate for Aedes pupae
    mu_A_E <- mortality(temperature_t)[10] # Mortality rate for additional Aedes pupal pupal (delta in paper)
    mu_A <- mortality(temperature_t)[11]  # Mortality rate for Aedes adult
    
    # Aedes dev rates (egg, larva & pupa) more complicated
    # Number of hatching eggs at time t-k will be null if:
    # k is less than the minimum desiccation period
    # Eggs were submerged before achieving min desiccation period
    
    Td <- 6 # 6 Minimum desiccation period in days

    delta_sp <- 0.197 #0.197 # proportion of spontaneous hatching without flooding

    # 1/rate = actual time for eggs to hatch, the bigger of:
    theta_A_O_inv <- pmax(Td, 1/Schoolfield(temperature_t, Aedes_Embryogenesis)) # eq. 15 GL2018 supplementary
    
    # Aedes developmental rates (egg, larva & pupa) dependent on temperature:
    
    theta_A_O <- 1/theta_A_O_inv
    theta_A_L <- Schoolfield(temperature_t, Aedes_Larval_development)
    theta_A_P <- Schoolfield(temperature_t, Aedes_Pupal_development)
    
    # Aedes developmental rates (adults) with infinite blood meal availability:
    
    theta_A_A1_no_lim <- 0.0173*(temperature_t - 273.16 - 9.60) # eq. 20 A1 GL2018 supplementary - from Fischer no matter Kelvin of celsius as we consider differences
    theta_A_A2_no_lim <- theta_A_A1_no_lim/reduction_A
    
    a_Culex <- biting_rate(temperature_t)[[1]]
    a_Aedes <- biting_rate(temperature_t)[[2]]
    
    ov_A <- oviposition_rate_Aedes(t, Pars)[[1]] # Oviposition rate Aedes
    b_A <- oviposition_rate_Aedes(t, Pars)[[2]] # Eggs per batch Aedes
    K_A <- oviposition_rate_Aedes(t, Pars)[[3]] # Carrying capacity eggs in pond
    tau_A_O2 <- oviposition_rate_Aedes(t, Pars)[[4]] # Hatching rate Aedes
    
    N_L1 <- livestocktotal
    
    # Total Aedes = only adults female A1 and A2 are biting
    N_Aedes <- (A1_a + A2_a)
    
    # If sum of eggs is less than zero, automatically set it to zero
    # Then, if sum of eggs exceeds carrying capacity, set to carrying capacity:
    Egg_Aedes <- min(max(O_a1 + O_a2, 0), K_A)
    
    # Vector to host ratio for Aedes:    
    m_A_L1 <- prop_find*(N_Aedes/N_L1) #eq 11
    
    # These are the absolute developmental rates for adults: Take into account the vector to host ratio
    theta_A_A1 <- theta_A_A1_no_lim/(1 + (m_A_L1)/q_divided) #eq 13
    theta_A_A2 <- theta_A_A2_no_lim/(1 + (m_A_L1)/q_divided) #eq 13
    
    # Eggs laid per batch for Aedes   
    b_A_dens <- b_A/(1 + (m_A_L1)/q_divided) #eq 10
    
    # DIFFERENTIAL EQUATIONS ---------------------------------------------------
    
    # Number entering compartment per day = positive
    # Number exiting compartment per day = negative
    
    # Immature eggs:
    dO_a1 <- b_A_dens*(1 - Egg_Aedes/K_A)*ov_A*(F_a) - mu_A_Oi*O_a1 - theta_A_O*O_a1
    # Mature eggs:
    dO_a2 <- (1 - delta_sp)*theta_A_O*O_a1 - mu_A_Om*O_a2 - tau_A_O2*O_a2
    # Larvae:
    dL_a <- tau_A_O2*O_a2 - mu_A_L*L_a - theta_A_L*L_a + delta_sp*theta_A_O*O_a1
    # Pupae:
    dP_a <- theta_A_L*L_a - mu_A_P*P_a - theta_A_P*P_a
    # Nulliparous adults:
    dA1_a <- theta_A_P*(mu_A_E/2)*P_a - mu_A*A1_a - theta_A_A1*A1_a
    # Flyers:
    dF_a <- theta_A_A1*A1_a + theta_A_A2*A2_a - ov_A*F_a - mu_A*F_a
    # Parous adults:
    dA2_a <- ov_A*F_a - mu_A*A2_a - theta_A_A2*A2_a
    
    return(list(c(dO_a1,
                  dO_a2, dL_a,
                  dP_a, dA1_a,
                  dF_a, dA2_a)))
    
  }
  
  )
  
}


#--------------------------SYSTEM OF RVF EQUATIONS - Culex Only Periodic------------------------------

system_of_RVF_equations_density_periodic_Culex <- function(t, State, Pars) {
  
  with(as.list(c(State, Pars)), { # Variables from State & Pars brought into function
    
    temperature_t <- temperature_periodic(t) # Converts temps to Kelvin for each day
    # temperature_t <- constanttemp + 273.15 # Converts temps to Kelvin for each day when there is constant temperature
    
    # Mortality rates of species & life stage are dependent on temperature:
    
    mu_C_O <- mortality(temperature_t)[1] # Mortality rate for Culex eggs
    mu_C_L <- mortality(temperature_t)[2] # Mortality rate for Culex larvae
    mu_C_P <- mortality(temperature_t)[3] # Mortality rate for Culex pupae
    mu_C_E <- mortality(temperature_t)[4] # Mortality rate for additional Culex pupal (delta in paper)
    mu_C <- mortality(temperature_t)[5]   # Mortality rate for Culex adult
    
    # Culex developmental rates (egg, larva & pupa) dependent on temperature:
    
    theta_C_O <- 2 * Schoolfield(temperature_t, Culex_Embryogenesis) #assumed to be half the time for the first instar
    theta_C_L <- Schoolfield(temperature_t, Culex_Larval_development)
    theta_C_P <- Schoolfield(temperature_t, Culex_Pupal_development)
    # based on the regression of Reisen, W. K., Milby, M. M., Presser, S. B., & Hardy, J. L. (1992). Ecology of mosquitoes and St. Louis encephalitis virus in the Los Angeles Basin of California, 1987-1990. Journal of medical entomology, 29(4), 582–98. table 3 figure 5 for Culex quinquefasciatus. Maybe you could use Culex tartalis or their mean (See their discussion about lab/field dicrepancy on page 588
    
    # Culex developmental rates (adults) with infinite blood meal availability:
    
    theta_C_A1_no_lim <- 0.0173*(temperature_t - 273.16 - 9.60) # eq. 12 GL2018 supplementary - from Fischer no matter Kelvin of celsius as we consider differences
    theta_C_A2_no_lim <- theta_C_A1_no_lim/reduction_C
    
    a_Culex <- biting_rate(temperature_t)[[1]]
    
    epsilon_C <- EIP(temperature_t)[[1]] # Culex extrinsic incubation period
    
    ov_C <- oviposition_rate_periodic_Culex(t, Pars)[[1]] # Oviposition rate Culex
    b_C <- oviposition_rate_periodic_Culex(t, Pars)[[2]] # Eggs per batch Culex
    K_C <- oviposition_rate_periodic_Culex(t, Pars)[[3]] # Carrying capacity eggs in pond
    
    # Total livestock = sum of susceptible, exposed, infected and recovered:    
    N_L1 <- livestocktotal
    
    # Total Culex = only adults female C1 and C2 are biting
    N_Culex <- (A1_c + A2_c)
    
    # Vector to host ratio for Culex:
    m_C_L1 <- prop_find*(N_Culex/N_L1)
    
    # Total Aedes = only adults female A1 and A2 are biting
    N_Aedes <- 0
    
    # These are the absolute developmental rates for adults: Take into account the vector to host ratio
    theta_C_A1 <- theta_C_A1_no_lim/(1 + ( m_C_L1)/q_divided)
    theta_C_A2 <- theta_C_A2_no_lim/(1 + ( m_C_L1)/q_divided)
    
    phi_L1_C <- 1 # proportion of meal, if multiple host need to use Sota model
    a_L1_C <- a_Culex*phi_L1_C # biting rate since phi = 1
    a_L1_C_dens <- a_L1_C/(1 + ( m_C_L1)/q_divided)
    
    # Eggs laid per batch for Culex
    b_C_dens <- b_C/(1 + ( m_C_L1)/q_divided)
    
    # DIFFERENTIAL EQUATIONS ---------------------------------------------------
    
    # Number entering compartment per day = positive
    # Number exiting compartment per day = negative
    
    # CULEX SUSCEPTIBLE:
    
    # Eggs:
    dO_c <- b_C_dens*(1 - min(O_c/K_C, 1))*ov_C*(F_c) - mu_C_O*O_c - theta_C_O*O_c
    # Larvae:
    dL_c <- theta_C_O*O_c - mu_C_L*L_c - theta_C_L*L_c
    # Pupae:
    dP_c <- theta_C_L*L_c - mu_C_P*P_c - theta_C_P*P_c
    # Nulliparous adults:
    dA1_c <- theta_C_P*(mu_C_E/2)*P_c - mu_C*A1_c - theta_C_A1*A1_c
    # Flyers:
    dF_c <- theta_C_A1*A1_c + theta_C_A2*A2_c - ov_C*F_c - mu_C*F_c
    # Parous adults:
    dA2_c <- ov_C*F_c - mu_C*A2_c - theta_C_A2*A2_c
    
    
    
    #--------------------------------------------------------------------------#
    
    return(list(c(dO_c, dL_c,
                  dP_c, dA1_c,
                  dF_c, dA2_c )))
    
  }
  
  )
  
}

#--------------------------SYSTEM OF RVF EQUATIONS - Aedes Only Periodic------------------------------

system_of_RVF_equations_density_periodic_Aedes <- function(t, State, Pars) {
  
  with(as.list(c(State, Pars)), { # Variables from State & Pars brought into function
    
    temperature_t <- temperature_periodic(t) # Converts temps to Kelvin for each day
    # temperature_t <- constanttemp + 273.15 # Converts temps to Kelvin for each day when there is constant temperature
    
    # Mortality rates of species & life stage are dependent on temperature:
    
    mu_A_Oi <- mortality(temperature_t)[6] # Mortality rate for immature Aedes eggs
    mu_A_Om <- mortality(temperature_t)[7] # Mortality rate for mature Aedes eggs
    mu_A_L <- mortality(temperature_t)[8] # Mortality rate for Aedes larvae
    mu_A_P <- mortality(temperature_t)[9] # Mortality rate for Aedes pupae
    mu_A_E <- mortality(temperature_t)[10] # Mortality rate for additional Aedes pupal pupal (delta in paper)
    mu_A <- mortality(temperature_t)[11]  # Mortality rate for Aedes adult
    
    
    # Aedes dev rates (egg, larva & pupa) more complicated
    # Number of hatching eggs at time t-k will be null if:
    # k is less than the minimum desiccation period
    # Eggs were submerged before achieving min desiccation period
    
    Td <- 6 # 6 Minimum desiccation period in days

    delta_sp <- 0.197 #0.197 # proportion of spontaneous hatching without flooding
    
    
    # 1/rate = actual time for eggs to hatch, the bigger of:
    theta_A_O_inv <- pmax(Td, 1/Schoolfield(temperature_t, Aedes_Embryogenesis)) # eq. 15 GL2018 supplementary
    
    # Aedes developmental rates (egg, larva & pupa) dependent on temperature:
    
    theta_A_O <- 1/theta_A_O_inv
    theta_A_L <- Schoolfield(temperature_t, Aedes_Larval_development)
    theta_A_P <- Schoolfield(temperature_t, Aedes_Pupal_development)
    
    # Aedes developmental rates (adults) with infinite blood meal availability:
    
    theta_A_A1_no_lim <- 0.0173*(temperature_t - 273.16 - 9.60) # eq. 20 A1 GL2018 supplementary - from Fischer no matter Kelvin of celsius as we consider differences
    theta_A_A2_no_lim <- theta_A_A1_no_lim/reduction_A
    
    a_Culex <- biting_rate(temperature_t)[[1]]
    a_Aedes <- biting_rate(temperature_t)[[2]]
    
    epsilon_C <- EIP(temperature_t)[[1]] # Culex extrinsic incubation period
    epsilon_A <- EIP(temperature_t)[[2]] # Aedes extrinsic incubation period
    
    ov_A <- oviposition_rate__periodic_Aedes(t, Pars)[[1]] # Oviposition rate Aedes
    b_A <- oviposition_rate__periodic_Aedes(t, Pars)[[2]] # Eggs per batch Aedes
    K_A <- oviposition_rate__periodic_Aedes(t, Pars)[[3]] # Carrying capacity eggs in pond
    tau_A_O2 <- oviposition_rate__periodic_Aedes(t, Pars)[[4]] # Hatching rate Aedes
    
    # Total livestock   
    # N_L1 <- livestocktotal
    
    # Total Aedes = only adults female A1 and A2 are biting
    N_Aedes <- (A1_a + A2_a)
    
    # If sum of eggs is less than zero, automatically set it to zero
    # Then, if sum of eggs exceeds carrying capacity, set to carrying capacity:
    Egg_Aedes <- min(max(O_a1 + O_a2, 0), K_A)
    
    # Vector to host ratio for Aedes:    
    m_A_L1 <- prop_find*(N_Aedes/N_L1)
    
    # These are the absolute developmental rates for adults: Take into account the vector to host ratio
    theta_A_A1 <- theta_A_A1_no_lim/(1 + (m_A_L1)/q_divided)
    theta_A_A2 <- theta_A_A2_no_lim/(1 + (m_A_L1)/q_divided)
    
    phi_L1_C <- 1 # proportion of meal, if multiple host need to use Sota model
    a_L1_C <- a_Culex*phi_L1_C # biting rate since phi = 1
    a_L1_C_dens <- a_L1_C/(1 + (m_A_L1)/q_divided)
    
    phi_L1_A <- 1 # proportion of meal, if multiple host need to use Sota model                             
    a_L1_A <- a_Aedes*phi_L1_A # biting rate since phi = 1                   
    a_L1_A_dens <- a_L1_A/(1 + (m_A_L1)/q_divided)                
    
    # Eggs laid per batch for Aedes   
    b_A_dens <- b_A/(1 + (m_A_L1)/q_divided)
    
    # DIFFERENTIAL EQUATIONS ---------------------------------------------------
    
    # Number entering compartment per day = positive
    # Number exiting compartment per day = negative
    
    # Immature eggs:
    dO_a1 <- b_A_dens*(1 - Egg_Aedes/K_A)*ov_A*(F_a) - mu_A_Oi*O_a1 - theta_A_O*O_a1
    # Mature eggs:
    dO_a2 <- (1 - delta_sp)*theta_A_O*O_a1 - mu_A_Om*O_a2 - tau_A_O2*O_a2
    # Larvae:
    dL_a <- tau_A_O2*O_a2 + theta_A_O*O_a1*delta_sp - mu_A_L*L_a - theta_A_L*L_a
    # Pupae:
    dP_a <- theta_A_L*L_a - mu_A_P*P_a - theta_A_P*P_a
    # Nulliparous adults:
    dA1_a <- theta_A_P*(mu_A_E/2)*P_a - mu_A*A1_a - theta_A_A1*A1_a
    # Flyers:
    dF_a <- theta_A_A1*A1_a + theta_A_A2*A2_a - ov_A*F_a - mu_A*F_a
    # Parous adults:
    dA2_a <- ov_A*F_a - mu_A*A2_a - theta_A_A2*A2_a  
    
    return(list(c(dO_a1,
                  dO_a2, dL_a,
                  dP_a, dA1_a,
                  dF_a, dA2_a)))
  }
  
  )
  
}
