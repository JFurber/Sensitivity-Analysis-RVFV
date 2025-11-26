#--------------------------------SCHOOLFIELD------------------------------------

# Developmental rates of mosquitoes from one stage to the next
# Note - larva moult four times before pupation - called instar stages
# Reference (13) in GL2018 supplementary 
# Also see eq (7) in Schoolfield 1981
# Stage parameters in Sharpe_DeMichele_parameters.R, also seen in Table S6 in GL2018 supplementary

Schoolfield <- function(temperature, stage_parameters) {
  
  # Maturation process is controlled by an enzyme
  # Enzyme is active in given temperature range (temperature K = Kelvin)
  # Deactivated only at high temperatures
  # Mean development rate formula involves:
  
    # Developmental rate assuming no enzyme inactivation at temperature 298K:
    RH025 <- stage_parameters[1]
    
    # Changes in thermodynamics enthalpies characteristic of the organism:
    HA <- stage_parameters[2] 
    HH <- stage_parameters[4]
    
    # Temperature when half of the enzyme is deactivated from high temperatures
    TH <- stage_parameters[3]
 
    # Universal Gas Constant:
    R_gas <- 1.987

    # Unsure about this following point:

  # See Equation (14) in Supplementary PNAS materials GL2018 for mean development rate formula
  # This is just breaking that into digestible chunks:
  
  term1 <- RH025*(temperature/298)
    
  arg_term2 <- (HA/R_gas)*(1/298 - 1/temperature)
  term2 <- exp(arg_term2)

  arg_term3 <- (HH/R_gas)*(1/TH - 1/temperature)
  term3 <- exp(arg_term3)

  # Mean development rate at temperature = T:
  rate <- term1*term2/(1 + term3)
  
  return(rate)

}

d_rate_d_temp <- function(temperature, stage_parameters) {
  
  RH025 <- stage_parameters[1]
  HA <- stage_parameters[2]
  TH <- stage_parameters[3]
  HH <- stage_parameters[4]
    
  R_gas <- 1.987
    
  term1 <- RH025*(temperature/298)

  arg_term2 <- (HA/R_gas)*(1/298 - 1/temperature)
  term2 <- exp(arg_term2)

  arg_term3 <- (HH/R_gas)*(1/TH - 1/temperature)
  term3 <- exp(arg_term3)

  d_rate_d_temp <-
    term1*term2/(temperature*(1 + term3)) +
    (term1*term2/term3**2)*(((HA/R_gas)/temperature**(2))*(1 + term3) -
                              ((HH/R_gas)/temperature**(2))*term3)
    
  return(d_rate_d_temp)

}
