# -----------------------------BITING RATE--------------------------------------

# Biting rate aka number of gonotrophic cycles per day (denoted as theta in the paper)

biting_rate <- function(temperature) {
    
  a_Aedes <- 0.0173*(temperature - 273.16 - 9.60) #GL2018 eq20
  a_Culex <- 0.0173*(temperature - 273.16 - 9.60) #GL2018 eq12
    
  d_a_Aedes_d_temp <- 0.0173 # derivative of above with respect to Temperature
  d_a_Culex_d_temp <- 0.0173 # derivative of above with respect to Temperature
    
  return(c(a_Culex, a_Aedes, d_a_Culex_d_temp, d_a_Aedes_d_temp))

}

# Extrinsic incubation period aka time from exposed to infectious

EIP <- function(temperature) {
  
  # epsilon_C <- 1/(-0.1038 + 0.0071*(temperature))
  # epsilon_A <- 1/(-0.1038 + 0.0071*(temperature))

  epsilon_C <- 1/(-0.1038 + 0.0071*(temperature - 273.15)) #GL2018 eq30 - adjusted by Jess to match eq30
  epsilon_A <- 1/(-0.1038 + 0.0071*(temperature - 273.15)) #GL2018 eq30 - adjusted by Jess to match eq30
    
  d_epsilon_C_d_temp <- -0.0071/(-0.1038 + 0.0071*(temperature - 273.15))**2 # derivative of above with respect to Temperature
  d_epsilon_A_d_temp <- -0.0071/(-0.1038 + 0.0071*(temperature - 273.15))**2 # derivative of above with respect to Temperature
  
  return(c(epsilon_C, epsilon_A, d_epsilon_C_d_temp, d_epsilon_A_d_temp))

}