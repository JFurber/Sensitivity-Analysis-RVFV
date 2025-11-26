#--------------------------------MORTALITY--------------------------------------

# Information found in Table S3 GL2018 
# input is temperature in Kelvin K

mortality <- function(temperature) {
  
  temperature_C <- temperature - 273.15
  # temperature_K <- temperature # measured in Kelvins
  
  # Culex eggs - in groups of <10, 10-13, 13-19, 19-30, 30-36, >36 degrees C:
  
  # if (temperature_K < 10 + 273.15) {
  #   mu_C_O <- 0.95 # <10 degrees C
  # } else if (temperature_K >= 10 + 273.15 & temperature_K < 13 + 273.15) {
  #   mu_C_O <- 1 - 0.97 # 10-13 degrees C
  # } else if (temperature_K >= 13 + 273.15 & temperature_K <= 19 + 273.15) {
  #   mu_C_O <- 1 - 54.259*exp(-0.3114*(temperature_K - 273.15)) # 13-19 degrees C
  # } else if (temperature_K > 19 + 273.15 & temperature_K <= 30 + 273.15) {
  #   mu_C_O <- 1 - 0.22 # 
  # } else if (temperature_K > 30 + 273.15 & temperature_K <= 36 + 273.15) {
  #   mu_C_O <- 1 - (0.0876*(temperature_K - 273.15) - 2.3577)
  # } else if (temperature_K > 36 + 273.15) {
  #   mu_C_O <- 0.95
  # }
  
  if (temperature <= 286.15) {
    mu_C_O <- 1 - 0.97 #temperature less than 13 degrees C
  } else if (temperature > 286.15 & temperature <= 292.15) { #temperature between 13 and 19 degrees C
    mu_C_O <- 1 - 54.259*exp(-0.3114*(temperature - 273.15)) 
  } else if (temperature > 292.15 & temperature <= 303.15) { # temperature between 19 and 30 degrees C
    mu_C_O <-1 - 0.22 
  } else if (temperature > 303.15) { #greater than 30 degrees C
    mu_C_O <- 1 - (0.0876*(temperature - 273.15) - 2.3577)
  }
  
  # Culex larvae: Based on reference (14) in GL2018 supplementary
  
  mu_C_L <- 37.9317808331 - 0.2573339304*(temperature) +  0.0004364566*((temperature)**2)
  
  # Culex pupae: Based on reference (14) in GL2018 supplementary
    
  mu_C_P <- 80.3113158804 - 0.5439116495*(temperature) + 0.0009210259*((temperature)**2)
  
  # Culex adult:
    
  mu_C_E <- 0.83  # additional pupal mortality associated with emergence of the adult (6) - symbol delta in the paper
  mu_C <- 0.16 # Daily adult mortality for Culex - based on survivorship of (25)
  
  # @Gianni - Q5 - Does this Aedes egg mortality (mu_A_O) is not dependent on temperature? - no - not according to your supplementary:
  # Aedes egg - immature eggs - crudely estimated as 1/4 years
  
  mu_A_Oi <- 0.0003913894
  
  # Aedes egg - mature eggs - reference (6) GL2018 supplementary
  
  mu_A_Om <- 0.011
  
  # Aedes larvae:
    
  mu_A_L <- 50.1205 - 0.3393650263*temperature + 0.0005747698*(temperature)**2   
  
  # Aedes pupae:
  
  mu_A_P <- 3.524873 - (2.394308e-02)*temperature + 
    (4.066735e-05)*(temperature)**2 
  
  # @Gianni - Q6 - Unsure on the following two (mu_A_E and mu_A) - What are they?:
  
  #Aedes adult:
    
  mu_A_E <- 0.83 # additional pupal mortality, same as mu_C_E - symbol delta in the paper
  mu_A <- 1/(25.8 - 0.45*(temperature - 273.16))  
  
  # Return all mortality rates:
  
  return(c(mu_C_O, mu_C_L, mu_C_P, mu_C_E, mu_C, mu_A_Oi, mu_A_Om,
           mu_A_L, mu_A_P, mu_A_E, mu_A))
  
}
