
#------------------------------TEMPERATURE--------------------------------------

temperature <- function(times) {
  
  index_time <- which.min(abs(database_thisloc$daynum - times))
    
  temperature <- database_thisloc$tmean[index_time] + 273.16

  return(temperature)

}

D_T_a_dt <- function(times) {
  
  index_time <- which.min(abs(database_thisloc$daynum - times))
  
  if (index_time == 1) {
    D_T_a_dt <- 0
  }
  
  if (index_time != 1) {
    D_T_a_dt <- (database_thisloc$tmean[index_time]- 
                   database_thisloc$tmean[index_time - 1])/delta_t
  }
        
  return(D_T_a_dt)
    
}

D_T2_a_dt2 <- function(times) {
  
  print("Warning, no explicit function yet")
        
  return(D_T2_a_dt2)
        
}
