#-----------------------------------PONDS---------------------------------------

# From Gianni's file:
# general functions: Pond dynamics based on rainfall data. It is based on the differential equation
# dS_p/dt= a*csi(t)-b*S_p where csi(t) is the rainfall data
# For the rainfall case the area of the pond is  solution of the equation:
# S_p(t)=exp(-b *t)*(S_p(0)+int( a*csi(t)*exp(b *t)))
# with periodic boundary conditions, i.e. S_p(t=0)=S_p(t=1 year) and by imposing the average size equal to a value of our choice


S_p <- function(times) {
  
  index_time <- which.min(abs(database_thisloc$daynum - times))
  S_p <- database_thisloc$wbarea[index_time]
  
  return(S_p)

}

D_S_p_dt <- function(times) {

  index_time <- which.min(abs(database_thisloc$daynum - times))

  if (index_time == 1) {
    D_S_p_dt <- 0
  }

  else {
    D_S_p_dt <- (database_thisloc$wbarea[index_time] -
                   database_thisloc$wbarea[index_time - 1])/delta_t
  }

  return(D_S_p_dt)

}

D_S2_p_dt2 <- function(times) {

  index_time <- which.min(abs(database_thisloc$daynum - times))

  if (index_time == 1) {
    D_S2_p_dt2 <- 0
  }

  if (index_time != 1) {
    D_S_p_dt_int <- (database_thisloc$wbarea[index_time] -
                       database_thisloc$wbarea[index_time - 1])/delta_t
    D_S2_p_dt2 <- (D_S_p_dt_int[index_time] -
                     D_S_p_dt_int[index_time - 1])/delta_t
  }

  return(D_S2_p_dt2)

}