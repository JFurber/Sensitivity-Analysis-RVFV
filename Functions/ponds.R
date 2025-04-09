#-----------------------------------PONDS---------------------------------------

# From Gianni's file:
# general functions: Pond dynamics based on rainfall data. It is based on the differential equation
# dS_p/dt= a*csi(t)-b*S_p where csi(t) is the rainfall data
# For the rainfall case the area of the pond is  solution of the equation:
# S_p(t)=exp(-b *t)*(S_p(0)+int( a*csi(t)*exp(b *t)))
# with periodic boundary conditions, i.e. S_p(t=0)=S_p(t=1 year) and by imposing the average size equal to a value of our choice


S_p <- function(times) {
  
  # index_time <- which.min(abs(database_thisloc$daynum - times))
  # S_p <- database_thisloc$wbarea[index_time]
  S_p <- constantwb
  
  return(S_p)

}


S_p_periodic<-function(times){
  
  # if (pond_cond=="periodic") {
    S_p<-A_S_p*cos(omega_S_p*times+phi_S_p)+C_S_p
  # }
  # else if (pond_cond=="periodic2") {
  #   S_p<-A_S_p*cos(2*omega_S_p*times+phi_S_p)+C_S_p
  # }
  # else if (pond_cond=="periodic_trend") {
  #   
  #   S_p1<-0
  #   #        time_half<-n_periods*period/2
  #   time_half<-2*period
  #   
  #   if((times>=0) && (times<time_half))
  #   {
  #     
  #     S_p1<- -lambda*times+2*lambda*time_half
  #     
  #   }
  #   else  if ((times>=time_half) &&  (times<2*time_half)){
  #     S_p1<- lambda*times
  #     
  #   }
  #   else if ((times>=2*time_half) &&  (times<3*time_half)){
  #     
  #     
  #     S_p1<- -lambda*times+4*lambda*time_half
  #     
  #   }
  #   else  if ((times>=3*time_half) &&  (times<4*time_half)){
  #     S_p1<- lambda*times-2*lambda*time_half
  #     
  #   }
  #   else if ((times>=4*time_half) &&  (times<5*time_half)){
  #     
  #     S_p1<- -lambda*times+6*lambda*time_half
  #     
  #   }
  #   else  if ((times>=5*time_half) &&  (times<=6*time_half)){
  #     S_p1<- lambda*times-4*lambda*time_half
  #     
  #   }
  #   
  #   
  #   S_p<-(A_S_p*cos(omega_S_p*times+phi_S_p)+C_S_p)+S_p1
  #   
  # }
  
  return(S_p)
  
}



D_S_p_dt <- function(times) {
  
  D_S_p_dt <- 0

  # index_time <- which.min(abs(database_thisloc$daynum - times))
  # 
  # if (index_time == 1) {
  #   D_S_p_dt <- 0
  # }
  # 
  # else {
  #   D_S_p_dt <- (database_thisloc$wbarea[index_time] -
  #                  database_thisloc$wbarea[index_time - 1])/delta_t
  # }

  return(D_S_p_dt)

}


D_S_p_dt_periodic<-function(times)
{
  
  # if (pond_cond=="periodic"){
    D_S_p_dt<- -omega_S_p*A_S_p*sin(omega_S_p*times+phi_S_p)
  # }
  # else if (pond_cond=="periodic2"){
  #   D_S_p_dt<- -2*omega_S_p*A_S_p*sin(2*omega_S_p*times+phi_S_p)
  # }
  # else if (pond_cond=="periodic_trend") {
  #   
  #   
  #   S_p1<-0
  #   
  #   time_half<-2*period
  #   
  #   if((times>=0) && (times<time_half))
  #   {
  #     
  #     S_p1<- -lambda
  #   }
  #   else  if ((times>=time_half) &&  (times<2*time_half)){
  #     S_p1<- lambda
  #   }
  #   else if ((times>=2*time_half) &&  (times<3*time_half)){
  #     
  #     
  #     S_p1<- -lambda
  #   }
  #   else  if ((times>=3*time_half) &&  (times<4*time_half)){
  #     S_p1<- lambda
  #   }
  #   else if ((times>=4*time_half) &&  (times<5*time_half)){
  #     
  #     S_p1<- -lambda
  #   }
  #   else  if ((times>=5*time_half) &&  (times<=6*time_half)){
  #     S_p1<- lambda
  #   }
  #   
  #   #  time_half<-n_periods*period/2
  #   
  #   D_S_p_dt<- -omega_S_p*A_S_p*sin(omega_S_p*times+phi_S_p)+S_p1
  #   
  # }
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


D_S2_p_dt2_periodic<-function(times)
{

  if (pond_cond=="periodic"){
    D_S2_p_dt2<- (omega_S_p)^2*A_S_p*cos(omega_S_p*times+phi_S_p)
    
  }
  else if (pond_cond=="periodic2"){
    D_S2_p_dt2<- (2*omega_S_p)^2*A_S_p*cos(2*omega_S_p*times+phi_S_p)
    
  }
  else if (pond_cond=="periodic_trend") {
    
    print("warning, no explicit function yet")
    
  }

  return(D_S2_p_dt2)
  
}