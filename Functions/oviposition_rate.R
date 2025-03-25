#---------------------------OVIPOSITION RATES-----------------------------------



oviposition_rate_Culex <- function(t, Pars) {
  
  with(as.list(c(Pars)), { # Bring variables from Pars into function
    
    b_C <- lambda_C # Eggs laid per batch
    
    #carrying capacity
    K_C <- rho_max_C*kappa_C*S_p(t) # Carrying capacity: max eggs in pond # eq(5) GL2018 supplementary
    # Max egg density x area of pond where eggs laid x surface area of pond
    # d_K_C_dt <- rho_max_C*kappa_C*D_S_p_dt(t)
        
    oviposition_rate_Culex <- N_ponds*S_p(t)/(A_C*t_dep)
    # Likelihood of finding pond x egg deposition rate
    # d_oviposition_rate_Culex_dt <- N_ponds*D_S_p_dt(t)/(A_C*t_dep)
       
    # return(list(oviposition_rate_Culex, b_C, K_C, 
                # d_oviposition_rate_Culex_dt, d_K_C_dt))
    return(list(oviposition_rate_Culex, b_C, K_C))
        
    }
    
  )
  
}

oviposition_rate_Aedes <- function(t, Pars_A) {
  
  with(as.list(c(Pars_A)), { # Bring variables from Pars_A to function
    
    b_A <- lambda_A # Eggs laid per batch
    
    K_A <- rho_max_A*kappa_A*S_p(t) # Carrying capacity: max eggs in pond
    # Max egg density x area of pond where eggs laid x surface area of pond
    # d_K_A_dt <- rho_max_A*kappa_A*D_S_p_dt(t)
    
    oviposition_rate_Aedes <- N_ponds*S_p(t)/(A_A*t_dep)
    # Likelihood of finding pond x egg deposition rate
    # d_oviposition_rate_Aedes_dt <- N_ponds*D_S_p_dt(t)/(A_A*t_dep)
    
    # If rate of change of pond area is negative, then eggs won't be submerged
    # If eggs are not submerged they will not hatch
    
    tau_A <- pmax(D_S_p_dt(t)/S_p(t), 0) # Hatching rate
    # d_tau_dS_long <- pmax(-(D_S_p_dt(t)**2/(S_p(t))**2) + D_S2_p_dt2(t)/S_p(t),
                          # 0)
    
    # return(list(oviposition_rate_Aedes, b_A, K_A, tau_A, 
    #             d_oviposition_rate_Aedes_dt, d_K_A_dt, d_tau_dS_long))
    return(list(oviposition_rate_Aedes, b_A, K_A, tau_A))
    
    }
  
  )
  
}